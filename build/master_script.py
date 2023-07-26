import os
import sys

def deal(config_dir):
    variable_dict = {}
    with open(config_dir) as f:
        for line in f:
            key = line.split('=')[0].strip()
            val = line.split('=')[1].strip()
            if (val.find(',') >= 0):
                val = list(x.strip() for x in list(val.split(',')))
            variable_dict[key] = val
    # for key in variable_dict:
        # print(key, variable_dict[key])
    return variable_dict

def generate_query(variable_dict, terrain_type):
    tested_dataset = variable_dict['tested_dataset']
    for dir in tested_dataset:
        input_dir = variable_dict['datasets_dir'] + dir + '/'
        query_dir = variable_dict['query_dir']
        run_script = variable_dict['scripts_dir'] + 'exp/generate_query.sh'
        gridnum = '16' if dir == 'small' else '256'
        if terrain_type == 'default':
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' 0 ' + gridnum + terrain_type + ' 0'
        elif terrain_type == 'weighted':
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' 1 ' + gridnum + terrain_type + ' 0'
        else:
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' 0 ' + gridnum + terrain_type + ' 1'
        # print(cmd)
        os.system(cmd)
    print(terrain_type + ' query generate finished')

def read_clean_file(file_path):
    clean_dict = {}
    dis = []
    for i in range(100):
        dis[i].append(float(nan))
    with open(file_path) as f:
        for line in f:
            if line.find('index time') >= 0:
                clean_dict['index_time'] = line.split(':')[1].strip()
            elif line.find('index size') >= 0:
                clean_dict['index_size'] = line.split(':')[1].strip()
            elif line.find('average_mixed_time') >= 0:
                amt = float(line.split(':')[1].strip()) / 1000
                clean_dict['query_time'] = str(amt)
            elif line.find('query_construction') >= 0:
                clean_dict['query_construction'] = line.split(':')[1].strip()
            elif line.find('query_dijkstra') >= 0:
                clean_dict['query_dijkstra'] = line.split(':')[1].strip()
                for i in range(100):
                    cur_line = f.readline()
                    dis[i] = float(cur_line)
                clean_dict['dis'] = dis
    return clean_dict

def read_disgap(file_path):
    dis = []
    for i in range(10):
        dis.append(float(nan))
    with open(file_path) as f:
        for line in f:
            if line.find('Query results begin') >= 0:
                for k in range(10):
                    t_query_time = 0.0
                    for i in range(100):
                        cur_line = f.readline()
                        t_query_time += float(cur_line.split(' ')[1].strip())
                    t_query_time /= 100000  # average time in million-seconds
                    dis[k] = t_query_time
    return dis

                    
# List, List -> float
# return the distance error of two list. Pivot is the second list.
def calc_error(dis_1, dis_2):
    error = 0.0
    for i in range(len(dis_1)):
        error += abs(dis_1[i] - dis_2[i]) / dis_2[i]
    error /= len(dis_1)
    return error

def plot_default(variable_dict):
    clean_dir = variable_dict['output_dir'] + 'default/'
    pos_1 = [3, 10, 17, 24, 31, 38, 45, 51]
    pos_2 = [5, 15, 25, 35, 45, 55, 65, 75]
    plot_data = {}
    for algorithm in variable_dict['algorithms']:
        if (len(algorithm) <= 0):
            continue
        plot_data[algorithm] = {}
        for dataset in variable_dict['dataset_list']:
            break
            file_name = dataset + '-' + algorithm + '-' + 'default' + '.cln'
            file_path = clean_dir + file_name
            clean_dict = read_clean_file(file_path)
            plot_data[algorithm][dataset] = clean_dict
    for dataset in variable_dict['dataset_list']:
        continue
        for algorithm in variable_dict['algorithms']:
            plot_data[algorithm][dataset]['error'] = float(nan)
            if algorithm == 'MMP':
                continue
            relative_error = calc_error(plot_data[algorithm][dataset][dis], plot_data['MMP'][dataset][dis])
            plot_data[algorithm][dataset]['error'] = relative_error
    plot_res_dir = variable_dict['scripts_dir'] + 'figures/data'
    os.system('mkdir -p ' + plot_res_dir)
    os.system('mkdir -p out')
    plot_default_res = plot_res_dir + '/default.res'
    f = open(plot_default_res, 'w')
    for i in range(len(variable_dict['dataset_list'])):
        cur_dataset = variable_dict['dataset_list'][i]
        print(cur_dataset)
        continue
        print(pos_1[i], pos_2[i], 
            plot_data['SE'][dataset]['index_time'],
            plot_data['EAR'][dataset]['index_time'],
            plot_data['SE'][dataset]['index_size'],
            plot_data['EAR'][dataset]['index_size'],
            plot_data['FixedS'][dataset]['query_time'],
            plot_data['UnfixedS'][dataset]['query_time'],
            plot_data['Kalgo'][dataset]['query_time'],
            plot_data['SE'][dataset]['query_time'],
            plot_data['EAR'][dataset]['query_time'],
            plot_data['FixedS'][dataset]['error'],
            plot_data['UfixedS'][dataset]['error'],
            plot_data['Kalgo'][dataset]['error'],
            plot_data['SE'][dataset]['error'],
            plot_data['EAR'][dataset]['error'],
            file=f)
    cmd = 'gnuplot -c ' + variable_dict['scripts_dir'] + 'figures/default2-2row.plot ' + plot_default_res
    print(cmd)
    # os.system(cmd)

def run_default(variable_dict):
    algorithm = variable_dict['algorithms']
    tested_dataset = variable_dict['tested_dataset']
    for dir in tested_dataset:
        if len(dir) == 0:
            continue
        input_dir = variable_dict['datasets_dir'] + dir + '/'
        output_dir = variable_dict['output_dir'] + 'default/'
        os.system('mkdir -p ' + output_dir)
        query_dir = variable_dict['query_dir']  + ' '
        run_script = variable_dict['scripts_dir'] + 'exp/exp_default.sh'
        cleaner = variable_dict['scripts_dir'] + 'exp/clean.py'
        gridnum = '16' if dir == 'small' else '256'
        # print(run_script_dir)
        for algo in algorithm:
            if len(algo) == 0:
                continue
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + output_dir + ' ' + algo + ' ' + gridnum + ' ' + cleaner
            os.system(cmd)
            # print(cmd)
        
def plot_weighted(variable_dict):
    clean_dir = variable_dict['output_dir'] + 'weighted/'
    pos_1 = [3, 10, 17, 24, 31, 38, 45, 51]
    pos_2 = [5, 15, 25, 35, 45, 55, 65, 75]
    plot_data = {}
    for algorithm in variable_dict['algorithms']:
        if len(algorithm) <= 0 or algorithm == 'MMP':
            continue
        plot_data[algorithm] = {}
        for dataset in variable_dict['dataset_list']:
            break
            file_name = dataset + '-' + algorithm + '-weighted' + '.cln'
            file_path = clean_dir + file_name
            clean_dict = read_clean_file(file_path)
            plot_data[algorithm][dataset] = clean_dict
    for dataset in variable_dict['dataset_list']:
        continue
        for algorithm in variable_dict['algorithms']:
            plot_data[algorithm][dataset]['error'] = float(nan)
            if algorithm == 'MMP' or algorithm == 'FixedS':
                continue
            relative_error = calc_error(plot_data[algorithm][dataset][dis], plot_data['FixedS'][dataset][dis])
            plot_data[algorithm][dataset]['error'] = relative_error
    plot_res_dir = variable_dict['scripts_dir'] + 'figures/data'
    os.system('mkdir -p ' + plot_res_dir)
    os.system('mkdir -p out')
    plot_weighted_res = plot_res_dir + '/weighted.res'
    f = open(plot_weighted_res, 'w')
    for i in range(len(variable_dict['dataset_list'])):
        cur_dataset = variable_dict['dataset_list'][i]
        print(cur_dataset)
        continue
        print(pos_1[i], pos_2[i], 
            plot_data['SE'][dataset]['index_time'],
            plot_data['EAR'][dataset]['index_time'],
            plot_data['SE'][dataset]['index_size'],
            plot_data['EAR'][dataset]['index_size'],
            plot_data['FixedS'][dataset]['query_time'],
            plot_data['UnfixedS'][dataset]['query_time'],
            plot_data['Kalgo'][dataset]['query_time'],
            plot_data['SE'][dataset]['query_time'],
            plot_data['EAR'][dataset]['query_time'],
            plot_data['FixedS'][dataset]['error'],
            plot_data['UfixedS'][dataset]['error'],
            plot_data['Kalgo'][dataset]['error'],
            plot_data['SE'][dataset]['error'],
            plot_data['EAR'][dataset]['error'],
            file=f)
    cmd = 'gnuplot -c ' + variable_dict['scripts_dir'] + 'figures/weighted-2row.plot ' + plot_weighted_res
    print(cmd)
    # os.system(cmd)

def run_weighted(variable_dict):
    algorithm = variable_dict['algorithms']
    tested_dataset = variable_dict['tested_dataset']
    for dir in tested_dataset:
        if len(dir) == 0:
            continue
        input_dir = variable_dict['datasets_dir'] + dir + '/'
        output_dir = variable_dict['output_dir'] + 'weighted/'
        os.system('mkdir -p ' + output_dir)
        query_dir = variable_dict['query_dir']  + ' '
        run_script = variable_dict['scripts_dir'] + 'exp/exp_weighted.sh'
        cleaner = variable_dict['scripts_dir'] + 'exp/clean.py'
        gridnum = '16' if dir == 'small' else '256'
        # print(run_script_dir)
        for algo in algorithm:
            if len(algo) == 0 or algo == 'MMP':
                print('[MMP] does not support weighted terrain. Experiment for ', algo, ' is skipped.') 
                continue
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + output_dir + ' ' + algo + ' ' + gridnum + ' ' + cleaner
            os.system(cmd)
            # print(cmd)

def plot_epsilon(variable_dict):
    clean_dir = variable_dict['output_dir'] + 'epsilon/'
    pos_1 = [0.05, 0.10, 0.15, 0.20, 0.25]
    plot_data = {}
    dataset = 'BigMount.off'
    for algorithm in variable_dict['algorithms']:
        if len(algorithm) <= 0:
            continue
        plot_data[algorithm] = {}
        for eps in variable_dict['epsilon_val']:
            file_name = dataset + '-' + algorithm + '-eps' + '.cln'
            file_path = clean_dir + file_name
            # use default results for MMP and FixedS
            if algorithm == 'FixedS' or algorithm == 'MMP':
                file_name = dataset + '-' + algorithm + '-default.cln'
                file_path = variable_dict['output_dir'] + 'default/' + file_name
            clean_dict = read_clean_file(file_path)
            plot_data[algorithm][eps] = clean_dict

    for eps in variable_dict['epsilon_val']:
        for algorithm in variable_dict['algorithms']:
            if len(algorithm) < 0 or algorithm == 'MMP':
                continue
            plot_data[algorithm][eps]['error'] = float(nan)
            relative_error = calc_error(plot_data[algorithm][eps][dis], plot_data['MMP'][eps][dis])
            plot_data[algorithm][eps]['error'] = relative_error

    plot_res_dir = variable_dict['scripts_dir'] + 'figures/data'
    os.system('mkdir -p ' + plot_res_dir)
    os.system('mkdir -p out')
    plot_epsilon_res = plot_res_dir + '/epsilon.res'
    f = open(plot_epsilon_res, 'w')

    for i in range(len(variable_dict['epsilon_val'])):
        cur_epsilon = variable_dict['epsilon_val'][i]
        print(cur_epsilon)
        continue
        print(pos_1[i], 
            plot_data['SE'][cur_epsilon]['index_time'],
            plot_data['SE'][cur_epsilon]['index_size'],
            plot_data['EAR'][cur_epsilon]['index_time'],
            plot_data['EAR'][cur_epsilon]['index_size'],
            plot_data['FixedS'][cur_epsilon]['query_time'],
            plot_data['UnfixedS'][cur_epsilon]['query_time'],
            plot_data['Kalgo'][cur_epsilon]['query_time'],
            plot_data['SE'][cur_epsilon]['query_time'],
            plot_data['EAR'][cur_epsilon]['query_time'],
            plot_data['FixedS'][cur_epsilon]['error'],
            plot_data['UfixedS'][cur_epsilon]['error'],
            plot_data['Kalgo'][cur_epsilon]['error'],
            plot_data['SE'][cur_epsilon]['error'],
            plot_data['EAR'][cur_epsilon]['error'],
            file=f)
    cmd = 'gnuplot -c ' + variable_dict['scripts_dir'] + 'figures/eps-2row.plot ' + plot_epsilon_res
    print(cmd)
    # os.system(cmd)


def run_epsilon(variable_dict):
    algorithm = variable_dict['algorithms']
    tested_dataset = variable_dict['tested_dataset']
    input_dir = variable_dict['datasets_dir'] + 'epsilon/'
    output_dir = variable_dict['output_dir'] + 'epsilon/'
    os.system('mkdir -p ' + output_dir)
    query_dir = variable_dict['query_dir']  + ' '
    run_script = variable_dict['scripts_dir'] + 'exp/exp_epsilon.sh'
    cleaner = variable_dict['scripts_dir'] + 'exp/clean.py'
    # print(run_script_dir)
    for algo in algorithm:
        if algo == 'FixedS' or algo == 'MMP':
            print('Grid number only influences [UnfixedS, KAlgo, SE-Oracle and EAR-Oracle]. Experiment for ', algo, ' is skipped.')
            continue
        for i in range(len(variable_dict['epsilon_val'])):
            eps = variable_dict['epsilon_val'][i]
            flag = variable_dict['epsilon_flag'][i]
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + output_dir + ' ' + algo + ' ' + eps + ' ' + flag + ' ' + cleaner
            os.system(cmd)
            # print(cmd)


def plot_gridnum(variable_dict):
    clean_dir = variable_dict['output_dir'] + 'girdnum/'
    pos_1 = [0.05, 0.10, 0.15, 0.20, 0.25]
    plot_data = {}
    dataset = 'Simply_GunnisonForest.off'
    for algorithm in variable_dict['algorithms']:
        if len(algorithm) <= 0:
            continue
        plot_data[algorithm] = {}
        for girdnum in variable_dict['gridnum_val']:
            file_name = dataset + '-' + algorithm + '-gridnum' + '.cln'
            file_path = clean_dir + file_name
            # use default results for all other algorithms. 
            if algorithm != 'EAR':
                file_name = dataset + '-' + algorithm + '-default.cln'
                file_path = variable_dict['output_dir'] + 'default/' + file_name
            clean_dict = read_clean_file(file_path)
            plot_data[algorithm][girdnum] = clean_dict

    for gridnum in variable_dict['gridnum_val']:
        for algorithm in variable_dict['algorithms']:
            if len(algorithm) < 0 or algorithm == 'MMP':
                continue
            plot_data[algorithm][gridnum]['error'] = float(nan)
            relative_error = calc_error(plot_data[algorithm][gridnum][dis], plot_data['MMP'][gridnum][dis])
            plot_data[algorithm][gridnum]['error'] = relative_error

    plot_res_dir = variable_dict['scripts_dir'] + 'figures/data'
    os.system('mkdir -p ' + plot_res_dir)
    os.system('mkdir -p out')
    plot_gridnum_res = plot_res_dir + '/gridnum.res'
    f = open(plot_gridnum_res, 'w')

    for i in range(len(variable_dict['gridnum_val'])):
        cur_gridnum = variable_dict['gridnum_val'][i]
        print(cur_gridnum)
        continue
        print(pos_1[i], 
            plot_data['EAR'][cur_gridnum]['index_time'],
            plot_data['EAR'][cur_gridnum]['index_size'],
            plot_data['EAR'][cur_gridnum]['query_time'],
            plot_data['EAR'][cur_gridnum]['error'],
            plot_data['FixedS'][cur_gridnum]['query_time'],
            plot_data['UnfixedS'][cur_gridnum]['query_time'],
            plot_data['Kalgo'][cur_gridnum]['query_time'],
            plot_data['FixedS'][cur_gridnum]['error'],
            plot_data['UfixedS'][cur_gridnum]['error'],
            plot_data['Kalgo'][cur_gridnum]['error'],
            file=f)
    cmd = 'gnuplot -c ' + variable_dict['scripts_dir'] + 'figures/grid-2row.plot ' + plot_gridnum_res
    print(cmd)
    # os.system(cmd)

def run_gridnum(variable_dict):
    algorithm = variable_dict['algorithms']
    tested_dataset = variable_dict['tested_dataset']
    input_dir = variable_dict['datasets_dir'] + 'gridnum/'
    output_dir = variable_dict['output_dir'] + 'gridnum/'
    os.system('mkdir -p ' + output_dir)
    query_dir = variable_dict['query_dir']  + ' '
    run_script = variable_dict['scripts_dir'] + 'exp/exp_gridnum.sh'
    cleaner = variable_dict['scripts_dir'] + 'exp/clean.py'
    # print(run_script_dir)
    for algo in algorithm:
        if algo != 'EAR':
            print('Grid number only influences [EAR-Oracle]. Experiment for ', algo, ' is skipped.')
            continue
        for i in range(len(variable_dict['gridnum_val'])):
            gridnum = variable_dict['gridnum_val'][i]
            flag = variable_dict['gridnum_flag'][i]
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + output_dir + ' ' + algo + ' ' + gridnum + ' ' + flag + ' ' + cleaner
            os.system(cmd)
            # print(cmd)

def plot_spnum(variable_dict):
    clean_dir = variable_dict['output_dir'] + 'spnum/'
    pos_1 = [0.05, 0.10, 0.15, 0.20, 0.25]
    plot_data = {}
    dataset = 'BigMount.off'
    for algorithm in variable_dict['algorithms']:
        if len(algorithm) <= 0:
            continue
        plot_data[algorithm] = {}
        for spnum in variable_dict['spnum_val']:
            file_name = dataset + '-' + algorithm + '-spnum' + '.cln'
            file_path = clean_dir + file_name
            # use default results for all other algorithms. 
            if algorithm == 'UnfixedS' or algorithm == 'KAlgo' or algorithm == 'MMP':
                file_name = dataset + '-' + algorithm + '-default.cln'
                file_path = variable_dict['output_dir'] + 'default/' + file_name
            clean_dict = read_clean_file(file_path)
            plot_data[algorithm][spnum] = clean_dict

    for spnum in variable_dict['spnum_val']:
        for algorithm in variable_dict['algorithms']:
            if len(algorithm) < 0 or algorithm == 'MMP':
                continue
            plot_data[algorithm][spnum]['error'] = float(nan)
            relative_error = calc_error(plot_data[algorithm][spnum][dis], plot_data['MMP'][spnum][dis])
            plot_data[algorithm][spnum]['error'] = relative_error

    plot_res_dir = variable_dict['scripts_dir'] + 'figures/data'
    os.system('mkdir -p ' + plot_res_dir)
    os.system('mkdir -p out')
    plot_spnum_res = plot_res_dir + '/spnum.res'
    f = open(plot_spnum_res, 'w')

    for i in range(len(variable_dict['spnum_val'])):
        cur_spnum = variable_dict['spnum_val'][i]
        print(cur_spnum)
        continue
        print(pos_1[i], 
            plot_data['SE'][cur_spnum]['index_time'],
            plot_data['SE'][cur_spnum]['index_size'],
            plot_data['EAR'][cur_spnum]['index_time'],
            plot_data['EAR'][cur_spnum]['index_size'],
            plot_data['FixedS'][cur_spnum]['query_time'],
            plot_data['UnfixedS'][cur_spnum]['query_time'],
            plot_data['Kalgo'][cur_spnum]['query_time'],
            plot_data['SE'][cur_spnum]['query_time'],
            plot_data['EAR'][cur_spnum]['query_time'],
            plot_data['FixedS'][cur_spnum]['error'],
            plot_data['UfixedS'][cur_spnum]['error'],
            plot_data['Kalgo'][cur_spnum]['error'],
            plot_data['SE'][cur_spnum]['error'],
            plot_data['EAR'][cur_spnum]['error'],
            file=f)
    cmd = 'gnuplot -c ' + variable_dict['scripts_dir'] + 'figures/spnum-2row.plot ' + plot_spnum_res
    print(cmd)
    # os.system(cmd)

def run_spnum(variable_dict):
    algorithm = variable_dict['algorithms']
    tested_dataset = variable_dict['tested_dataset']
    input_dir = variable_dict['datasets_dir'] + 'spnum/'
    output_dir = variable_dict['output_dir'] + 'spnum/'
    os.system('mkdir -p ' + output_dir)
    query_dir = variable_dict['query_dir']  + ' '
    run_script = variable_dict['scripts_dir'] + 'exp/exp_spnum.sh'
    cleaner = variable_dict['scripts_dir'] + 'exp/clean.py'
    # print(run_script_dir)
    for algo in algorithm:
        if algo == 'UnfixedS' or algo == 'KAlgo' or algo == 'MMP':
            print('Steiner points number only influences [FixedS, SE-Oracle and EAR-Oracle]. Experiment for ', algo, ' is skipped.')
            continue
        for i in range(len(variable_dict['spnum_val'])):
            spnum = variable_dict['spnum_val'][i]
            flag = variable_dict['spnum_flag'][i]
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + output_dir + ' ' + algo + ' ' + spnum + ' ' + flag + ' ' + cleaner
            os.system(cmd)
            # print(cmd)

def plot_disgap_breakdown(variable_dict):
    disgap_dir = variable_dict['output_dir'] + 'disgap/'
    pos_1 = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    pos_2 = [5, 15, 25, 35, 45, 55, 65, 75]
    plot_data = {}
    dataset = 'Simply_GunnisonForest.off'
    for algorithm in variable_dict['algorithms']:
        if len(algorithm) <= 0 or algorithm == 'MMP':
            continue
        file_name = dataset + '-' + algorithm + '-disgap' + '.log'
        file_path = clean_dir + file_name
        query_time_list = read_disgap(file_path)
        plot_data[algorithm] = query_time_list

    plot_res_dir = variable_dict['scripts_dir'] + 'figures/data/'
    os.system('mkdir -p ' + plot_res_dir)
    os.system('mkdir -p out')
    plot_disgap_res = plot_res_dir + 'disgap.res'
    f = open(plot_disgap_res, 'w')

    for i in range(10):
        print(pos_1[i], 
            plot_data['EAR'][i],
            plot_data['FixedS'][i],
            plot_data['UnfixedS'][i],
            plot_data['Kalgo'][i],
            plot_data['SE'][i],
            file=f)

    breakdown_dir = variable_dict['output_dir'] + 'default/'
    breakdown_data = {}
    for dataset in variable_dict['dataset_list']:
        file_name = dataset + '-EAR-default.cln' 
        file_path = breakdown_dir + file_name
        clean_dict = read_clean_file(file_path)
        breakdown_data[dataset] = clean_dict
    
    plot_breakdown_res = plot_res_dir + 'breakdown.res'
    f2 = open(plot_breakdown_res, 'w')
    for i in range(len(variable_dict['dataset_list'])):
        cur_dataset = variable_dict['dataset_list'][i]
        print(pos_2[i],
            breakdown_data[cur_dataset]['query_construction'],
            breakdown_data[cur_dataset]['query_dijkstra'],
            file=f2)

    cmd = 'gnuplot -c ' + variable_dict['scripts_dir'] + 'figures/breakdown-disgap.plot ' + plot_breakdown_res + ' ' + plot_disgap_res
    print(cmd)
    # os.system(cmd)

def run_disgap(variable_dict):
    algorithm = variable_dict['algorithms']
    tested_dataset = variable_dict['tested_dataset']
    input_dir = variable_dict['datasets_dir'] + 'disgap/'
    output_dir = variable_dict['output_dir'] + 'disgap/'
    os.system('mkdir -p ' + output_dir)
    query_dir = variable_dict['query_dir']  + ' '
    run_script = variable_dict['scripts_dir'] + 'exp/exp_disgap.sh'
    cleaner = variable_dict['scripts_dir'] + 'exp/clean.py'
    gridnum = '256'
    # print(run_script_dir)
    for algo in algorithm:
        if len(algo) == 0:
            continue
        cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + output_dir + ' ' + algo + ' ' + gridnum + ' ' + cleaner
        os.system(cmd)
        # print(cmd)

def plot_scalability(variable_dict):
    disgap_dir = variable_dict['output_dir'] + 'scalability/'
    pos_1 = [0.05, 0.10, 0.15, 0.20, 0.25]
    plot_data = {}
    for algorithm in variable_dict['algorithms']:
        if len(algorithm) <= 0:
            continue
        plot_data[algorithm] = {}
        for dataset in variable_dict['scalability_list']:
            file_name = dataset + '-' + algorithm + '-scalability' + '.cln'
            file_path = clean_dir + file_name
            clean_dict = read_clean_file(file_path)
            plot_data[algorithm][dataset] = clean_dict

    for dataset in variable_dict['scalability_list']:
        for algorithm in variable_dict['algorithms']:
            plot_data[algorithm][dataset]['error'] = float(nan)
            if algorithm == 'MMP':
                continue
            relative_error = calc_error(plot_data[algorithm][dataset][dis], plot_data['MMP'][dataset][dis])
            plot_data[algorithm][dataset]['error'] = relative_error

    plot_res_dir = variable_dict['scripts_dir'] + 'figures/data/'
    os.system('mkdir -p ' + plot_res_dir)
    os.system('mkdir -p out')

    Averagerr_dir = plot_res_dir + 'EP_high_Averagerr.res'
    f1 = open(Averagerr_dir, 'w')

    for i in range(len(variable_dict['scalability_list'])):
        cur_dataset = variable_dict['scalability_list'][i]
        print(pos_1[i],
            plot_data['FixedS'][cur_dataset]['error'],
            plot_data['KAlgo'][cur_dataset]['error'],
            plot_data['EAR'][cur_dataset]['error'],
            plot_data['SE'][cur_dataset]['error'],
            plot_data['UnfixedS'][cur_dataset]['error'],
            file=f1)

    cmd = 'gnuplot -c ' + variable_dict['scripts_dir'] + 'figures/n_effect-2row.plot ' + plot_res_dir + 'EP_high'
    print(cmd)
    # os.system(cmd)

def run_scalability(variable_dict):
    algorithm = variable_dict['algorithms']
    tested_dataset = variable_dict['tested_dataset']
    input_dir = variable_dict['datasets_dir'] + 'scalability/'
    output_dir = variable_dict['output_dir'] + 'scalability/'
    os.system('mkdir -p ' + output_dir)
    query_dir = variable_dict['query_dir']  + ' '
    run_script = variable_dict['scripts_dir'] + 'exp/exp_scalability.sh'
    cleaner = variable_dict['scripts_dir'] + 'exp/clean.py'
    gridnum = '16' if dir == 'small' else '256'
    # print(run_script_dir)
    for algo in algorithm:
        if len(algo) == 0:
            continue
        cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + output_dir + ' ' + algo + ' ' + gridnum + ' ' + cleaner
        os.system(cmd)
        # print(cmd)

if __name__ == '__main__':
    argc = len(sys.argv)
    if argc != 2:
        print('[Invalide], Usage: python3 master_script.py [config path]')
    else:
        config_dir = sys.argv[1]
        variable_dict = deal(config_dir)
        # generate_query(variable_dict, 'default')
        # generate_query(variable_dict, 'weighted')
        # generate_query(variable_dict, 'disgap')
        run_scalability(variable_dict)
        run_disgap(variable_dict)
        run_spnum(variable_dict)
        run_gridnum(variable_dict)
        run_epsilon(variable_dict)
        run_weighted(variable_dict)
        run_default(variable_dict)
        # plot_default(variable_dict)
        # plot_weighted(variable_dict)
        # plot_epsilon(variable_dict)
