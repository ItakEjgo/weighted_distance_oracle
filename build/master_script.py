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
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' 0 ' + gridnum + terrain_type
        else:
            cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' 1 ' + gridnum + terrain_type
        # print(cmd)
        os.system(cmd)
    print(terrain_type + ' query generate finished')

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

if __name__ == '__main__':
    argc = len(sys.argv)
    if argc != 2:
        print('[Invalide], Usage: python3 master_script.py [config path]')
    else:
        config_dir = sys.argv[1]
        variable_dict = deal(config_dir)
        generate_query(variable_dict, 'default')
        generate_query(variable_dict, 'weighted')
        run_default(variable_dict)
        run_weighted(variable_dict)
        run_epsilon(variable_dict)
        run_gridnum(variable_dict)
        run_spnum(variable_dict)
