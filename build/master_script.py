import os
import sys

def deal(config_dir):
    variable_dict = {}
    with open(config_dir) as f:
        for line in f:
            key = line.split('=')[0].strip()
            val = line.split('=')[1].strip()
            variable_dict[key] = val
    return variable_dict

def generate_query(variable_dict, gridnum, query_num, terrain_type):
    input_dir = variable_dict['datasets_dir'] + 'small/'
    query_dir = variable_dict['query_dir']
    run_script = variable_dict['scripts_dir'] + 'exp/generate_query.sh'
    cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + gridnum + ' ' + query_num + ' ' + terrain_type
    # print(cmd)
    os.system(cmd)
    print(terrain_type + ' query generate finished')

def run_default_small(variable_dict, gridnum, query_num, algorithm):
    input_dir = variable_dict['datasets_dir'] + 'small/'
    output_dir = variable_dict['output_dir'] + 'default/'
    os.system('mkdir -p ' + output_dir)
    query_dir = variable_dict['query_dir']  + ' '
    run_script = variable_dict['scripts_dir'] + 'exp/exp_default.sh'
    cleaner = variable_dict['scripts_dir'] + 'exp/clean.py'
    # print(run_script_dir)
    cmd = 'bash ' + run_script + ' ' + input_dir + ' ' + query_dir + ' ' + output_dir + ' ' + gridnum + ' ' + query_num + ' ' + algorithm + ' ' + cleaner
    os.system(cmd)
    # print(cmd)

if __name__ == '__main__':
    argc = len(sys.argv)
    if argc != 2:
        print('[Invalide], Usage: python3 master_script.py [config path]')
    else:
        config_dir = sys.argv[1]
        variable_dict = deal(config_dir)
        # generate_query(variable_dict, '16', '100', 'default')
        run_default_small(variable_dict, '16', '100', 'EAR')