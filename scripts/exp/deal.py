import os
import sys

out_dir = '../exp'

if __name__ == '__main__':
    q_num = int(sys.argv[4])
    print('target file:' + sys.argv[1] + '-' + sys.argv[3] + '.log' + " \npartition size:" + sys.argv[2] + " \nflag:" + sys.argv[3])
    partition_num = int(sys.argv[2])
    dis = [0.0] * q_num * 3
    time = [0.0] * q_num * 3
    index_time = float('nan')
    index_size = float('nan')
    breakdown_construction = float('nan')
    breakdown_dijkstra = float('nan')
    for i in range(partition_num):
        cur_file = sys.argv[1] + '-part' + str(i) + '.log'
        print('deal ' + cur_file + '...')
        with open(out_dir + '/' + cur_file, 'r') as f:
            for line in f:
                if line.find('Query results begin') >= 0:
                    for x in range(3):
                        for j in range(int(q_num / partition_num)):
                            cur_line = f.readline()
                            id = int(x * q_num) + int(q_num / partition_num * i) + j
                            dis[id] = float(cur_line.split(' ')[0])
                            time[id] = float(cur_line.split(' ')[1])
                elif line.find('Index Time') >= 0:
                    index_time = float(line.split(' ')[3])
                elif line.find('Index memory usage') >= 0:
                    index_size = float(line.split(' ')[3])
                elif line.find('Breakdown_construction') >= 0:
                    breakdown_construction = float(line.split(' ')[1])
                elif line.find('Breakdown_dijkstra') >= 0:
                    breakdown_dijkstra = float(line.split(' ')[1])


    inner_time = 0.0
    inter_time = 0.0
    mixed_time = 0.0
    for i in range(q_num):
        inner_time += time[i]
        inter_time += time[q_num + i]
        mixed_time += time[2 * q_num + i]
    print('output file is written to:' + out_dir + '/' + sys.argv[1] + '-' + sys.argv[3] + '.log')
    f = open(out_dir + '/' + sys.argv[1] + '-' + sys.argv[3] + '.log', 'w')
    print('index time(ms): ', index_time, file=f)
    print('index size(MB): ', index_size, file=f)
    print('inner_time(ms): ', inner_time, file=f)
    print('inter_time(ms): ', inter_time, file=f)
    print('mixed_time(ms): ', mixed_time, file=f)
    print('query_construction(ms): ', breakdown_construction, file=f)
    print('query_dijkstra(ms): ', breakdown_dijkstra, file=f)
    print('clean files...')
    os.system('rm ' + out_dir + '/' + sys.argv[1] + '-part*')
    for i in range(q_num * 3):
        print(dis[i], file=f)
    
        
    

