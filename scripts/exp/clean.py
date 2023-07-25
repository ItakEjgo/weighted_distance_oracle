import os
import sys


if __name__ == '__main__':
    argc = len(sys.argv)
    if argc < 2 or argc > 3:
        print('[Invalid], Usage: ./clean.py [data_dir] (with default 100 queries) or ./clean.py [data_dir] [query_num]')
    else:
        if (argc == 2):
            query_num = 100 # default query num
        else:
            query_num = int(sys.argv[2]) # user defined query num

        dis = []
        time = []
        for i in range(query_num):
            dis.append(float('nan'))
            time.append(float('nan'))
        index_time = float('nan')
        index_size = float('nan')
        breakdown_construction = float('nan')
        breakdown_dijkstra = float('nan')

        data_dir = sys.argv[1]
        with open(data_dir) as f:
            for line in f:
                if line.find('Query results begin') >= 0:
                    # for x in range(3):
                    for j in range(query_num):
                        cur_line = f.readline()
                        dis[j] = (float(cur_line.split(' ')[0]))
                        time[j] = (float(cur_line.split(' ')[1]))
                elif line.find('Index Time') >= 0:
                    index_time = float(line.split(' ')[3])
                elif line.find('Index memory usage') >= 0:
                    index_size = float(line.split(' ')[3])
                elif line.find('Breakdown_construction') >= 0:
                    breakdown_construction = float(line.split(' ')[1])
                elif line.find('Breakdown_dijkstra') >= 0:
                    breakdown_dijkstra = float(line.split(' ')[1])

        # inner_time = 0.0
        # inter_time = 0.0
        mixed_time = 0.0
        for i in range(query_num):
            # inner_time += time[i]
            # inter_time += time[query_num + i]
            # mixed_time += time[2 * query_num + i]
            mixed_time += time[i]
        mixed_time /= query_num
        output_dir = data_dir[:-4] + '.cln'
        print('output file is written as: ' + output_dir)
        f = open(output_dir, 'w')
        print('index_time(ms): ', index_time, file=f)
        print('index_size(MB): ', index_size, file=f)
        # print('average_inner_time(micro-s): ', inner_time, file=f)
        # print('average_inter_time(micro-s): ', inter_time, file=f)
        print('average_mixed_time(micro-s): ', mixed_time, file=f)
        print('query_construction(ms): ', breakdown_construction, file=f)
        print('query_dijkstra(ms): ', breakdown_dijkstra, file=f)
        # for i in range(query_num * 3):
        for i in range(query_num):
            print(dis[i], file=f)
        print('clean [' + data_dir + '] finished')
        