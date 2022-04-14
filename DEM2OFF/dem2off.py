from re import S
import sys
import os
import pymeshlab
import random


def DEMprojection(data_path, data_name):
    projection_cmd = 'gdalwarp -t_srs EPSG:3857 ' + data_path + ' tifs/' + data_name + '.tif'
    # print(projection_cmd)
    os.system(projection_cmd)
    print('GDAL project to EPSG:3857 finished.')


def DEM2PLY(data_name):
    in_data_path = 'tifs/' + data_name + '.tif'
    dem2ply_cmd = './dem2mesh -verbose --aggressiveness 0 -inputFile ' + in_data_path + ' -outputFile ' + 'plys/' + data_name + '.ply'
    os.system(dem2ply_cmd)
    print('DEM2PLY finished.')
    # print(dem2ply_cmd)


def PLY2OBJ(data_name):
    out_data_path = 'plys/' + data_name + '.ply'
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(out_data_path)
    ms.save_current_mesh(data_name+'.obj')
    # clean_cmd = './manifold --input --depth 8 ' + data_name + '.obj --output ' + data_name + '.obj'
    # os.system(clean_cmd)



def DEM2OFF(data_path, data_name, area):
    DEMprojection(data_path, data_name)
    DEM2PLY(data_name)
    PLY2OBJ(data_name)
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(data_name + '.obj')

    facenum = area * 100 
    ms.apply_filter('simplification_quadric_edge_collapse_decimation', targetfacenum=facenum, preservenormal=True)

    v_num = ms.current_mesh().vertex_number()
    while True:
        ms.apply_filter('remove_duplicate_faces')
        ms.apply_filter('remove_duplicate_vertices')
        ms.apply_filter('remove_unreferenced_vertices')
        ms.apply_filter('repair_non_manifold_edges') 
        v_num -= ms.current_mesh().vertex_number()
        if v_num == 0:
            break
        else:
            v_num = ms.current_mesh().vertex_number() 

    ms.save_current_mesh('offs/' + data_name + '.off')

    # ms.save_current_mesh(data_name + '.off')
    # os.system('rm *.obj')
    print('DEM2OFF finished. The output file is: ', 'offs/' + data_name + '.off')


if __name__ == '__main__':
    print('origin_data: ', sys.argv[1])
    data_name = sys.argv[1].split('/')[-1].split('.')[0]
    area = int(data_name.split('-')[-1])
    data_name = data_name.split('-')[0]
    print('data_name=', data_name, 'area=', area)
    DEM2OFF(sys.argv[1], data_name, area)
