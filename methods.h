//
// Created by huang on 2022/2/7.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_METHODS_H
#define WEIGHTED_DISTANCE_ORACLE_METHODS_H

#include "base.h"
#include "weighted_distance_oracle.h"
#include "greedySpanner.h"
#include "quad.h"
#include "getopt.hpp"

namespace Methods{

    using namespace std;
    using namespace Base;

    pair<vector<float>, tuple<float, float, float> > A2A_SE(ofstream &fout, const Mesh &mesh, const float &err,
                                                       const unsigned &sp_num, const unsigned &q_num, const vector<float> &face_weight){
        srand((unsigned int)time(0));

        map<unsigned, vector<unsigned> > edge_bisector_map = {}, bisector_point_map = {},  face_point_map = {};
        map<unsigned, unsigned> point_face_map = {};
        map<unsigned, Point> point_location_map = {};

        auto index_start = chrono::_V2::system_clock::now();
        auto memory_begin = physical_memory_used_by_process();

//    vector<float> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsJACM(surface_mesh, err, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);

        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(mesh, sp_num, edge_bisector_map,
                                                                                bisector_point_map, point_face_map,
                                                                                point_location_map, face_point_map);
        unsigned num_base_graph_vertices = ret_place_points.second;
        fout << "base graph |V| = " << num_base_graph_vertices << endl;

        vector<unsigned> pid_list = {};  //  partition tree will index all vertices
        for (auto i = 0; i < num_base_graph_vertices; i++){
            pid_list.push_back(i);
        }

        WeightedDistanceOracle::constructBaseGraph(mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
        float l = 8 / err + 6;
        l *= 2; // Enhanced Edges
        fout << "Begin to build partition tree..." << endl;
        tree.constructTree(fout, pid_list, l);
        fout << "Partition tree finished." << endl;

        set<WeightedDistanceOracle::nodePair> node_pairs;
        tree.generateNodePairSet(err, 1, sp_num, node_pairs); // use approximation based on SGP'2013

        fout << "Node pair set size: " << node_pairs.size() << endl;
        fout << "Node pair percentage: " << fixed << setprecision(3) << 1.0 * node_pairs.size() / num_base_graph_vertices / num_base_graph_vertices << endl;

        auto index_end = chrono::_V2::system_clock::now();
        auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
        float index_time = index_duration.count();
        auto memory_end = physical_memory_used_by_process();
        fout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

        vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

        map<unsigned, unsigned> new_id;
        for (auto i = 0; i < leaf_nodes.size(); i++){
            new_id[i] = i;
        }

        Surface_mesh_shortest_path shortest_paths(mesh);
        vector<float> V2V_results = {};
        float inner_query_time = 0.0, inter_query_time = 0.0;

        unsigned percent = 1;
        for (auto i = 0; i < q_num * 3; i++){
            if (i > percent * q_num / 10){
                fout << percent++ << "0% queries finished." << endl;
            }

            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            auto fid_s = A2A_fid[i].first;
            auto fid_t = A2A_fid[i].second;

            auto q_start = chrono::_V2::system_clock::now();
            float oracle_distance = WeightedDistanceOracle::distanceQueryA2A(s, fid_s, t, fid_t,tree,
                                                                              face_point_map, point_location_map, node_pairs);
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            V2V_results.push_back(oracle_distance);
            if (i < q_num){
                inner_query_time += static_cast<float>(q_duration.count());
            }
            else{
                inter_query_time += static_cast<float>(q_duration.count());
            }

        }

        return make_pair(V2V_results, make_tuple(index_time, inner_query_time, inter_query_time));
    }

    pair<vector<float>, tuple<float, float, float> > A2A_LQT(ofstream &fout, const Mesh &mesh, const float &err,
                                                        const unsigned &sp_num, const unsigned &level, const unsigned &q_num,
                                                        const vector<float> &face_weight) {
        srand((int)time(0));

        map<unsigned, vector<unsigned> > edge_bisector_map, bisector_point_map, face_point_map;
        map<unsigned, unsigned> point_face_map;
        map<unsigned, Point> point_location_map;

        auto index_start = chrono::_V2::system_clock::now();
        auto memory_begin = physical_memory_used_by_process();

//    vector<float> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsJACM(surface_mesh, err, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(mesh, sp_num, edge_bisector_map,
                                                                                bisector_point_map, point_face_map,
                                                                                point_location_map, face_point_map);

        fout << "base graph |V| = " << ret_place_points.second << endl;
        Quad::quadTree quad_tree(mesh, face_point_map);
        float rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
        fout << "quad tree level = " << level << endl;
//    fout << "root.Len / Lmin = " << rootLen / min_len << endl;
        for (auto i = 0; i < level; i++){
            quad_tree.buildLevel(mesh, face_point_map);
        }
        set<unsigned> pids;
        auto node = quad_tree.level_nodes[quad_tree.level][0];
        fout << "Leaf maximum side length: " << max(node->x_max - node->x_min, node->y_max - node->y_min) << endl;
        int leaf_boundary_vertices = -1;
        for (auto node: quad_tree.level_nodes[quad_tree.level]){
            leaf_boundary_vertices = max(leaf_boundary_vertices, static_cast<int>(node->boundary_points_id.size()));
            for (auto pid: node->boundary_points_id){
                pids.insert(pid);
            }
        }
        fout << "V = " << mesh.num_vertices() << " | ";
        fout << "LQT-Leaf = " << pids.size() << endl;
        fout << "maximum leaf node boundary vertices = " << leaf_boundary_vertices << endl;
        vector<unsigned> pid_list(pids.begin(), pids.end());

        int num_base_graph_vertices = ret_place_points.second;

        WeightedDistanceOracle::constructBaseGraph(mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
        float l = 8 / err + 10;
        tree.constructTree(fout, pid_list, l);

        set<WeightedDistanceOracle::nodePair> node_pairs;
        tree.generateNodePairSet(err, 1, sp_num, node_pairs); // use approximations based on SGP'2013
        fout << "Node pair set size: " << node_pairs.size() << endl;
        fout << "Node pair percentage: " << fixed << setprecision(2) << 1.0 * node_pairs.size() / ret_place_points.second / ret_place_points.second << endl;

        vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

        map<unsigned, unsigned> new_id;
        auto spanner = Quad::generateSpanner(pid_list, node_pairs, new_id);

        Quad::distancePreprocessing(quad_tree, kSkip::my_base_graph, face_point_map);

        auto index_end = chrono::_V2::system_clock::now();
        auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
        float index_time = index_duration.count();
        auto memory_end = physical_memory_used_by_process();

        fout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

        vector<float> A2A_result = {};
        float inner_query_time = 0.0, inter_query_time = 0.0;
        int same_box = 0, diff_box = 0;

        int percent = 1;
        for (auto i = 0; i < q_num * 3; i++){

            if (i > percent * q_num * 3 / 10){
                fout << percent++ << "0% queries finished." << endl;
            }

            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            auto fid_s = A2A_fid[i].first;
            auto fid_t = A2A_fid[i].second;

//        fout << "s = " << s << " t = " << t << endl;

            auto q_start = chrono::_V2::system_clock::now();
            auto ret = Quad::queryA2A(spanner, kSkip::my_base_graph, face_point_map, tree, node_pairs,
                                      point_location_map,
                                      s, fid_s, t, fid_t,quad_tree, new_id);
            float spanner_distance = ret.first;
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            if (i < q_num){
                inner_query_time += static_cast<float>(q_duration.count());
            }
            else{
                inter_query_time += static_cast<float>(q_duration.count());
            }
            if (ret.second){
                same_box++;
            }
            else{
                diff_box++;
            }
            A2A_result.push_back(spanner_distance);
        }
        fout << fixed << setprecision(3) << "Quad query distribution: same/diff = " << 1.0 * same_box / q_num << " / " << 1.0 * diff_box / q_num << endl;
        return make_pair(A2A_result, make_tuple(index_time, inner_query_time, inter_query_time));
    }

//    pair<vector<float>, pair<float, float> > New_A2A(ofstream &fout, string &file_name, float eps, int point_num, int level, int q_num, int weight_flag, vector<float> &face_weight) {
//        srand((int)time(0));
//        Mesh surface_mesh;
//        ifstream fin(file_name);
//        fin >> surface_mesh;
//
//        map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
//        map<int, int> point_face_map;
//        map<int, Point> point_location_map;
//
//        auto index_start = chrono::_V2::system_clock::now();
//        auto memory_begin = physical_memory_used_by_process();
//
////    vector<float> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
////    auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsJACM(surface_mesh, eps, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
//        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(surface_mesh, point_num, edge_bisector_map,
//                                                                                bisector_point_map, point_face_map,
//                                                                                point_location_map, face_point_map);
//
//        fout << "base graph |V| = " << ret_place_points.second << endl;
//        Quad::quadTree quad_tree(surface_mesh, face_point_map);
//        float rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
//        fout << "quad tree level = " << level << endl;
////    fout << "root.Len / Lmin = " << rootLen / min_len << endl;
//        for (auto i = 0; i < level; i++){
//            quad_tree.buildLevel(surface_mesh, face_point_map);
//        }
//        set<int> pids;
//        auto node = quad_tree.level_nodes[quad_tree.level][0];
//        fout << "Leaf maximum side length: " << max(node->x_max - node->x_min, node->y_max - node->y_min) << endl;
//        int leaf_boundary_vertices = -1;
//        for (auto node: quad_tree.level_nodes[quad_tree.level]){
//            leaf_boundary_vertices = max(leaf_boundary_vertices, static_cast<int>(node->boundary_points_id.size()));
//            for (auto pid: node->boundary_points_id){
//                pids.insert(pid);
//            }
//        }
//        fout << "V = " << surface_mesh.num_vertices() << " | ";
//        fout << "LQT-Leaf = " << pids.size() << endl;
//        fout << "maximum leaf node boundary vertices = " << leaf_boundary_vertices << endl;
//        vector<int> pid_list(pids.begin(), pids.end());
//
//        int num_base_graph_vertices = ret_place_points.second;
//
//        WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
//                                                   bisector_point_map, point_face_map, point_location_map);
//        fout << "Base graph construction finished." << endl;
//        float s_value = 2 / eps + 1;
//
//        map<int, int> new_id;
//        for (auto i = 0; i < pid_list.size(); i++){
//            new_id[pid_list[i]] = i;
//        }
//        kSkip::Graph spanner;
//        spanner.init(static_cast<int>(pid_list.size()));
//        GreedySpanner::setSeparationValue(2 / eps + 1);
//        GreedySpanner::setBoundaryVertexFlag(surface_mesh.num_vertices());
//        GreedySpanner::generateGreedySpanner(pid_list, spanner, new_id);
//        Quad::distancePreprocessing(quad_tree, kSkip::my_base_graph, face_point_map);
//
//
//        auto index_end = chrono::_V2::system_clock::now();
//        auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
//        float index_time = index_duration.count();
//        auto memory_end = physical_memory_used_by_process();
//
//        fout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;
//
//        vector<float> A2A_result = {};
//        float query_time = 0.0;
//        int same_box = 0, diff_box = 0;
//
//        int percent = 1;
//        for (auto i = 0; i < q_num; i++){
//
//            if (i > percent * q_num / 10){
//                fout << percent++ << "0% queries finished." << endl;
//            }
//
//            auto s = A2A_query[i].first;
//            auto t = A2A_query[i].second;
//            auto fid_s = A2A_fid[i].first;
//            auto fid_t = A2A_fid[i].second;
//
////        fout << "s = " << s << " t = " << t << endl;
//
//            auto q_start = chrono::_V2::system_clock::now();
//            auto ret = Quad::queryA2A(surface_mesh, spanner, kSkip::my_base_graph, face_point_map, point_location_map,
//                                      s, fid_s, t, fid_t,quad_tree, new_id);
//            float spanner_distance = ret.first;
//            auto q_end = chrono::_V2::system_clock::now();
//            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
//            query_time += static_cast<float>(q_duration.count());
//            if (ret.second){
//                same_box++;
//            }
//            else{
//                diff_box++;
//            }
//            A2A_result.push_back(spanner_distance);
//        }
//        fout << fixed << setprecision(3) << "Quad query distribution: same/diff = " << 1.0 * same_box / q_num << " / " << 1.0 * diff_box / q_num << endl;
//        return make_pair(A2A_result, make_pair(index_time, query_time));
//    }

    pair<vector<float>, tuple<float, float, float> > A2A_MMP(ofstream &fout, const Mesh &mesh, const AABB_tree &aabb_tree, const unsigned &q_num){
        srand((int)time(0));

        Surface_mesh_shortest_path shortest_paths(mesh);
        vector<float> A2A_result = {};
        float inner_query_time = 0.0, inter_query_time = 0.0;

        unsigned percent = 1;
        for (auto i = 0; i < q_num * 3; i++){
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;

            if (i > percent * q_num * 3 / 10){
                fout << percent++ << "0% queries finished." << endl;
            }

            auto q_start = chrono::_V2::system_clock::now();

            auto location_s = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(s, aabb_tree, mesh);
            auto location_t = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(t, aabb_tree, mesh);
            shortest_paths.add_source_point(location_s.first, location_s.second);
            auto ret_pair = shortest_paths.shortest_distance_to_source_points(location_t.first, location_t.second);
            float mmp_distance = ret_pair.first;
            shortest_paths.remove_all_source_points();

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);

            A2A_result.push_back(mmp_distance);
            if (i < q_num){
                inner_query_time += static_cast<float>(q_duration.count());
            }
            else{
                inter_query_time += static_cast<float>(q_duration.count());
            }
        }

        return make_pair(A2A_result, make_tuple(0.0, inner_query_time, inter_query_time));
    }

    pair<vector<float>, tuple<float, float, float> > A2A_FixedS(ofstream &fout, const Mesh &mesh,
                                                           const AABB_tree &aabb_tree,
                                                           const unsigned &q_num, const unsigned &sp_num,
                                                           const vector<float> &face_weight){

        srand((int)time(0));

        map<unsigned, vector<unsigned> > edge_bisector_map = {}, bisector_point_map = {},  face_point_map = {};
        map<unsigned, unsigned> point_face_map = {};
        map<unsigned, Point> point_location_map = {};

        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(mesh, sp_num, edge_bisector_map,
                                                                                bisector_point_map, point_face_map,
                                                                                point_location_map, face_point_map);
        unsigned num_base_graph_vertices = ret_place_points.second;

        fout << "base graph |V| = " << num_base_graph_vertices << endl;

        WeightedDistanceOracle::constructBaseGraph(mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        vector<float> A2A_result = {};

        float inner_query_time = 0.0, inter_query_time = 0.0;
        unsigned percent = 1;
        for (auto i = 0; i < q_num * 3; i++){
            if (i > percent * q_num * 3 / 10){
                fout << percent++ << "0% queries finished." << endl;
            }
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            unsigned fid_s = A2A_fid[i].first;
            unsigned fid_t = A2A_fid[i].second;
            auto q_start = chrono::_V2::system_clock::now();

            float dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map, point_location_map);

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            if (i < q_num){
                inner_query_time += static_cast<float>(q_duration.count());
            }
            else{
                inter_query_time += static_cast<float>(q_duration.count());
            }

            A2A_result.push_back(dijk_distance);
        }
        return make_pair(A2A_result, make_tuple(0, inner_query_time, inter_query_time));
    }

    pair<vector<float>, tuple<float, float, float> > A2A_UnfixedS(ofstream &fout, const Mesh &mesh, const AABB_tree &aabb_tree,
                                                             const float &err, const unsigned &q_num,
                                                             const vector<float> &face_weight){

        srand((int)time(0));


        map<unsigned, vector<unsigned> > edge_bisector_map = {}, bisector_point_map = {},  face_point_map = {};
        map<unsigned, unsigned> point_face_map = {};
        map<unsigned, Point> point_location_map = {};

        vector<float> gama = WeightedDistanceOracle::getVertexGamma(mesh, face_weight);
        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsJACM(mesh, err, gama, edge_bisector_map,
                                                                               bisector_point_map, point_face_map,
                                                                               point_location_map, face_point_map);

        unsigned num_base_graph_vertices = ret_place_points.second;
        fout << "base graph |V| = " << num_base_graph_vertices << endl;

        WeightedDistanceOracle::constructBaseGraph(mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        vector<float> A2A_result = {};

        float inner_query_time = 0.0, inter_query_time = 0.0;
        unsigned percent = 1;
        for (auto i = 0; i < q_num * 3; i++){
            if (i > percent * q_num * 3 / 10){
                fout << percent++ << "0% queries finished." << endl;
            }
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            unsigned fid_s = A2A_fid[i].first;
            unsigned fid_t = A2A_fid[i].second;
            auto q_start = chrono::_V2::system_clock::now();

            float dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map, point_location_map);

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            if (i < q_num){
                inner_query_time += static_cast<float>(q_duration.count());
            }
            else{
                inter_query_time += static_cast<float>(q_duration.count());
            }

            A2A_result.push_back(dijk_distance);
        }
        return make_pair(A2A_result, make_tuple(0, inner_query_time, inter_query_time));
    }

    pair<vector<float>, tuple<float, float, float> > A2A_KAlgo(ofstream &fout, Mesh &mesh, const unsigned &q_num, const unsigned &K,
                                                          const vector<float> &face_weight){

        srand((int)time(0));

        float l_min = 1e60, theta_m = -1.0;
        for (auto fd: mesh.faces()){
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))){
                auto p0 = mesh.points()[mesh.source(hed)],
                        p1 = mesh.points()[mesh.target(hed)];
                float e_len = sqrt(CGAL::squared_distance(p0, p1));
                if (floatCmp(e_len - l_min) < 0){
                    l_min = e_len;
                }
            }
            if (fd != Mesh::null_face()) {
                vector<int> vids = {};
                auto hed = mesh.halfedge(fd);
                for (int i = 0; i != 3; i++) {
                    int uid = mesh.source(hed).idx();
                    vids.push_back(uid);
                    hed = mesh.next(hed);
                }
                for (int i = 0; i != 3; i++) {
                    auto t_u = CGAL::SM_Vertex_index(vids[i]),
                            t_v = CGAL::SM_Vertex_index(vids[(i + 1) % 3]),
                            t_w = CGAL::SM_Vertex_index(vids[(i + 2) % 3]);
                    float angle = CGAL::approximate_angle(mesh.points()[t_u], mesh.points()[t_v], mesh.points()[t_w]);
//                cout << "angle = " << angle << endl;
                    if (floatCmp(theta_m) < 0 || floatCmp(angle - theta_m) < 0) {
                        theta_m = angle;
                    }
                }
            }
        }
        fout << "theta_m = " << theta_m << endl;

        vector<float> A2A_result = {};
        kSkip::Graph g;
        kSkip::constructMeshGraph(mesh, g);

        float inner_query_time = 0.0, inter_query_time = 0.0;
        unsigned percent = 1;
        for (auto i = 0; i < q_num * 3; i++){
            if (i > percent * q_num * 3 / 10){
                fout << percent++ << "0% queries finished." << endl;
            }
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            unsigned fid_s = A2A_fid[i].first;
            unsigned fid_t = A2A_fid[i].second;
            auto q_start = chrono::_V2::system_clock::now();

//        float dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map, point_location_map);
            float kAlgo_distance = kSkip::computeDistanceBound(mesh, g, K, s, fid_s, t, fid_t, l_min);

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            if (i < q_num){
                inner_query_time += static_cast<float>(q_duration.count());
            }
            else{
                inter_query_time += static_cast<float>(q_duration.count());
            }

            A2A_result.push_back(kAlgo_distance);
        }
        return make_pair(A2A_result, make_tuple(0, inner_query_time, inter_query_time));
    }

//    int findProperQuadLevel(string &file_name, int point_num){
//
//        srand((int)time(0));
//        Mesh surface_mesh;
//        ifstream fin(file_name);
//        fin >> surface_mesh;
//        vector<float> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
//
//        int num_leaf_nodes = 0;
//        int level = 0;
//        while (num_leaf_nodes < floor(0.2 * surface_mesh.num_vertices() )){
//            level++;
//
//            map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
//            map<int, int> point_face_map;
//            map<int, Point> point_location_map;
//
//            auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(surface_mesh, point_num, edge_bisector_map,
//                                                                                    bisector_point_map, point_face_map,
//                                                                                    point_location_map, face_point_map);
//
//            Quad::quadTree quad_tree(surface_mesh, face_point_map);
//            float rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
//            for (auto i = 0; i < level; i++){
//                quad_tree.buildLevel(surface_mesh, face_point_map);
//            }
//            set<int> pids;
//            auto node = quad_tree.level_nodes[quad_tree.level][0];
//            for (auto node: quad_tree.level_nodes[quad_tree.level]){
//                for (auto pid: node->boundary_points_id){
//                    pids.insert(pid);
//                }
//            }
//            num_leaf_nodes = pids.size();
//        }
//        if (num_leaf_nodes >= floor(0.3 * surface_mesh.num_vertices())) level--;
//        cout << "proper level = " << level << endl;
//        return level;
//    }

    void run_new(int argc, char* argv[]){
        bool generate_flag = getarg(0, "--generate");
        string input = getarg("", "--input"),
            output = getarg("", "--output");
        unsigned grid_num = getarg(4, "--grid-num");
        unsigned q_num = getarg(1000, "--query-num");
        bool weighted_flag = getarg(0, "--weighted");
        string method_type = getarg("", "--method");
        float err = getarg(0.2, "--eps");
        unsigned sp_num = getarg(5, "--sp-num");


        Mesh mesh;
        ifstream fin(input);
        fin >> mesh;
        AABB_tree aabb_tree;
        CGAL::Polygon_mesh_processing::build_AABB_tree(mesh, aabb_tree);
        vector<float> mesh_boundary = retrieveMeshBoundary(mesh);

        if (generate_flag){
            generateQueriesA2A(mesh, mesh_boundary, aabb_tree, q_num, grid_num, 1);  // The first q_num queries are inner-box queries
            generateQueriesA2A(mesh, mesh_boundary, aabb_tree, q_num, grid_num, 0);  // The second q_num queries are inter-box queries
            float inner_ratio = generateQueriesA2ANoFlag(mesh, mesh_boundary, aabb_tree, q_num, grid_num); // The last q_num queries are random generated.
            ofstream fout("A2A.query");
            for (auto i = 0; i < A2A_query.size(); i++){
                fout << fixed << setprecision(6) << A2A_query[i].first << " " << A2A_fid[i].first << " " << A2A_query[i].second << " " << A2A_fid[i].second << endl;
            }
            cout << q_num << " A2A queries generate finished." << endl;

            vector<float> face_weight(mesh.num_faces(), 1.0);
            if (weighted_flag) face_weight = generateFaceWeight(mesh.num_faces());
            ofstream fout2("face_weight.query");
            for (auto i = 0; i < face_weight.size(); i++){
                fout2 << fixed << setprecision(6) << face_weight[i] << endl;
            }
            cout << "face weight generate finished." << endl;
            cout <<  fixed << setprecision(3) << q_num << " random queries, inner ratio = " << inner_ratio << ", inter ratio = " << 1 - inner_ratio << endl;
        }
        else{
            ofstream fout(output);
            loadQueriesA2A(A2A_query, A2A_fid);
            fout << "Load A2A queries finished." << endl;
            vector<float> face_weight = {};
            loadFaceWeight(face_weight);
            fout << "Load face weight finished." << endl;
            float inner = 0;
            for (auto i = 2 * q_num; i < A2A_query.size(); i++){
                inner += (A2A_fid[i].first == A2A_fid[i].second) ? 1.0 : 0.0;
            }
            inner = inner / q_num;
            fout <<  fixed << setprecision(3) << q_num << " random queries, inner ratio = " << inner << ", inter ratio = " << 1 - inner << endl;

            if (method_type == "FixedS"){
                fout << "Run Algorithm 0: Bisector-Fixed-Scheme" << endl;
                auto res_bisector_fixed_S = A2A_FixedS(fout, mesh, aabb_tree, q_num, sp_num, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_bisector_fixed_S.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                auto res_time = res_bisector_fixed_S.second;
                fout << fixed << setprecision(3) << "[Inner Query Running Time] Bisector-Fixed-Scheme: " << get<1>(res_time) << " ms" << endl;
                fout << fixed << setprecision(3) << "[Inter Query Running Time] Bisector-Fixed-Scheme: " << get<2>(res_time) << " ms" << endl;
            }
            else if (method_type == "UnfixedS"){
                fout << "Run Algorithm 1: Bisector-Unfixed-Scheme" << endl;
                auto res_bisector_unfixed_S = A2A_UnfixedS(fout, mesh, aabb_tree, err, q_num, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_bisector_unfixed_S.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                auto res_time = res_bisector_unfixed_S.second;
                fout << fixed << setprecision(3) << "[Inner Query Running Time] Bisector-Unfixed-Scheme: " << get<1>(res_time) << " ms" << endl;
                fout << fixed << setprecision(3) << "[Inter Query Running Time] Bisector-Unfixed-Scheme: " << get<2>(res_time) << " ms" << endl;

            }
            else if (method_type == "KAlgo"){
                fout << "Run Algorithm 2: K-Algo" << endl;
                unsigned K = floor(1 / err + 1);
                auto res_KAlgo = A2A_KAlgo(fout, mesh, q_num, K, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_KAlgo.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                auto res_time = res_KAlgo.second;
                fout << fixed << setprecision(3) << "[Inner Query Running Time] K-Algo: " << get<1>(res_time) << " ms" << endl;
                fout << fixed << setprecision(3) << "[Inter Query Running Time] K-Algo: " << get<2>(res_time) << " ms" << endl;

            }
            else if (method_type == "SE"){
                fout << "Run Algorithm 3: SE-Oracle" << endl;
                auto res_SE = A2A_SE(fout, mesh, err, sp_num, q_num, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_SE.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                auto res_time = res_SE.second;
                fout << fixed << setprecision(3) << "[Index Time] SE-oracle: " << get<0>(res_time) << " ms" << endl;
                fout << fixed << setprecision(3) << "[Inner Query Running Time] SE-oracle: " << get<1>(res_time) << " ms" << endl;
                fout << fixed << setprecision(3) << "[Inter Query Running Time] SE-oracle: " << get<2>(res_time) << " ms" << endl;
            }
            else if (method_type == "LQT"){
                fout << "Run Algorithm 4: LQT-oracle" << endl;
                unsigned level = floor(log2(1.0 * grid_num) * 0.5 + eps);
                auto res_LQT = A2A_LQT(fout, mesh, err, sp_num, level, q_num, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_LQT.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                auto res_time = res_LQT.second;
                fout << fixed << setprecision(3) << "[Index Time] LQT-oracle: " << get<0>(res_time) << " ms" << endl;
                fout << fixed << setprecision(3) << "[Inner Query Running Time] LQT-oracle: " << get<1>(res_time) << " ms" << endl;
                fout << fixed << setprecision(3) << "[Inter Query Running Time] LQT-oracle: " << get<2>(res_time) << " ms" << endl;

            }
            else if (method_type == "MMP"){
                fout << "Run Algorithm 5: MMP-Algorithm" << endl;
                auto res_MMP = A2A_MMP(fout, mesh, aabb_tree, q_num);

                fout << "Query results begin: " << endl;
                for (auto dis: res_MMP.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                auto res_time = res_MMP.second;
                fout << fixed << setprecision(3) << "[Inner Query Running Time] MMP-Algo: " << get<1>(res_time) << " ms" << endl;
                fout << fixed << setprecision(3) << "[Inter Query Running Time] MMP-Algo: " << get<2>(res_time) << " ms" << endl;

            }
            else{
                cout << "Method should between 0 and 5." << endl;
            }

        }
    }

    void run(int argc, char* argv[]){
//        int generate_queries, algo_type, q_num, sp_num, lqt_lev, weight_flag;
//        float eps;
//        string mesh_file, output_file;
//        getOpt2(argc, argv, generate_queries, mesh_file, weight_flag, q_num, eps, sp_num, algo_type, lqt_lev, output_file);
//
//        if (lqt_lev < 0 && algo_type == 4){
//            lqt_lev = findProperQuadLevel(mesh_file, sp_num);
//        }
//
//        if (generate_queries){
//            cout << "Generate A2A queries start..." << endl;
//            A2A_query = generateQueriesA2A(mesh_file, q_num, A2A_fid);
//            ofstream fout("A2A.query");
//            for (auto i = 0; i < A2A_query.size(); i++){
//                fout << fixed << setprecision(6) << A2A_query[i].first << " " << A2A_fid[i].first << " " << A2A_query[i].second << " " << A2A_fid[i].second << endl;
//            }
//            cout << q_num << " A2A queries generate finished." << endl;
//            auto face_weight = generateFaceWeight(mesh_file);
//            ofstream fout2("face_weight.query");
//            for (auto i = 0; i < face_weight.size(); i++){
//                fout2 << fixed << setprecision(6) << face_weight[i] << endl;
//            }
//            cout << "face weight generate finished." << endl;
//        }
//        else{
//            string prefix;
//
//            if (weight_flag){
//                switch (algo_type) {
//                    case 0: prefix = "../results/weighted/fixedS/"; break;
//                    case 1: prefix = "../results/weighted/unfixedS/"; break;
//                    case 2: prefix = "../results/weighted/KAlgo/"; break;
//                    case 3: prefix = "../results/weighted/SE/"; break;
//                    case 4: prefix = "../results/weighted/LQT/"; break;
//                    case 5: prefix = "../results/weighted/MMP/"; break;
//                    case 6: prefix = "../results/weighted/greedy/"; break;
//                    default: break;
//                }
//            }
//            else{
//                switch (algo_type) {
//                    case 0: prefix = "../results/unweighted/fixedS/"; break;
//                    case 1: prefix = "../results/unweighted/unfixedS/"; break;
//                    case 2: prefix = "../results/unweighted/KAlgo/"; break;
//                    case 3: prefix = "../results/unweighted/SE/"; break;
//                    case 4: prefix = "../results/unweighted/LQT/"; break;
//                    case 5: prefix = "../results/unweighted/MMP/"; break;
//                    case 6: prefix = "../results/unweighted/greedy/"; break;
//                    default: break;
//                }
//            }
//
//
//            output_file = prefix + output_file;
//            ofstream fout(output_file);
//            fout << "Load A2A queries..." << endl;
//            loadQueriesA2A(A2A_query, A2A_fid);
//            fout << "Load A2A queries finished." << endl;
//
//            vector<float> face_weight = {};
//            fout << "Load face weights..." << endl;
//            loadFaceWeight(face_weight);
//            fout << "Load face weights finished." << endl;
//
//            fout << fixed << setprecision(3) << "eps = " << eps << endl;
//            if (weight_flag){
//                fout << "Terrain type is: Weighted." << endl;
//            }
//            else{
//                fout << "Terrain type is: Unweighted" << endl;
//                for (auto i = 0; i < face_weight.size(); i++){
//                    face_weight[i] = 1.000000;
//                }
//                fout << "Face weight are set to be 1.0" << endl;
//            }
//
//            fout << "Run algorithm " << algo_type;
//            // Algo 0: Bisector Fixed Scheme
//            // Algo 1: Bisector Unfixed Scheme
//            // Algo 2: K-Algo
//            // Algo 3: SE-Oracle
//            // Algo 4: LQT-Oracle
//            // Algo 5: MMP-Algo # approximate construction on unweighted terrain
//            // Algo 6: greedy-Algo # generate spanner by the greedy algorithm
//            if (algo_type == 0){
//                fout << ": Bisector-Fixed-Scheme" << endl;
//                auto res_bisector_fixed_S = A2A_FixedS(fout, mesh_file, q_num, sp_num, weight_flag, face_weight);
//                fout << "Query results begin: " << endl;
//                for (auto dis: res_bisector_fixed_S.first){
//                    fout << fixed << setprecision(3) << dis << endl;
//                }
//                fout << "Query results end..." << endl;
//                fout << fixed << setprecision(3) << "[Running Time] Bisector-Fixed-Scheme: " << res_bisector_fixed_S.second.second << " ms" << endl;
//            }
//            else if (algo_type == 1){
//                fout << ": Bisector-Unfixed-Scheme" << endl;
//                auto res_bisector_unfixed_S = A2A_UnfixedS(fout, mesh_file, eps, q_num, weight_flag, face_weight);
//                fout << "Query results begin: " << endl;
//                for (auto dis: res_bisector_unfixed_S.first){
//                    fout << fixed << setprecision(3) << dis << endl;
//                }
//                fout << "Query results end..." << endl;
//                fout << fixed << setprecision(3) << "[Running Time] Bisector-Unfixed-Scheme: " << res_bisector_unfixed_S.second.second << " ms" << endl;
//            }
//            else if (algo_type == 2){
//                fout << ": K-Algo on the fly" << endl;
//                int K = floor(1 / eps + 1);
//                auto res_KAlgo = A2A_KAlgo(fout, mesh_file, q_num, K);
//                fout << "Query results begin: " << endl;
//                for (auto dis: res_KAlgo.first){
//                    fout << fixed << setprecision(3) << dis << endl;
//                }
//                fout << "Query results end..." << endl;
//                fout << fixed << setprecision(3) << "[Running Time] K-Algo on the fly: " << res_KAlgo.second.second << " ms" << endl;
//            }
//            else if (algo_type == 3){
//                fout << ": SE-oracle" << endl;
//                auto res_SE = A2A_SE(fout, mesh_file, eps, sp_num, q_num, weight_flag);
//                fout << "Query results begin: " << endl;
//                for (auto dis: res_SE.first){
//                    fout << fixed << setprecision(3) << dis << endl;
//                }
//                fout << "Query results end..." << endl;
//                fout << fixed << setprecision(3) << "[Index Time] SE-oracle: " << res_SE.second.first << " ms" << endl;
//                fout << fixed << setprecision(3) << "[Running Time] SE-oracle: " << res_SE.second.second << " ms" << endl;
//            }
//            else if (algo_type == 4){
//                fout << ": LQT-oracle" << endl;
//                auto res_LQT = A2A_LQT(fout, mesh_file, eps, sp_num, lqt_lev, q_num, weight_flag, face_weight);
//                fout << "Query results begin: " << endl;
//                for (auto dis: res_LQT.first){
//                    fout << fixed << setprecision(3) << dis << endl;
//                }
//                fout << "Query results end..." << endl;
//                fout << fixed << setprecision(3) << "[Index Time] LQT-oracle: " << res_LQT.second.first << " ms" << endl;
//                fout << fixed << setprecision(3) << "[Running Time] LQT-oracle: " << res_LQT.second.second << " ms" << endl;
//            }
//            else if (algo_type == 5){
//                fout << ": MMP-Algorithm" << endl;
//                auto res_MMP = A2A_MMP(fout, mesh_file, q_num);
//                fout << "Query results begin: " << endl;
//
//                for (auto dis: res_MMP.first){
//                    fout << fixed << setprecision(3) << dis << endl;
//                }
//                fout << "Query results end..." << endl;
//                fout << fixed << setprecision(3) << "[Running Time] MMP-Algo: " << res_MMP.second.second << " ms" << endl;
//            }
//            else if (algo_type == 6){
//                fout << ": GreedySpanner Algorithm" << endl;
//                auto res_greedy = New_A2A(fout, mesh_file, eps, sp_num, lqt_lev, q_num, weight_flag, face_weight);
//                fout << "Query results begin: " << endl;
//                for (auto dis: res_greedy.first){
//                    fout << fixed << setprecision(3) << dis << endl;
//                }
//                fout << "Query results end..." << endl;
//                fout << fixed << setprecision(3) << "[Index Time] Greedy-Spanner: " << res_greedy.second.first << " ms" << endl;
//                fout << fixed << setprecision(3) << "[Running Time] Greedy-Spanner: " << res_greedy.second.second << " ms" << endl;
//            }
//            else{
//                cout << "Algorithm from 0 to 5!" << endl;
//            }
//        }
    }

}

#endif //WEIGHTED_DISTANCE_ORACLE_METHODS_H
