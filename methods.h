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

    // Returned distance ; index_time, query_time for inner-box, inter-box, mixed.
    pair<vector<float>, vector<float> > A2A_SE(ofstream &fout, const Mesh &mesh, const float &err,
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

        fout << fixed << setprecision(6) << "[Index Time] SE-oracle: " << index_time << " ms" << endl;
        fout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

        vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

        map<unsigned, unsigned> new_id;
        for (auto i = 0; i < leaf_nodes.size(); i++){
            new_id[i] = i;
        }

        Surface_mesh_shortest_path shortest_paths(mesh);
        vector<float> V2V_results = {}, res_time = {};
//        fout << "Query results begin: " << endl;

        for (auto i = 0; i < A2A_query.size(); i++){
//            if (i > percent * A2A_query.size() / 10){
//                fout << percent++ << "0% queries finished." << endl;
//            }
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            auto fid_s = A2A_fid[i].first;
            auto fid_t = A2A_fid[i].second;

            auto q_start = chrono::_V2::system_clock::now();
            float oracle_distance = WeightedDistanceOracle::distanceQueryA2A(s, fid_s, t, fid_t,tree,
                                                                              face_point_map, point_location_map, node_pairs);
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);

            V2V_results.emplace_back(oracle_distance);
            res_time.emplace_back(static_cast<float>(q_duration.count()));
            fout << fixed << setprecision(6) << oracle_distance << " " << static_cast<float>(q_duration.count()) << endl;
        }
//        fout << "Query results end. " << endl;
        return make_pair(V2V_results, res_time);
    }

    pair<vector<float>, vector<float> > A2A_EAR(ofstream &fout, const Mesh &mesh, const float &err,
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
        fout << "EAR-Leaf = " << pids.size() << endl;
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

        fout << fixed << setprecision(6) << "[Index Time] EAR-oracle: " << index_time << " ms" << endl;
        fout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

        vector<float> A2A_result = {}, res_time = {};
        float tot_query_construction_time = 0.0, tot_dijkstra_running_time = 0.0;

        float cur_query_time = 0.0;
        unsigned percent = 1;
        unsigned  kappa = floor(2.309 * (sp_num + 1));

//        fout << "Query results begin: " << endl;
        for (auto i = 0; i < A2A_query.size(); i++){
            cout << "query-" << i << " is finished." << endl;

//            if (i > percent * A2A_query.size() / 10){
//                fout << percent++ << "0% queries finished." << endl;
//            }

            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            auto fid_s = A2A_fid[i].first;
            auto fid_t = A2A_fid[i].second;

//        cout << "s = " << s << " t = " << t << endl;
            Quad::dijkstra_running_time = 0.0;

            auto q_start = chrono::_V2::system_clock::now();
            auto ret = Quad::queryA2A(spanner, kSkip::my_base_graph, face_point_map, tree, node_pairs,
                                      point_location_map,
                                      s, fid_s, t, fid_t,quad_tree, new_id, kappa);
            float spanner_distance = ret.first;
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);

            float entire_cur_query_time = static_cast<float>(q_duration.count());
//            if (i >= q_num * 2){
                tot_dijkstra_running_time += Quad::dijkstra_running_time;
                tot_query_construction_time += entire_cur_query_time - Quad::dijkstra_running_time;
//            }
            A2A_result.emplace_back(spanner_distance);
            res_time.emplace_back(entire_cur_query_time);
            // fout << fixed << setprecision(6) << spanner_distance << " " << static_cast<float>(q_duration.count()) << endl;

        }
//        fout << "Query results end. " << endl;
        fout << "Breakdown_construction: " << fixed << setprecision(6) << tot_query_construction_time << " ms" << endl;
        fout << "Breakdown_dijkstra: " << fixed << setprecision(6) << tot_dijkstra_running_time << " ms" << endl;
        return make_pair(A2A_result, res_time);
    }

    pair<vector<float>, vector<float> > A2A_MMP(ofstream &fout, const Mesh &mesh, const AABB_tree &aabb_tree, const unsigned &q_num,
                                                const unsigned &parallel_num, const unsigned &parallel_id){
        srand((int)time(0));

        Surface_mesh_shortest_path shortest_paths(mesh);
        vector<float> A2A_result = {}, res_time = {};

//        fout << "Query results begin: " << endl;
//        for (auto x = 0; x < 3; x++) {
//            unsigned start_id = x * q_num + q_num / parallel_num * parallel_id;
//            for (auto i = start_id; i < start_id + q_num / parallel_num; i++) {
        for (auto i = 0; i < A2A_query.size(); i++){
                auto s = A2A_query[i].first;
                auto t = A2A_query[i].second;

                auto q_start = chrono::_V2::system_clock::now();

                auto location_s = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(s, aabb_tree, mesh);
                auto location_t = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(t, aabb_tree, mesh);
                shortest_paths.add_source_point(location_s.first, location_s.second);
                auto ret_pair = shortest_paths.shortest_distance_to_source_points(location_t.first, location_t.second);
                float mmp_distance = ret_pair.first;
                shortest_paths.remove_all_source_points();

                auto q_end = chrono::_V2::system_clock::now();
                auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);

                A2A_result.emplace_back(mmp_distance);
                res_time.emplace_back(static_cast<float>(q_duration.count()));
                // fout << fixed << setprecision(6) << mmp_distance << " " << static_cast<float>(q_duration.count()) << endl;
//            }
        }
//        fout << "Query results end. " << endl;


        return make_pair(A2A_result, res_time);
    }

    pair<vector<float>, vector<float> > A2A_FixedS(ofstream &fout, const Mesh &mesh,
                                                           const AABB_tree &aabb_tree,
                                                           const unsigned &q_num, const unsigned &sp_num,
                                                           const vector<float> &face_weight,
                                                           const unsigned &parallel_num,
                                                           const unsigned &parallel_id){

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

        vector<float> A2A_result = {}, res_time = {};

//        fout << "Query results begin: " << endl;
        // for (auto x = 0; x < 1; x++){
        //     unsigned start_id = x * q_num + q_num / parallel_num * parallel_id;
        //     for (auto i = start_id; i < start_id + q_num / parallel_num; i++){
       for (auto i = 0; i < A2A_query.size(); i++){
                auto s = A2A_query[i].first;
                auto t = A2A_query[i].second;
                unsigned fid_s = A2A_fid[i].first;
                unsigned fid_t = A2A_fid[i].second;
                auto q_start = chrono::_V2::system_clock::now();

                float dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map, point_location_map);

                auto q_end = chrono::_V2::system_clock::now();
                auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);

                A2A_result.emplace_back(dijk_distance);
                res_time.emplace_back(static_cast<float>(q_duration.count()));
                // fout << fixed << setprecision(6) << A2A_result[i] << " " << res_time[i] << endl;
            // }
        }
//        fout << "Query results end. " << endl;

        return make_pair(A2A_result, res_time);
    }

    pair<vector<float>, vector<float> > A2A_UnfixedSOnTheFly(ofstream &fout, const Mesh &mesh, const AABB_tree &aabb_tree,
                                                     const float &err, const unsigned &q_num,
                                                     const vector<float> &face_weight,
                                                     const unsigned &parallel_num,
                                                     const unsigned &parallel_id){

        srand((int)time(0));

        vector<float> gama = WeightedDistanceOracle::getVertexGamma(mesh, face_weight);

        kSkip::Graph g;
        kSkip::constructMeshGraph(mesh, g);

        vector<float> A2A_result = {}, res_time = {};

//        fout << "Query results begin: " << endl;
        // for (auto x = 0; x < 1; x++) {
        //     unsigned start_id = x * q_num + q_num / parallel_num * parallel_id;
        //     for (auto i = start_id; i < start_id + q_num / parallel_num; i++) {
       for (auto i = 0; i < A2A_query.size(); i++){
                auto s = A2A_query[i].first;
                auto t = A2A_query[i].second;
                unsigned fid_s = A2A_fid[i].first;
                unsigned fid_t = A2A_fid[i].second;
                auto q_start = chrono::_V2::system_clock::now();

//                float dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map,
//                                                           point_location_map);
                float dijk_distance = kSkip::unfixedOnTheFly(mesh, g, err, gama, face_weight, s, fid_s, t, fid_t);

                auto q_end = chrono::_V2::system_clock::now();
                auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);

                A2A_result.emplace_back(dijk_distance);
                res_time.emplace_back(static_cast<float>(q_duration.count()));
                // fout << fixed << setprecision(6) << dijk_distance << " " << static_cast<float>(q_duration.count()) << endl;
            // }
        }
//        fout << "Query results end. " << endl;

        return make_pair(A2A_result, res_time);
    }

    pair<vector<float>, vector<float> > A2A_UnfixedS(ofstream &fout, const Mesh &mesh, const AABB_tree &aabb_tree,
                                                             const float &err, const unsigned &q_num,
                                                             const vector<float> &face_weight,
                                                             const unsigned &parallel_num,
                                                             const unsigned &parallel_id){

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

        vector<float> A2A_result = {}, res_time = {0.0};



        // fout << "Query results begin: " << endl;
        // for (auto x = 0; x < 3; x++) {
        //     unsigned start_id = x * q_num + q_num / parallel_num * parallel_id;
            // for (auto i = start_id; i < start_id + q_num / parallel_num; i++) {
        for (auto i = 0; i < A2A_query.size(); i++){

                auto s = A2A_query[i].first;
                auto t = A2A_query[i].second;
                unsigned fid_s = A2A_fid[i].first;
                unsigned fid_t = A2A_fid[i].second;
                auto q_start = chrono::_V2::system_clock::now();

                float dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map,
                                                           point_location_map);

                auto q_end = chrono::_V2::system_clock::now();
                auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);

                A2A_result.push_back(dijk_distance);
                // fout << fixed << setprecision(6) << dijk_distance << " " << static_cast<float>(q_duration.count()) << endl;
        }
        // }
        // fout << "Query results end. " << endl;

        return make_pair(A2A_result, res_time);
    }

    pair<vector<float>, vector<float> > A2A_KAlgo(ofstream &fout, Mesh &mesh, const unsigned &q_num, const unsigned &K,
                                                          const vector<float> &face_weight, const unsigned &parallel_num,
                                                          const unsigned &parallel_id){

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

        vector<float> A2A_result = {}, res_time = {};
        kSkip::Graph g;
        kSkip::constructMeshGraph(mesh, g);

//        fout << "Query results begin: " << endl;

        // for (auto x = 0; x < 1; x++) {
        //     unsigned start_id = x * q_num + q_num / parallel_num * parallel_id;
            // for (auto i = start_id; i < start_id + q_num / parallel_num; i++) {
       for (auto i = 0; i < A2A_query.size(); i++){
                auto s = A2A_query[i].first;
                auto t = A2A_query[i].second;
                unsigned fid_s = A2A_fid[i].first;
                unsigned fid_t = A2A_fid[i].second;
                auto q_start = chrono::_V2::system_clock::now();

                float kAlgo_distance = kSkip::computeDistanceBound(mesh, g, K, s, fid_s, t, fid_t, l_min);

                auto q_end = chrono::_V2::system_clock::now();
                auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);

                A2A_result.emplace_back(kAlgo_distance);
                res_time.emplace_back(static_cast<float>(q_duration.count()));
                // fout << fixed << setprecision(6) << kAlgo_distance << " " << static_cast<float>(q_duration.count()) << endl;
            // }
        }
//        fout << "Query results end. " << endl;

        return make_pair(A2A_result, res_time);
    }

    unsigned printHighwayNodeInfo(Base::Mesh &mesh, unsigned point_num, unsigned level){

        vector<float> face_weight(mesh.num_faces(), 1.0); // face weight for each face.

        int num_leaf_nodes = 0;


        map<unsigned, vector<unsigned> > edge_bisector_map, bisector_point_map, face_point_map;
        map<unsigned, unsigned> point_face_map;
        map<unsigned, Point> point_location_map;

//        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(mesh, point_num, edge_bisector_map,
//                                                                                bisector_point_map, point_face_map,
//                                                                                point_location_map, face_point_map);

        Quad::quadTree quad_tree(mesh, face_point_map);
//        float rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
        for (auto i = 0; i < level; i++){
            quad_tree.buildLevel(mesh, face_point_map);
        }
        set<int> pids;
        for (auto node: quad_tree.level_nodes[quad_tree.level]){
            for (auto pid: node->boundary_points_id){
                pids.insert(pid);
            }
        }
        num_leaf_nodes = pids.size();


        cout << "zeta = " << pow(4, level) << endl;
        cout << "highway nodes = " << num_leaf_nodes << endl;
        cout << "ratio = " << 1.0 * num_leaf_nodes / point_num << endl;
        return level;
    }

    void distanceGapQueryGeneration(ofstream &fout, const unsigned &pnum, Base::Mesh &mesh, const AABB_tree &aabb_tree,
                                    const float &err, const unsigned &sp_num, const unsigned &level, const unsigned &q_num,
                                    const vector<float> &face_weight){
        srand((int)time(0));

        map<unsigned, vector<unsigned> > edge_bisector_map, bisector_point_map, face_point_map;
        map<unsigned, unsigned> point_face_map;
        map<unsigned, Point> point_location_map;

        auto index_start = chrono::_V2::system_clock::now();
        auto memory_begin = physical_memory_used_by_process();

        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(mesh, sp_num, edge_bisector_map,
                                                                                bisector_point_map, point_face_map,
                                                                                point_location_map, face_point_map);

        cout << "base graph |V| = " << ret_place_points.second << endl;
        Quad::quadTree quad_tree(mesh, face_point_map);
        float rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
        cout << "quad tree level = " << level << endl;

        for (auto i = 0; i < level; i++){
            quad_tree.buildLevel(mesh, face_point_map);
        }
        set<unsigned> pids;
        auto node = quad_tree.level_nodes[quad_tree.level][0];
        cout << "Leaf maximum side length: " << max(node->x_max - node->x_min, node->y_max - node->y_min) << endl;
        int leaf_boundary_vertices = -1;
        for (auto node: quad_tree.level_nodes[quad_tree.level]){
            leaf_boundary_vertices = max(leaf_boundary_vertices, static_cast<int>(node->boundary_points_id.size()));
            for (auto pid: node->boundary_points_id){
                pids.insert(pid);
            }
        }
        cout << "V = " << mesh.num_vertices() << " | ";
        cout << "EAR-Leaf = " << pids.size() << endl;
        cout << "maximum leaf node boundary vertices = " << leaf_boundary_vertices << endl;
        vector<unsigned> pid_list(pids.begin(), pids.end());

        int num_base_graph_vertices = ret_place_points.second;

        WeightedDistanceOracle::constructBaseGraph(mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
        float l = 8 / err + 10;
        tree.constructTree(fout, pid_list, l);

        set<WeightedDistanceOracle::nodePair> node_pairs;
        tree.generateNodePairSet(err, 1, sp_num, node_pairs); // use approximations based on SGP'2013
        cout << "Node pair set size: " << node_pairs.size() << endl;
        cout << "Node pair percentage: " << fixed << setprecision(2) << 1.0 * node_pairs.size() / ret_place_points.second / ret_place_points.second << endl;

        vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

        map<unsigned, unsigned> new_id;
        auto spanner = Quad::generateSpanner(pid_list, node_pairs, new_id);

        Quad::distancePreprocessing(quad_tree, kSkip::my_base_graph, face_point_map);

        auto index_end = chrono::_V2::system_clock::now();
        auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
        float index_time = index_duration.count();
        auto memory_end = physical_memory_used_by_process();

        cout << fixed << setprecision(6) << "[Index Time] EAR-oracle: " << index_time << " ms" << endl;
        cout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

        vector<float> A2A_result = {}, res_time = {};
        float tot_query_construction_time = 0.0, tot_dijkstra_running_time = 0.0;

        float cur_query_time = 0.0;
        unsigned percent = 1;
        unsigned  kappa = floor(2.309 * (sp_num + 1));


        vector<unsigned> sids = {};
        for (auto i = 0; i != pnum; i++){
            sids.emplace_back(rand() % mesh.num_vertices());
        }

        float dis_gap = 2000;
        struct Q{
            Base::Point p;
            unsigned fid;
            Q(Base::Point a, unsigned b){
                p = a;
                fid = b;
            }
        };
        vector<vector<pair<Q, Q> > > queries(10);

        for (auto i = 0; i != pnum; i++) {
            auto vds = *(mesh.vertices_begin() + sids[i]);
            for (auto j = 0; j < mesh.num_vertices(); j++) {

                unsigned t_sum = 0;
                for (auto k = 0; k < 10; k++){
                    t_sum += queries[k].size();
                }
                if (t_sum == (i + 1) * 100){
                    break;
                }

                auto vdt = *(mesh.vertices_begin() + j);
                auto s = mesh.points()[vds];
                auto t = mesh.points()[vdt];
                auto location_s = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(s, aabb_tree, mesh);
                auto location_t = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(t, aabb_tree, mesh);

                auto ret = Quad::queryA2A(spanner, kSkip::my_base_graph, face_point_map, tree, node_pairs,
                                          point_location_map,
                                          s, location_s.first.idx(), t, location_t.first.idx(),quad_tree, new_id, kappa);
                float spanner_distance = ret.first;


                for (auto k = 0; k < 10; k++){
//                    cout << "query-" << k << " size = " << queries[k].size() << endl;
                    if (queries[k].size() < (i + 1) * 10 && Base::floatCmp(spanner_distance - (k + 1) * dis_gap) <= 0 && Base::floatCmp(spanner_distance - k * dis_gap) >= 0){
                        queries[k].emplace_back(Q(s, location_s.first.idx()), Q(t, location_t.first.idx()));
                        cout << "add 1 query to query-" << k << ", current size = " << queries[k].size() << endl;
                        break;
                    }
                }
            }
        }

        for (auto i = 0; i < 10; i++){
            while (queries[i].size() < 100){
                queries[i].emplace_back(queries[i][0]);
            }
            for (auto j = 0; j < queries[i].size(); j++){
                fout << fixed << setprecision(3) << queries[i][j].first.p << " " << queries[i][j].first.fid << " " << queries[i][j].second.p << " " << queries[i][j].second.fid << endl;
            }
        }

        ofstream fout2("face_weight.query");
        for (auto i = 0; i < face_weight.size(); i++){
            fout2 << fixed << setprecision(6) << face_weight[i] << endl;
        }
        cout << "face weight generate finished." << endl;
    }

    void run_new(int argc, char* argv[]){
        bool generate_flag = getarg(0, "--generate");
        string input = getarg("", "--input"),
            output = getarg("", "--output");
        unsigned grid_num = getarg(4, "--grid-num");
        unsigned q_num = getarg(100, "--query-num");
        bool weighted_flag = getarg(0, "--weighted");
        string method_type = getarg("EAR", "--method");
        float err = getarg(0.2, "--eps");
        unsigned sp_num = getarg(5, "--sp-num");
        unsigned parallel_num = getarg(1, "--parallel-num");
        unsigned parallel_id = getarg(0, "--parallel-id");

        Mesh mesh;
        ifstream fin(input);
        fin >> mesh;
        cout << mesh.num_vertices() << " " << mesh.num_edges() << " " << mesh.num_faces() << endl;
        AABB_tree aabb_tree;
        CGAL::Polygon_mesh_processing::build_AABB_tree(mesh, aabb_tree);
        vector<float> mesh_boundary = retrieveMeshBoundary(mesh);
        bool distance_gap_flag = false;
        if (generate_flag){
            ofstream fout("A2A.query");
            float inner_ratio;
            if (!distance_gap_flag){
                generateQueriesA2A(mesh, mesh_boundary, aabb_tree, q_num, grid_num, 1);  // The first q_num queries are inner-box queries
                generateQueriesA2A(mesh, mesh_boundary, aabb_tree, q_num, grid_num, 0);  // The second q_num queries are inter-box queries
                inner_ratio = generateQueriesA2ANoFlag(mesh, mesh_boundary, aabb_tree, q_num, grid_num); // The last q_num queries are random generated.
            }
            else{
                unsigned level = floor(log2(1.0 * grid_num) * 0.5 + eps);
                vector<float> face_weight(mesh.num_faces(), 1.0);
                distanceGapQueryGeneration(fout, 10, mesh, aabb_tree, err, sp_num, level, q_num, face_weight);
                return ;
            }
//            Base::generateQueriesA2A(mesh, q_num);
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
//            cout <<  fixed << setprecision(3) << q_num << " random queries, inner ratio = " << inner_ratio << ", inter ratio = " << 1 - inner_ratio << endl;
        }
        else{
            ofstream fout(output);
            loadQueriesA2A(A2A_query, A2A_fid);
            cout << "Load A2A queries finished." << endl;
            fout << "Load A2A queries finished." << endl;
            vector<float> face_weight = {};
            loadFaceWeight(face_weight);
            fout << "Load face weight finished." << endl;

            unsigned level = floor(log2(1.0 * grid_num) * 0.5 + eps);
//            printHighwayNodeInfo(mesh, mesh.num_vertices(), level);
//            return;

            float inner = 0;
            for (auto i = 2 * q_num; i < A2A_query.size(); i++){
                inner += (A2A_fid[i].first == A2A_fid[i].second) ? 1.0 : 0.0;
            }
            inner = inner / q_num;
            fout <<  fixed << setprecision(3) << q_num << " random queries, inner ratio = " << inner << ", inter ratio = " << 1 - inner << endl;
            pair<vector<float>, vector<float> > res;
            if (method_type == "FixedS"){
                fout << "Run Algorithm 0: Bisector-Fixed-Scheme" << endl;
                res = A2A_FixedS(fout, mesh, aabb_tree, q_num, sp_num, face_weight, parallel_num, parallel_id);
            }
            else if (method_type == "UnfixedS"){
                fout << "Run Algorithm 1: Bisector-Unfixed-Scheme" << endl;
//                res = A2A_UnfixedS(fout, mesh, aabb_tree, err, q_num, face_weight, parallel_num, parallel_id);
                res = A2A_UnfixedSOnTheFly(fout, mesh, aabb_tree, err, q_num, face_weight, parallel_num, parallel_id);
            }
            else if (method_type == "KAlgo"){
                fout << "Run Algorithm 2: K-Algo" << endl;
                unsigned K = ceil(1.0 / err + 1);
                cout << fixed << setprecision(2) << "eps = " << err << ", K = " << K << endl;
                res = A2A_KAlgo(fout, mesh, q_num, K, face_weight, parallel_num, parallel_id);
            }
            else if (method_type == "SE"){
                fout << "Run Algorithm 3: SE-Oracle" << endl;
                res = A2A_SE(fout, mesh, err, sp_num, q_num, face_weight);
            }
            else if (method_type == "EAR"){
                fout << "Run Algorithm 4: EAR-oracle" << endl;
                unsigned level = floor(log2(1.0 * grid_num) * 0.5 + eps);
                res = A2A_EAR(fout, mesh, err, sp_num, level, q_num, face_weight);
            }
            else if (method_type == "MMP"){
                fout << "Run Algorithm 5: MMP-Algorithm" << endl;
                res = A2A_MMP(fout, mesh, aabb_tree, q_num, parallel_num, parallel_id);
            }
            else{
                cout << "Method should between 0 and 5." << endl;
            }
           fout << "Query results begin: " << endl;
           for (auto i = 0; i != res.first.size(); i++){
               fout << fixed << setprecision(6) << res.first[i] << " " << res.second[i] << endl;
           }
           fout << "Query results end." << endl;
        }
    }
}

#endif //WEIGHTED_DISTANCE_ORACLE_METHODS_H
