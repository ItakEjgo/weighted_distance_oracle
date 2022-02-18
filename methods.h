//
// Created by huang on 2022/2/7.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_METHODS_H
#define WEIGHTED_DISTANCE_ORACLE_METHODS_H

#include "base.h"
#include "weighted_distance_oracle.h"
#include "greedySpanner.h"
#include "quad.h"

namespace Methods{

    using namespace std;

    vector<pair<int, int> > V2V_query;
    vector<pair<Base::Point, Base::Point> > A2A_query;
    vector<pair<int, int> > A2A_fid;

    pair<vector<double>, pair<double, double> > SE_A2A(ofstream &fout, string &file_name, double eps, int sp_num, int q_num, int weight_flag){
        srand((int)time(0));
        Base::Mesh surface_mesh;
        ifstream fin(file_name);
        fin >> surface_mesh;
        Base::AABB_tree aabb_tree;
        CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

        vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.

        map<int, vector<int> > edge_bisector_map, bisector_point_map,  face_point_map;
        map<int, int> point_face_map;
        map<int, Base::Point> point_location_map;

        auto index_start = chrono::_V2::system_clock::now();
        auto memory_begin = Base::physical_memory_used_by_process();

//    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsJACM(surface_mesh, eps, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);

        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(surface_mesh, sp_num, edge_bisector_map,
                                                                                bisector_point_map, point_face_map,
                                                                                point_location_map, face_point_map);
        int num_base_graph_vertices = ret_place_points.second;
        fout << "base graph |V| = " << num_base_graph_vertices << endl;

        vector<int> pid_list = {};  //  partition tree will index all vertices
        for (auto i = 0; i < num_base_graph_vertices; i++){
            pid_list.push_back(i);
        }

        WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
        double l = 8 / eps + 6;
        l *= 2;
        fout << "Begin to build partition tree..." << endl;
        tree.constructTree(fout, pid_list, l);
        fout << "Partition tree finished." << endl;

        set<WeightedDistanceOracle::nodePair> node_pairs;
        tree.generateNodePairSet(eps, 1, sp_num, node_pairs); // use approximation based on SGP'2013
        fout << "Node pair set size: " << node_pairs.size() << endl;
        fout << "Node pair percentage: " << fixed << setprecision(3) << 1.0 * node_pairs.size() / num_base_graph_vertices / num_base_graph_vertices << endl;

        auto index_end = chrono::_V2::system_clock::now();
        auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
        double index_time = index_duration.count();
        auto memory_end = Base::physical_memory_used_by_process();
        fout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

        vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

        map<int, int> new_id;
        for (auto i = 0; i < leaf_nodes.size(); i++){
            new_id[i] = i;
        }

        Base::Surface_mesh_shortest_path shortest_paths(surface_mesh);
        vector<double> V2V_results = {};
        double query_time = 0.0;

        int percent = 1;
        for (auto i = 0; i < q_num; i++){
            if (i > percent * q_num / 10){
                fout << percent++ << "0% queries finished." << endl;
            }

            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            auto fid_s = A2A_fid[i].first;
            auto fid_t = A2A_fid[i].second;

            auto q_start = chrono::_V2::system_clock::now();
            double oracle_distance = WeightedDistanceOracle::distanceQueryA2A(s, fid_s, t, fid_t,tree,
                                                                              face_point_map, point_location_map, node_pairs);
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            V2V_results.push_back(oracle_distance);
            query_time += static_cast<double>(q_duration.count());
        }


        return make_pair(V2V_results, make_pair(index_time, query_time));
    }

    pair<vector<double>, pair<double, double> > LQT_A2A(ofstream &fout, string &file_name, double eps, int point_num, int level, int q_num, int weight_flag, vector<double> &face_weight) {
        srand((int)time(0));
        Base::Mesh surface_mesh;
        ifstream fin(file_name);
        fin >> surface_mesh;
        Base::AABB_tree aabb_tree;
        CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

//    vector<double> face_weight = Base::generateFaceWeight(surface_mesh.num_faces());

        map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
        map<int, int> point_face_map;
        map<int, Base::Point> point_location_map;

        auto index_start = chrono::_V2::system_clock::now();
        auto memory_begin = Base::physical_memory_used_by_process();

//    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsJACM(surface_mesh, eps, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(surface_mesh, point_num, edge_bisector_map,
                                                                                bisector_point_map, point_face_map,
                                                                                point_location_map, face_point_map);

        fout << "base graph |V| = " << ret_place_points.second << endl;
        Quad::quadTree quad_tree(surface_mesh, face_point_map);
        double rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
        fout << "quad tree level = " << level << endl;
//    fout << "root.Len / Lmin = " << rootLen / min_len << endl;
        for (auto i = 0; i < level; i++){
            quad_tree.buildLevel(surface_mesh, face_point_map);
        }
        set<int> pids;
        auto node = quad_tree.level_nodes[quad_tree.level][0];
        fout << "Leaf maximum side length: " << max(node->x_max - node->x_min, node->y_max - node->y_min) << endl;
        int leaf_boundary_vertices = -1;
        for (auto node: quad_tree.level_nodes[quad_tree.level]){
            leaf_boundary_vertices = max(leaf_boundary_vertices, static_cast<int>(node->boundary_points_id.size()));
            for (auto pid: node->boundary_points_id){
                pids.insert(pid);
            }
        }
        fout << "V = " << surface_mesh.num_vertices() << " | ";
        fout << "LQT-Leaf = " << pids.size() << endl;
        fout << "maximum leaf node boundary vertices = " << leaf_boundary_vertices << endl;
        vector<int> pid_list(pids.begin(), pids.end());

        int num_base_graph_vertices = ret_place_points.second;

        WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
        double l = 8 / eps + 10;
        tree.constructTree(fout, pid_list, l);

        set<WeightedDistanceOracle::nodePair> node_pairs;
        tree.generateNodePairSet(eps, 1, point_num, node_pairs); // use approximations based on SGP'2013
        fout << "Node pair set size: " << node_pairs.size() << endl;
        fout << "Node pair percentage: " << fixed << setprecision(2) << 1.0 * node_pairs.size() / ret_place_points.second / ret_place_points.second << endl;

        vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

        map<int, int> new_id;
        auto spanner = Quad::generateSpanner(pid_list, node_pairs, new_id);

        Quad::distancePreprocessing(quad_tree, kSkip::my_base_graph, face_point_map);

        auto index_end = chrono::_V2::system_clock::now();
        auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
        double index_time = index_duration.count();
        auto memory_end = Base::physical_memory_used_by_process();

        fout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

        vector<double> A2A_result = {};
        double query_time = 0.0;
        int same_box = 0, diff_box = 0;

        int percent = 1;
        for (auto i = 0; i < q_num; i++){

            if (i > percent * q_num / 10){
                fout << percent++ << "0% queries finished." << endl;
            }

            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            auto fid_s = A2A_fid[i].first;
            auto fid_t = A2A_fid[i].second;

//        fout << "s = " << s << " t = " << t << endl;

            auto q_start = chrono::_V2::system_clock::now();
            auto ret = Quad::queryA2A(surface_mesh, spanner, kSkip::my_base_graph, face_point_map, tree, node_pairs,
                                      point_location_map,
                                      s, fid_s, t, fid_t,quad_tree, new_id);
            double spanner_distance = ret.first;
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            query_time += static_cast<double>(q_duration.count());
            if (ret.second){
                same_box++;
            }
            else{
                diff_box++;
            }
            A2A_result.push_back(spanner_distance);
        }
        fout << fixed << setprecision(3) << "Quad query distribution: same/diff = " << 1.0 * same_box / q_num << " / " << 1.0 * diff_box / q_num << endl;
        return make_pair(A2A_result, make_pair(index_time, query_time));
    }

    pair<vector<double>, pair<double, double> > New_A2A(ofstream &fout, string &file_name, double eps, int point_num, int level, int q_num, int weight_flag, vector<double> &face_weight) {
        srand((int)time(0));
        Base::Mesh surface_mesh;
        ifstream fin(file_name);
        fin >> surface_mesh;

        map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
        map<int, int> point_face_map;
        map<int, Base::Point> point_location_map;

        auto index_start = chrono::_V2::system_clock::now();
        auto memory_begin = Base::physical_memory_used_by_process();

//    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsJACM(surface_mesh, eps, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(surface_mesh, point_num, edge_bisector_map,
                                                                                bisector_point_map, point_face_map,
                                                                                point_location_map, face_point_map);

        fout << "base graph |V| = " << ret_place_points.second << endl;
        Quad::quadTree quad_tree(surface_mesh, face_point_map);
        double rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
        fout << "quad tree level = " << level << endl;
//    fout << "root.Len / Lmin = " << rootLen / min_len << endl;
        for (auto i = 0; i < level; i++){
            quad_tree.buildLevel(surface_mesh, face_point_map);
        }
        set<int> pids;
        auto node = quad_tree.level_nodes[quad_tree.level][0];
        fout << "Leaf maximum side length: " << max(node->x_max - node->x_min, node->y_max - node->y_min) << endl;
        int leaf_boundary_vertices = -1;
        for (auto node: quad_tree.level_nodes[quad_tree.level]){
            leaf_boundary_vertices = max(leaf_boundary_vertices, static_cast<int>(node->boundary_points_id.size()));
            for (auto pid: node->boundary_points_id){
                pids.insert(pid);
            }
        }
        fout << "V = " << surface_mesh.num_vertices() << " | ";
        fout << "LQT-Leaf = " << pids.size() << endl;
        fout << "maximum leaf node boundary vertices = " << leaf_boundary_vertices << endl;
        vector<int> pid_list(pids.begin(), pids.end());

        int num_base_graph_vertices = ret_place_points.second;

        WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);
        fout << "Base graph construction finished." << endl;
        double s_value = 2 / eps + 1;

        map<int, int> new_id;
        for (auto i = 0; i < pid_list.size(); i++){
            new_id[pid_list[i]] = i;
        }
        kSkip::Graph spanner;
        spanner.init(static_cast<int>(pid_list.size()));
        GreedySpanner::setSeparationValue(2 / eps + 1);
        GreedySpanner::setBoundaryVertexFlag(surface_mesh.num_vertices());
        GreedySpanner::generateGreedySpanner(pid_list, spanner, new_id);
        Quad::distancePreprocessing(quad_tree, kSkip::my_base_graph, face_point_map);


        auto index_end = chrono::_V2::system_clock::now();
        auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
        double index_time = index_duration.count();
        auto memory_end = Base::physical_memory_used_by_process();

        fout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

        vector<double> A2A_result = {};
        double query_time = 0.0;
        int same_box = 0, diff_box = 0;

        int percent = 1;
        for (auto i = 0; i < q_num; i++){

            if (i > percent * q_num / 10){
                fout << percent++ << "0% queries finished." << endl;
            }

            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            auto fid_s = A2A_fid[i].first;
            auto fid_t = A2A_fid[i].second;

//        fout << "s = " << s << " t = " << t << endl;

            auto q_start = chrono::_V2::system_clock::now();
            auto ret = Quad::queryA2A(surface_mesh, spanner, kSkip::my_base_graph, face_point_map, point_location_map,
                                      s, fid_s, t, fid_t,quad_tree, new_id);
            double spanner_distance = ret.first;
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            query_time += static_cast<double>(q_duration.count());
            if (ret.second){
                same_box++;
            }
            else{
                diff_box++;
            }
            A2A_result.push_back(spanner_distance);
        }
        fout << fixed << setprecision(3) << "Quad query distribution: same/diff = " << 1.0 * same_box / q_num << " / " << 1.0 * diff_box / q_num << endl;
        return make_pair(A2A_result, make_pair(index_time, query_time));
    }

    pair<vector<double>, pair<double, double> > MMP_A2A(ofstream &fout, string &mesh_file, int q_num){
        srand((int)time(0));
        Base::Mesh surface_mesh;
        ifstream fin(mesh_file);
        fin >> surface_mesh;
        Base::AABB_tree aabb_tree;
        CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

        Base::Surface_mesh_shortest_path shortest_paths(surface_mesh);
        vector<double> A2A_result = {};
        double query_time = 0.0;


        int percent = 1;
        for (auto i = 0; i < q_num; i++){
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;

            if (i > percent * q_num / 10){
                fout << percent++ << "0% queries finished." << endl;
            }

            auto q_start = chrono::_V2::system_clock::now();

            auto location_s = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(s, aabb_tree, surface_mesh);
            auto location_t = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(t, aabb_tree, surface_mesh);
            shortest_paths.add_source_point(location_s.first, location_s.second);
            auto ret_pair = shortest_paths.shortest_distance_to_source_points(location_t.first, location_t.second);
            double mmp_distance = ret_pair.first;
            shortest_paths.remove_all_source_points();

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);

            A2A_result.push_back(mmp_distance);
            query_time += static_cast<double>(q_duration.count());
        }

        return make_pair(A2A_result, make_pair(0.0, query_time));
    }

    pair<vector<double>, pair<double, double> > bisectorFixedScheme(ofstream &fout, string &file_name, int q_num, int point_num, int weight_flag, vector<double> &face_weight){

        srand((int)time(0));
        Base::Mesh surface_mesh;
        ifstream fin(file_name);
        fin >> surface_mesh;
        Base::AABB_tree aabb_tree;
        CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

//    vector<double> face_weight = Base::generateFaceWeight(surface_mesh.num_faces());

        map<int, vector<int> > edge_bisector_map, bisector_point_map,  face_point_map;
        map<int, int> point_face_map;
        map<int, Base::Point> point_location_map;

        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(surface_mesh, point_num, edge_bisector_map,
                                                                                bisector_point_map, point_face_map,
                                                                                point_location_map, face_point_map);
        int num_base_graph_vertices = ret_place_points.second;
        fout << "base graph |V| = " << num_base_graph_vertices << endl;

        WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        vector<double> A2A_result = {};

        double query_time = 0.0;
        int percent = 1;
        for (auto i = 0; i < q_num; i++){
            if (i > percent * q_num / 10){
                fout << percent++ << "0% queries finished." << endl;
            }
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            int fid_s = A2A_fid[i].first;
            int fid_t = A2A_fid[i].second;
            auto q_start = chrono::_V2::system_clock::now();

            double dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map, point_location_map);

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            query_time += static_cast<double>(q_duration.count());

            A2A_result.push_back(dijk_distance);
        }
        return make_pair(A2A_result, make_pair(0, query_time));
    }

    pair<vector<double>, pair<double, double> > bisectorUnfixedScheme(ofstream &fout, string &file_name, double eps, int q_num, int weight_flag, vector<double> &face_weight){

        srand((int)time(0));
        Base::Mesh surface_mesh;
        ifstream fin(file_name);
        fin >> surface_mesh;
        Base::AABB_tree aabb_tree;
        CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

//    vector<double> face_weight = Base::generateFaceWeight(surface_mesh.num_vertices());

        map<int, vector<int> > edge_bisector_map, bisector_point_map,  face_point_map;
        map<int, int> point_face_map;
        map<int, Base::Point> point_location_map;

        vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
        auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsJACM(surface_mesh, eps, gama, edge_bisector_map,
                                                                               bisector_point_map, point_face_map,
                                                                               point_location_map, face_point_map);

        int num_base_graph_vertices = ret_place_points.second;
        fout << "base graph |V| = " << num_base_graph_vertices << endl;

        WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                   bisector_point_map, point_face_map, point_location_map);

        vector<double> A2A_result = {};

        double query_time = 0.0;
        int percent = 1;
        for (auto i = 0; i < q_num; i++){
            if (i > percent * q_num / 10){
                fout << percent++ << "0% queries finished." << endl;
            }
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            int fid_s = A2A_fid[i].first;
            int fid_t = A2A_fid[i].second;
            auto q_start = chrono::_V2::system_clock::now();

            double dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map, point_location_map);

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            query_time += static_cast<double>(q_duration.count());

            A2A_result.push_back(dijk_distance);
        }
        return make_pair(A2A_result, make_pair(0, query_time));
    }

    pair<vector<double>, pair<double, double> > KAlgo_A2A(ofstream &fout, string &file_name, int q_num, int K){

        srand((int)time(0));
        Base::Mesh surface_mesh;
        ifstream fin(file_name);
        fin >> surface_mesh;

        vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.


        double l_min = 1e60, theta_m = -1.0;
        for (auto fd: surface_mesh.faces()){
            for (auto hed: surface_mesh.halfedges_around_face(surface_mesh.halfedge(fd))){
                auto p0 = surface_mesh.points()[surface_mesh.source(hed)],
                        p1 = surface_mesh.points()[surface_mesh.target(hed)];
                double e_len = sqrt(CGAL::squared_distance(p0, p1));
                if (Base::doubleCmp(e_len - l_min) < 0){
                    l_min = e_len;
                }
            }
            if (fd != Base::Mesh::null_face()) {
                vector<int> vids = {};
                auto hed = surface_mesh.halfedge(fd);
                for (int i = 0; i != 3; i++) {
                    int uid = surface_mesh.source(hed).idx();
                    vids.push_back(uid);
                    hed = surface_mesh.next(hed);
                }
                for (int i = 0; i != 3; i++) {
                    auto t_u = CGAL::SM_Vertex_index(vids[i]),
                            t_v = CGAL::SM_Vertex_index(vids[(i + 1) % 3]),
                            t_w = CGAL::SM_Vertex_index(vids[(i + 2) % 3]);
                    double angle = CGAL::approximate_angle(surface_mesh.points()[t_u], surface_mesh.points()[t_v], surface_mesh.points()[t_w]);
//                cout << "angle = " << angle << endl;
                    if (Base::doubleCmp(theta_m) < 0 || Base::doubleCmp(angle - theta_m) < 0) {
                        theta_m = angle;
                    }
                }
            }
        }
        fout << "theta_m = " << theta_m << endl;

        vector<double> A2A_result = {};
        kSkip::Graph g;
        kSkip::constructMeshGraph(surface_mesh, g);

        double query_time = 0.0;
        int percent = 1;
        for (auto i = 0; i < q_num; i++){
            if (i > percent * q_num / 10){
                fout << percent++ << "0% queries finished." << endl;
            }
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            int fid_s = A2A_fid[i].first;
            int fid_t = A2A_fid[i].second;
            auto q_start = chrono::_V2::system_clock::now();

//        double dijk_distance = kSkip::queryGraphA2A(kSkip::my_base_graph, s, fid_s, t, fid_t, face_point_map, point_location_map);
            double kAlgo_distance = kSkip::computeDistanceBound(surface_mesh, g, K, s, fid_s, t, fid_t, l_min);

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            query_time += static_cast<double>(q_duration.count());

            A2A_result.push_back(kAlgo_distance);
        }
        return make_pair(A2A_result, make_pair(0, query_time));
    }

    int findProperQuadLevel(string &file_name, int point_num){

        srand((int)time(0));
        Base::Mesh surface_mesh;
        ifstream fin(file_name);
        fin >> surface_mesh;
        vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.

        int num_leaf_nodes = 0;
        int level = 0;
        while (num_leaf_nodes < floor(0.2 * surface_mesh.num_vertices() )){
            level++;

            map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
            map<int, int> point_face_map;
            map<int, Base::Point> point_location_map;

            auto ret_place_points = WeightedDistanceOracle::placeSteinerPointsFixed(surface_mesh, point_num, edge_bisector_map,
                                                                                    bisector_point_map, point_face_map,
                                                                                    point_location_map, face_point_map);

            Quad::quadTree quad_tree(surface_mesh, face_point_map);
            double rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
            for (auto i = 0; i < level; i++){
                quad_tree.buildLevel(surface_mesh, face_point_map);
            }
            set<int> pids;
            auto node = quad_tree.level_nodes[quad_tree.level][0];
            for (auto node: quad_tree.level_nodes[quad_tree.level]){
                for (auto pid: node->boundary_points_id){
                    pids.insert(pid);
                }
            }
            num_leaf_nodes = pids.size();
        }
        if (num_leaf_nodes >= floor(0.3 * surface_mesh.num_vertices())) level--;
        cout << "proper level = " << level << endl;
        return level;
    }

    void run(int argc, char* argv[]){
        int generate_queries, algo_type, q_num, sp_num, lqt_lev, weight_flag;
        double eps;
        string mesh_file, output_file;
        Base::getOpt2(argc, argv, generate_queries, mesh_file, weight_flag, q_num, eps, sp_num, algo_type, lqt_lev, output_file);

        if (lqt_lev < 0 && algo_type == 4){
            lqt_lev = findProperQuadLevel(mesh_file, sp_num);
        }

        if (generate_queries){
            cout << "Generate A2A queries start..." << endl;
            A2A_query = Base::generateQueriesA2A(mesh_file, q_num, A2A_fid);
            ofstream fout("A2A.query");
            for (auto i = 0; i < A2A_query.size(); i++){
                fout << fixed << setprecision(6) << A2A_query[i].first << " " << A2A_fid[i].first << " " << A2A_query[i].second << " " << A2A_fid[i].second << endl;
            }
            cout << q_num << " A2A queries generate finished." << endl;
            auto face_weight = Base::generateFaceWeight(mesh_file);
            ofstream fout2("face_weight.query");
            for (auto i = 0; i < face_weight.size(); i++){
                fout2 << fixed << setprecision(6) << face_weight[i] << endl;
            }
            cout << "face weight generate finished." << endl;
        }
        else{
            string prefix;

            if (weight_flag){
                switch (algo_type) {
                    case 0: prefix = "../results/weighted/fixedS/"; break;
                    case 1: prefix = "../results/weighted/unfixedS/"; break;
                    case 2: prefix = "../results/weighted/KAlgo/"; break;
                    case 3: prefix = "../results/weighted/SE/"; break;
                    case 4: prefix = "../results/weighted/LQT/"; break;
                    case 5: prefix = "../results/weighted/MMP/"; break;
                    case 6: prefix = "../results/weighted/greedy/"; break;
                    default: break;
                }
            }
            else{
                switch (algo_type) {
                    case 0: prefix = "../results/unweighted/fixedS/"; break;
                    case 1: prefix = "../results/unweighted/unfixedS/"; break;
                    case 2: prefix = "../results/unweighted/KAlgo/"; break;
                    case 3: prefix = "../results/unweighted/SE/"; break;
                    case 4: prefix = "../results/unweighted/LQT/"; break;
                    case 5: prefix = "../results/unweighted/MMP/"; break;
                    case 6: prefix = "../results/unweighted/greedy/"; break;
                    default: break;
                }
            }


            output_file = prefix + output_file;
            ofstream fout(output_file);
            fout << "Load A2A queries..." << endl;
            Base::loadQueriesA2A(A2A_query, A2A_fid);
            fout << "Load A2A queries finished." << endl;

            vector<double> face_weight = {};
            fout << "Load face weights..." << endl;
            Base::loadFaceWeight(face_weight);
            fout << "Load face weights finished." << endl;

            fout << fixed << setprecision(3) << "eps = " << eps << endl;
            if (weight_flag){
                fout << "Terrain type is: Weighted." << endl;
            }
            else{
                fout << "Terrain type is: Unweighted" << endl;
                for (auto i = 0; i < face_weight.size(); i++){
                    face_weight[i] = 1.000000;
                }
                fout << "Face weight are set to be 1.0" << endl;
            }

            fout << "Run algorithm " << algo_type;
            // Algo 0: Bisector Fixed Scheme
            // Algo 1: Bisector Unfixed Scheme
            // Algo 2: K-Algo
            // Algo 3: SE-Oracle
            // Algo 4: LQT-Oracle
            // Algo 5: MMP-Algo # approximate construction on unweighted terrain
            // Algo 6: greedy-Algo # generate spanner by the greedy algorithm
            if (algo_type == 0){
                fout << ": Bisector-Fixed-Scheme" << endl;
                auto res_bisector_fixed_S = bisectorFixedScheme(fout, mesh_file, q_num, sp_num, weight_flag, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_bisector_fixed_S.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                fout << fixed << setprecision(3) << "[Running Time] Bisector-Fixed-Scheme: " << res_bisector_fixed_S.second.second << " ms" << endl;
            }
            else if (algo_type == 1){
                fout << ": Bisector-Unfixed-Scheme" << endl;
                auto res_bisector_unfixed_S = bisectorUnfixedScheme(fout, mesh_file, eps, q_num, weight_flag, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_bisector_unfixed_S.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                fout << fixed << setprecision(3) << "[Running Time] Bisector-Unfixed-Scheme: " << res_bisector_unfixed_S.second.second << " ms" << endl;
            }
            else if (algo_type == 2){
                fout << ": K-Algo on the fly" << endl;
                int K = floor(1 / eps + 1);
                auto res_KAlgo = KAlgo_A2A(fout, mesh_file, q_num, K);
                fout << "Query results begin: " << endl;
                for (auto dis: res_KAlgo.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                fout << fixed << setprecision(3) << "[Running Time] K-Algo on the fly: " << res_KAlgo.second.second << " ms" << endl;
            }
            else if (algo_type == 3){
                fout << ": SE-oracle" << endl;
                auto res_SE = SE_A2A(fout, mesh_file, eps, sp_num, q_num, weight_flag);
                fout << "Query results begin: " << endl;
                for (auto dis: res_SE.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                fout << fixed << setprecision(3) << "[Index Time] SE-oracle: " << res_SE.second.first << " ms" << endl;
                fout << fixed << setprecision(3) << "[Running Time] SE-oracle: " << res_SE.second.second << " ms" << endl;
            }
            else if (algo_type == 4){
                fout << ": LQT-oracle" << endl;
                auto res_LQT = LQT_A2A(fout, mesh_file, eps, sp_num, lqt_lev, q_num, weight_flag, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_LQT.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                fout << fixed << setprecision(3) << "[Index Time] LQT-oracle: " << res_LQT.second.first << " ms" << endl;
                fout << fixed << setprecision(3) << "[Running Time] LQT-oracle: " << res_LQT.second.second << " ms" << endl;
            }
            else if (algo_type == 5){
                fout << ": MMP-Algorithm" << endl;
                auto res_MMP = MMP_A2A(fout, mesh_file, q_num);
                fout << "Query results begin: " << endl;

                for (auto dis: res_MMP.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                fout << fixed << setprecision(3) << "[Running Time] MMP-Algo: " << res_MMP.second.second << " ms" << endl;
            }
            else if (algo_type == 6){
                fout << ": GreedySpanner Algorithm" << endl;
                auto res_greedy = New_A2A(fout, mesh_file, eps, sp_num, lqt_lev, q_num, weight_flag, face_weight);
                fout << "Query results begin: " << endl;
                for (auto dis: res_greedy.first){
                    fout << fixed << setprecision(3) << dis << endl;
                }
                fout << "Query results end..." << endl;
                fout << fixed << setprecision(3) << "[Index Time] Greedy-Spanner: " << res_greedy.second.first << " ms" << endl;
                fout << fixed << setprecision(3) << "[Running Time] Greedy-Spanner: " << res_greedy.second.second << " ms" << endl;
            }
            else{
                cout << "Algorithm from 0 to 5!" << endl;
            }
        }
    }

}

#endif //WEIGHTED_DISTANCE_ORACLE_METHODS_H
