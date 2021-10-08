#include "base.h"
#include "weighted_distance_oracle.h"
#include "highway_well_separated.h"
#include "k_skip.h"
#include "quad.h"

using namespace std;
int K;
vector<pair<int, int> > V2V_query;
vector<pair<Base::Point, Base::Point> > A2A_query;

void generate_queriesA2A(string &file_name, double eps, int type, int point_num, int q_num) {
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
    Base::AABB_tree aabb_tree;
    CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);
    for (auto i = 0; i < q_num; i++){
        auto p1_pair = Base::generateArbitrarySurfacePoint(surface_mesh, aabb_tree);
        auto p2_pair = Base::generateArbitrarySurfacePoint(surface_mesh, aabb_tree);
        A2A_query.emplace_back(p1_pair.first, p2_pair.first);
    }
}

void generate_queriesV2V(string &file_name, double eps, int type, int point_num, int q_num){
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;

    for (auto i = 0; i < q_num; i++){
        int s = rand() % surface_mesh.num_vertices(),
                t = rand() % surface_mesh.num_vertices();
        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
        V2V_query.push_back(make_pair(s, t));
    }
}

pair<vector<double>, pair<double, double> > SE_A2A(string &file_name, double eps, int type, int point_num, int q_num, int q_type){
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
    Base::AABB_tree aabb_tree;
    CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
    vector<double> face_max_length = {}; // maximum edge length for each face.
    WeightedDistanceOracle::getFaceMaxLength(surface_mesh, face_max_length); // get the maximum edge length for each face.

    map<int, vector<int> > edge_bisector_map, bisector_point_map,  face_point_map;
    map<int, int> point_face_map;
    map<int, Base::Point> point_location_map;

    auto index_start = chrono::_V2::system_clock::now();
    auto memory_begin = Base::physical_memory_used_by_process();
//    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsJACM(surface_mesh, eps, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);

    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
    int num_base_graph_vertices = ret_place_points.second;
    cout << "base graph |V| = " << num_base_graph_vertices << endl;

    vector<int> pid_list = {};  //  partition tree will index all vertices
    for (auto i = 0; i < num_base_graph_vertices; i++){
        pid_list.push_back(i);
    }

    WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                                     bisector_point_map, point_face_map, point_location_map);

    WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
    double l = 8 / eps + 10;
    tree.constructTree(pid_list, l);

    set<WeightedDistanceOracle::nodePair> node_pairs;
    tree.generateNodePairSet(eps, type, point_num, node_pairs);
    cout << "Node pair set size: " << node_pairs.size() << endl;
    cout << "Node pair percentage: " << fixed << setprecision(3) << 1.0 * node_pairs.size() / num_base_graph_vertices / num_base_graph_vertices << endl;

    auto index_end = chrono::_V2::system_clock::now();
    auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
    double index_time = index_duration.count();
    auto memory_end = Base::physical_memory_used_by_process();
    cout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

    vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());
//    for (auto i = 0; i < leaf_nodes.size(); i++){
//        cout << "leaf i: " << i << "-" << leaf_nodes[i]->center_idx << endl;
//    }
    map<int, int> new_id;
    for (auto i = 0; i < leaf_nodes.size(); i++){
        new_id[i] = i;
    }

    Base::Surface_mesh_shortest_path shortest_paths(surface_mesh);
    vector<double> V2V_results = {};
    double query_time = 0.0;
    if (!q_type){
        for (auto i = 0; i < q_num; i++){
            int s = V2V_query[i].first;
            int t = V2V_query[i].second;
            assert(leaf_nodes[s]->center_idx == s);
            assert(leaf_nodes[t]->center_idx == t);

            auto q_start = chrono::_V2::system_clock::now();
            vector<WeightedDistanceOracle::PartitionTreeNode*> As, At;
            tree.getPathToRoot(leaf_nodes[s], As);
            tree.getPathToRoot(leaf_nodes[t], At);
            double oracle_distance = WeightedDistanceOracle::distanceQueryBf(node_pairs, As, At);
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            V2V_results.push_back(oracle_distance);
            query_time += static_cast<double>(q_duration.count());
        }
    }
    else{
        for (auto i = 0; i < q_num; i++){
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;

            auto q_start = chrono::_V2::system_clock::now();
            double oracle_distance = WeightedDistanceOracle::distanceQueryA2A(s, t, surface_mesh, tree, kSkip::my_base_graph,
                                                                              face_point_map, point_location_map, node_pairs, aabb_tree, new_id);
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
            V2V_results.push_back(oracle_distance);
            query_time += static_cast<double>(q_duration.count());
        }
    }

    return make_pair(V2V_results, make_pair(index_time, query_time));
}

pair<vector<double>, pair<double, double> > quad_A2A(string &file_name, double eps, int point_num, int type, int level, int q_num, int q_type) {
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
    Base::AABB_tree aabb_tree;
    CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);
    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
    vector<double> face_max_length = {}; // maximum edge length for each face.
    vector<double> face_min_length = {}; // minimum edge length for each face.

    WeightedDistanceOracle::getFaceMaxMinLength(surface_mesh, face_max_length, face_min_length); // get the maximum edge length for each face.

    double max_len = -1.0, min_len = 1e60;
    for (auto val: face_max_length){
        max_len = max(max_len, val);
        min_len = min(min_len, val);
    }
    cout << "max edge length: " << max_len << endl;
    cout << "min edge length: " << min_len << endl;


    WeightedDistanceOracle::getFaceMaxLength(surface_mesh, face_max_length); // get the maximum edge length for each face.
    map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
    map<int, int> point_face_map;
    map<int, Base::Point> point_location_map;

    auto index_start = chrono::_V2::system_clock::now();

    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsJACM(surface_mesh, eps, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
    auto memory_begin = Base::physical_memory_used_by_process();
//    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
    cout << "base graph |V| = " << ret_place_points.second << endl;
    Quad::quadTree quad_tree(surface_mesh, face_point_map);
    double rootLen = max(quad_tree.root->x_max - quad_tree.root->x_min, quad_tree.root->y_max - quad_tree.root->y_min);
    cout << "quad tree level = " << level << endl;
    cout << "root.Len / Lmin = " << rootLen / min_len << endl;
    for (auto i = 0; i < level; i++){
        quad_tree.buildLevel(surface_mesh, face_point_map);
    }
    set<int> pids;
    auto node = quad_tree.level_nodes[quad_tree.level][0];
    cout << "Leaf maximum side length: " << max(node->x_max - node->x_min, node->y_max - node->y_min) << endl;
    for (auto node: quad_tree.level_nodes[quad_tree.level]){
        for (auto pid: node->boundary_points_id){
            pids.insert(pid);
        }
    }
    cout << "V = " << surface_mesh.num_vertices() << endl;
    cout << "Leaf = " << pids.size() << endl;
    vector<int> pid_list(pids.begin(), pids.end());

    int num_base_graph_vertices = ret_place_points.second;

    WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                                     bisector_point_map, point_face_map, point_location_map);

    WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
    double l = 8 / eps + 10;
    tree.constructTree(pid_list, l);

    set<WeightedDistanceOracle::nodePair> node_pairs;
    tree.generateNodePairSet(eps, type, point_num, node_pairs);
    cout << "Node pair set size: " << node_pairs.size() << endl;
    cout << "Node pair percentage: " << fixed << setprecision(2) << 1.0 * node_pairs.size() / ret_place_points.second / ret_place_points.second << endl;

    vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

    map<int, int> new_id;
    auto spanner = Quad::generateSpanner(pid_list, node_pairs, new_id);

    Quad::distancePreprocessing(quad_tree, kSkip::my_base_graph, face_point_map);

    auto index_end = chrono::_V2::system_clock::now();
    auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
    double index_time = index_duration.count();

    auto memory_end = Base::physical_memory_used_by_process();

    cout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;
    vector<double> V2V_result0 = {}, V2V_result1 = {};
    double query_time = 0.0;
    int same_box = 0, diff_box = 0;
    double check_first_time = 0;

    if (!q_type){
        for (auto i = 0; i < q_num; i++){
            int s = V2V_query[i].first;
            int t = V2V_query[i].second;
            auto q_start = chrono::_V2::system_clock::now();
            auto ret = Quad::querySpanner(surface_mesh, spanner, s, t, quad_tree, new_id, tree, aabb_tree, node_pairs);
//        auto ret = Quad::queryWSPD(surface_mesh, quad_tree, s, t, node_pairs, tree, new_id);
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
            V2V_result0.push_back(spanner_distance);
        }
    }
    else{
        Quad::WSPD_hit = 0;
        for (auto i = 0; i < q_num; i++){
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;

            auto q_start = chrono::_V2::system_clock::now();
            auto ret = Quad::queryA2A(surface_mesh, spanner, kSkip::my_base_graph, face_point_map, point_location_map, s, t,
                                      quad_tree, new_id, tree, node_pairs, aabb_tree, 0);
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
            V2V_result0.push_back(spanner_distance);
        }
    }
    cout << fixed << setprecision(3) << "check first quad query: " << check_first_time << " ms. Hit ratio = " << 1.0 * Quad::WSPD_hit / q_num << endl;
    cout << fixed << setprecision(3) << "Quad query distribution: same/diff = " << 1.0 * same_box / q_num << " / " << 1.0 * diff_box / q_num << endl;
    return make_pair(V2V_result0, make_pair(index_time, query_time));
}

pair<vector<double>, pair<double, double> > base_graph_run(int q_num, int q_type){
    vector<double> V2V_results = {};
    if (q_type){
        V2V_results.resize(q_num, -1.0);
        return make_pair(V2V_results, make_pair(0, 0.0));
    }
    double query_time = 0.0;
    for (auto i = 0; i < q_num; i++){
        int s = V2V_query[i].first;
        int t = V2V_query[i].second;
        auto q_start = chrono::_V2::system_clock::now();
        double dijk_distance = kSkip::dijkstra(kSkip::my_base_graph, s, t).first;
        auto q_end = chrono::_V2::system_clock::now();
        auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
        query_time += static_cast<double>(q_duration.count());

        V2V_results.push_back(dijk_distance);

    }
    return make_pair(V2V_results, make_pair(0, query_time));
}

pair<vector<double>, pair<double, double> > MMP_A2A(string &file_name, double eps, int type, int point_num, int q_num, int q_type){
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
    Base::AABB_tree aabb_tree;
    CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

    Base::Surface_mesh_shortest_path shortest_paths(surface_mesh);
    vector<double> V2V_results = {};
    double query_time = 0.0;

    if (!q_type){
        cout << "V2V not considered" << endl;
    }
    else{
        int percent = 1;
        for (auto i = 0; i < q_num; i++){
            auto s = A2A_query[i].first;
            auto t = A2A_query[i].second;
            if (i > percent * q_num / 10){
                cout << percent++ << "0% queries finished." << endl;
            }

            auto location_s = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(s, aabb_tree, surface_mesh);
            auto location_t = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(t, aabb_tree, surface_mesh);
//            CGAL::Polygon_mesh_processing::barycentric_coordinates()

            auto q_start = chrono::_V2::system_clock::now();
            shortest_paths.add_source_point(location_s.first, location_s.second);
            auto ret_pair = shortest_paths.shortest_distance_to_source_points(location_t.first, location_t.second);
            double mmp_distance = ret_pair.first;
            shortest_paths.remove_all_source_points();
            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);

//            surface_mesh.remove_vertex(svd);
//            surface_mesh.remove_vertex(tvd);

            V2V_results.push_back(mmp_distance);
            query_time += static_cast<double>(q_duration.count());
        }
    }

    return make_pair(V2V_results, make_pair(0.0, query_time));
}

void distanceQueryTest(string &file_name, double eps, int type, int point_num, int q_num, int q_type){
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
//    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsJACM(surface_mesh, eps, gama, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);

    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
    int num_base_graph_vertices = ret_place_points.second;
    cout << "base graph |V| = " << num_base_graph_vertices << endl;

    vector<int> pid_list = {};  //  partition tree will index all vertices
    for (auto i = 0; i < num_base_graph_vertices; i++){
        pid_list.push_back(i);
    }

    WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                               bisector_point_map, point_face_map, point_location_map);

    WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
    double l = 8 / eps + 10;
    tree.constructTree(pid_list, l);

    set<WeightedDistanceOracle::nodePair> node_pairs;
    tree.generateNodePairSet(eps, type, point_num, node_pairs);
    cout << "Node pair set size: " << node_pairs.size() << endl;
    cout << "Node pair percentage: " << fixed << setprecision(3) << 1.0 * node_pairs.size() / num_base_graph_vertices / num_base_graph_vertices << endl;

    auto index_end = chrono::_V2::system_clock::now();
    auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
    double index_time = index_duration.count();
    auto memory_end = Base::physical_memory_used_by_process();
    cout << "Index memory usage: " << (memory_end - memory_begin) / 1000 << " MB" << endl;

    vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());
//    for (auto i = 0; i < leaf_nodes.size(); i++){
//        cout << "leaf i: " << i << "-" << leaf_nodes[i]->center_idx << endl;
//    }
    map<int, int> new_id;
    for (auto i = 0; i < leaf_nodes.size(); i++){
        new_id[i] = i;
    }

    Base::Surface_mesh_shortest_path shortest_paths(surface_mesh);
    vector<double> V2V_results = {};
    double query_time = 0.0;

    for (auto i = 0; i < q_num; i++){
        int s = rand() % leaf_nodes.size();
        int t = rand() % leaf_nodes.size();
        assert(leaf_nodes[s]->center_idx == s);
        assert(leaf_nodes[t]->center_idx == t);

        auto q_start = chrono::_V2::system_clock::now();
        vector<WeightedDistanceOracle::PartitionTreeNode*> As, At;
        tree.getPathToRoot(leaf_nodes[s], As);
        tree.getPathToRoot(leaf_nodes[t], At);
        double oracle_distance = WeightedDistanceOracle::distanceQueryBf(node_pairs, As, At);
        double oracle_distance_eff = WeightedDistanceOracle::distanceQueryEfficient(node_pairs, As, At);
        auto q_end = chrono::_V2::system_clock::now();
        auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
        query_time += static_cast<double>(q_duration.count());
        double err = fabs(oracle_distance_eff - oracle_distance) / oracle_distance;
        if (Base::doubleCmp(err) > 0){
            cout << fixed << setprecision(3) << "bf = " << oracle_distance << "| eff = " << oracle_distance_eff << endl;
        }
    }



//    return make_pair(V2V_results, make_pair(index_time, query_time));
}


int main(int argc, char* argv[]) {
    string file_name;
    int type, point_num, level, q_num, q_type;
    double eps;
    Base::getOpt(argc, argv, file_name, eps, type, point_num, level, q_num, q_type);
    if (q_type){
        cout << "Query type: A2A..." << endl;
    }
    else{
        cout << "Query type: V2V..." << endl;
    }
    cout << fixed << setprecision(3) << "Epsilon: " << eps << endl;
//    highway_test(file_name, eps, point_num);

//    distanceQueryTest(file_name, eps, type, point_num, q_num, q_type);
//    return 0;

    if (!q_type){
        generate_queriesV2V(file_name, eps, type, point_num, q_num);
    }
    else{
        generate_queriesA2A(file_name, eps, type, point_num, q_num);
    }

    auto res_quad = quad_A2A(file_name, eps, point_num, type, level, q_num, q_type);
    cout << "LQT-oracle finished." << endl;

//    auto res_oracle = SE_A2A(file_name, eps, type, point_num, q_num, q_type);
//    cout << "SE-oracle finished." << endl;

//    evaluateBaseGraphEffect(file_name, eps, point_num);
//    kSkipTest(file_name, eps, point_num, type);
//    kSkipGraphTest();

//    auto res_dijk = base_graph_run(q_num, q_type);
//    cout << "Base graph dijkstra finished." << endl;

    auto res_mmp = MMP_A2A(file_name, eps, type, point_num, q_num, q_type);
    cout << "mmp exact finished." << endl;

    cout << "index time cmp: " << endl;
//    cout << fixed << setprecision(3) << "oracle index: " << res_oracle.second.first << " ms." << endl;
    cout << fixed << setprecision(3) << "quad index: " << res_quad.second.first << " ms." << endl;
    cout << "query time cmp: " << endl;
//    cout << fixed << setprecision(3) << "oracle query: " << res_oracle.second.second << " ms." << endl;
//    cout << fixed << setprecision(3) << "dijk query: " << res_dijk.second.second << " ms." << endl;
    cout << fixed << setprecision(3) << "quad query: " << res_quad.second.second << " ms." << endl;
    cout << fixed << setprecision(3) << "MMP query: " << res_mmp.second.second << " ms." << endl;


    cout << "totally " << q_num << " queries" << endl;
    double oracle_tot_err, quad_tot_err,min_err, max_err;
    for (auto i = 0; i < q_num; i++){
//        double relative_error = fabs(res_oracle.first[i] - res_quad.first[i]) / res_oracle.first[i];
//        double relative_error = fabs(res_dijk.first[i] - res_quad.first[i]) / res_dijk.first[i];
        double oracle_relative_error = 0.0;
//        double oracle_relative_error = fabs(res_oracle.first[i] - res_mmp.first[i]) / res_mmp.first[i];
        double quad_relative_error = fabs(res_quad.first[i] - res_mmp.first[i]) / res_mmp.first[i];



        if (Base::doubleCmp(oracle_relative_error - eps) > 0 || Base::doubleCmp(quad_relative_error - eps) > 0){
            cout << "Point 1 = " << A2A_query[i].first << endl;
            cout << "Point 2 = " << A2A_query[i].second << endl;
//            cout << fixed << setprecision(2) << "large error for: " << res_oracle.first[i] << " | " << res_quad.first[i] << " | " <<  << oracle_relative_error << endl;//
//            cout << fixed << setprecision(3) << "oracle error for: " << res_oracle.first[i] << " | " << res_mmp.first[i] << " | " << oracle_relative_error << endl;
            cout << fixed << setprecision(3) << "oracle error for: " << res_quad.first[i] << " | " << res_mmp.first[i] << " | " << quad_relative_error << endl;
            // cout << fixed << setprecision(2) << "large error for: " << res_dijk.first[i] << " | " << res_quad.first[i] << " | " << relative_error << endl;
        }
//        cout << fixed << setprecision(3) << res_oracle.first[i] << " | " << res_quad.first[i] << " | " << relative_error << endl;

        oracle_tot_err += oracle_relative_error;
        quad_tot_err += quad_relative_error;
//        if (Base::doubleCmp(relative_error - min_err) < 0) min_err = relative_error;
//        if (Base::doubleCmp(relative_error - max_err) > 0) max_err = relative_error;
    }
    oracle_tot_err /= q_num;
    quad_tot_err /= q_num;
//    cout << "[QUALITY] min relative = " << fixed << setprecision(3) << min_err << endl;
//    cout << "[QUALITY] max relative = " << fixed << setprecision(3) << max_err << endl;
//    cout << "[QUALITY] oracle avg relative = " << fixed << setprecision(3) << oracle_tot_err << endl;
    cout << "[QUALITY] quad avg relative = " << fixed << setprecision(3) << quad_tot_err << endl;
    return 0;
}
