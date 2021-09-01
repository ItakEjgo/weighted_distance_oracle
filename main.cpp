#include "base.h"
#include "weighted_distance_oracle.h"
#include "highway_well_separated.h"
#include "k_skip.h"
#include "quad.h"

using namespace std;
int K;
vector<pair<int, int> > V2V_query;

pair<vector<double>, pair<double, double> > run(string &file_name, double eps, int type, int point_num, int q_num){
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
//    CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh);
    Base::AABB_tree aabb_tree;
    CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

    for (auto i = 0; i < q_num; i++){
        int s = rand() % surface_mesh.num_vertices(),
                t = rand() % surface_mesh.num_vertices();
        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
        V2V_query.push_back(make_pair(s, t));
    }

    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
    vector<double> face_max_length = {}; // maximum edge length for each face.
    WeightedDistanceOracle::getFaceMaxLength(surface_mesh, face_max_length); // get the maximum edge length for each face.

    map<int, vector<int> > edge_bisector_map, bisector_point_map,  face_point_map;
    map<int, int> point_face_map;
    map<int, Base::Point> point_location_map;
    vector<int> pid_list = {};
    for (auto i = 0; i < surface_mesh.num_vertices(); i++){
        pid_list.push_back(i);
    }

    auto index_start = chrono::_V2::system_clock::now();

    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
    int num_base_graph_vertices = ret_place_points.second;
    auto ret_base_graph = WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                                     bisector_point_map, point_face_map, point_location_map);

    auto ret_preprocessing_pair = WeightedDistanceOracle::distanceBoundMapPreprocessing(surface_mesh, type,
                                                                                       point_face_map,
                                                                                       face_max_length,
                                                                                       point_location_map,
                                                                                       pid_list);
    WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
    auto ret_build_tree = tree.constructTree(pid_list);

    set<WeightedDistanceOracle::nodePair> node_pairs;
    auto ret_node_pair_generate = tree.generateNodePairSet(eps, type, point_num, node_pairs);
    cout << "Node pair set size: " << node_pairs.size() << endl;
    cout << "Node pair percentage: " << fixed << setprecision(2) << 1.0 * node_pairs.size() / surface_mesh.num_vertices() / surface_mesh.num_vertices() << endl;

    auto index_end = chrono::_V2::system_clock::now();
    auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
    double index_time = index_duration.count();


    vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

//    for (auto i = 0; i < 10000; i++){
//        auto p1_pair = Base::generateArbitrarySurfacePoint(surface_mesh, aabb_tree);
//        auto p2_pair = Base::generateArbitrarySurfacePoint(surface_mesh, aabb_tree);
//        auto p1 = p1_pair.first, p2 = p2_pair.first;
//        double query_distance = WeightedDistanceOracle::distanceQueryA2A(p1, p2, surface_mesh, tree, node_pairs, aabb_tree);
//        if (Base::doubleCmp(query_distance) >= 0){
//            cnt++;
////            cout << "p1 = " << p1 << " ; p2 = " << p2 << endl;
////            cout << "query_distance = " << query_distance << endl;
//            cout << "oracle distance: " << query_distance << endl;
//
//        }
//        else{
//            cout << "[INFO] node found in oracles!" << endl;
//        }
////        else{
//            int s = WeightedDistanceOracle::addVertexEdgeToBaseGraph(surface_mesh, p1, p1_pair.second, face_point_map, point_location_map);
//            int t = WeightedDistanceOracle::addVertexEdgeToBaseGraph(surface_mesh, p2, p2_pair.second, face_point_map, point_location_map);
//            vector<double> tmp_dall;
//            vector<int> tmp_fa;
//            Base::dijkstra_SSAD(WeightedDistanceOracle::base_graph, s, tmp_dall, tmp_fa);
//            double real_distance = tmp_dall[t];
//            cout << "real distance: " << real_distance << endl;
//
////        }
//        if (Base::doubleCmp(query_distance) < 0) continue;
//        double relative_error = fabs(query_distance - real_distance) / real_distance;
//        tot_err += relative_error;
//        if (Base::doubleCmp(relative_error - min_err) < 0) min_err = relative_error;
//        if (Base::doubleCmp(relative_error - max_err) > 0) max_err = relative_error;
//        if (relative_error > eps){
//            cout << "[WARNING] large error found in (" << s << "," << t << "): " << relative_error << " > " << eps << endl;
////            cout << "oracle distance of (" << s << "," << t << "): " << fixed << setprecision(2) << query_distance << endl;
////            cout << "real distance of (" << s << "," << t << "): " << fixed << setprecision(2) << real_distance << endl;
//        }
//    }
//    cout << fixed << setprecision(2) << (double)cnt / 10000 << " percent queries in the node_pair_set." << endl;
//    cout << "[QUALITY] min relative = " << fixed << setprecision(2) << min_err << endl;
//    cout << "[QUALITY] max relative = " << fixed << setprecision(2) << max_err << endl;
//    cout << "[QUALITY] avg relative = " << fixed << setprecision(2) << tot_err / cnt << endl;
//    return;

//    tot_err = 0.0, min_err = 1e20, max_err = -1.0;
    Base::Surface_mesh_shortest_path shortest_paths(surface_mesh);
    vector<double> V2V_results = {};
    double query_time = 0.0;
    for (auto i = 0; i < q_num; i++){
        int s = V2V_query[i].first;
        int t = V2V_query[i].second;
        assert(leaf_nodes[s]->center_idx == s);
        assert(leaf_nodes[t]->center_idx == t);
//        cout << s << "->" << t << endl;

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
    return make_pair(V2V_results, make_pair(index_time, query_time));
}
//
//void evaluateBaseGraphEffect(string &file_name, double eps, int point_num){
//    srand((int)time(0));
//
//    Base::Mesh surface_mesh;
//    ifstream fin(file_name);
//    fin >> surface_mesh;
//    int num_queries = 1000;
//    vector<pair<int, int> > query_pairs = {};   //  pre-computed query pairs
//    for (auto i = 0; i != num_queries; i++){
//        int s = rand() % surface_mesh.num_vertices(),
//                t = rand() % surface_mesh.num_vertices();
//        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
//        query_pairs.emplace_back(s, t);
//    }
//
//    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
//    vector<double> face_max_length = {}; // maximum edge length for each face.
//
//    map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
//    map<int, int> point_face_map;
//    map<int, Base::Point> point_location_map;
//    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num,
//                                                                             edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
//    double time_place_points_fixed = ret_place_points.first;
//    int number_placed_points_fixed = ret_place_points.second;
//
//    auto ret_base_graph = WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, number_placed_points_fixed, edge_bisector_map,
//                                                                     bisector_point_map, point_face_map,
//                                                                     point_location_map);
//    double time_base_graph_fixed = ret_base_graph;
//
//    vector<double> fixed_results = {};
//    for (auto query: query_pairs){
//        vector<double> dall;
//        vector<int> fa;
//        Base::dijkstra_SSAD(WeightedDistanceOracle::base_graph, query.first, dall, fa);
//        fixed_results.push_back(dall[query.second]);
//    }
//
//    edge_bisector_map.clear(); bisector_point_map.clear();
//    point_face_map.clear(); point_location_map.clear(); face_point_map.clear();
//    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    ret_place_points = WeightedDistanceOracle::placeBisectorPointsJACM(surface_mesh, eps, gama,
//                                                                            edge_bisector_map, bisector_point_map,
//                                                                            point_face_map, point_location_map, face_point_map);
//    double time_place_points_jacm = ret_place_points.first;
//    int number_placed_points_jacm = ret_place_points.second;
//    ret_base_graph = WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, number_placed_points_jacm, edge_bisector_map,
//                                                                     bisector_point_map, point_face_map, point_location_map);
//    double time_base_graph_jacm = ret_base_graph;
//
//    cout << "[TIME] Place Steiner points: " << fixed << setprecision(2) << time_place_points_fixed << " vs. " << time_place_points_jacm << " ms" << endl;
//    cout << "Number of Steiner points: " << number_placed_points_fixed << " vs. " << number_placed_points_jacm  << " nodes" << endl;
//    cout << "[TIME] Base graph construction: " << fixed << setprecision(2) << time_base_graph_fixed << " vs. " << time_base_graph_jacm << " ms" << endl;
//
//    vector<double> jacm_results = {};
//    for (auto query: query_pairs){
//        vector<double> dall;
//        vector<int> fa;
//        Base::dijkstra_SSAD(WeightedDistanceOracle::base_graph, query.first, dall, fa);
//        jacm_results.push_back(dall[query.second]);
//    }
//
//    double tot_err = 0.0;
//    double max_err = -1.0;
//    double min_err = 1e100;
//    for (auto i = 0; i != num_queries; i++){
//        double relative_err = fabs(fixed_results[i] - jacm_results[i]) / jacm_results[i];
//        cout << fixed_results[i] << " vs. " << jacm_results[i] << " error = " << relative_err << endl;
//        tot_err += relative_err;
//        if (Base::doubleCmp(relative_err - max_err) > 0) max_err = relative_err;
//        if (Base::doubleCmp(relative_err - min_err) < 0) min_err = relative_err;
//    }
//    cout << "min error = " << min_err << endl;
//    cout << "max error = " << max_err << endl;
//    cout << "average error = " << tot_err / num_queries << endl;
//}
//
//void highway_test(string &file_name, double eps, int point_num){
//
//    srand((int)time(0));
//
//    Base::Mesh surface_mesh;
//    ifstream fin(file_name);
//    fin >> surface_mesh;
//    int num_queries = 1000;
//    vector<pair<int, int> > query_pairs = {};   //  pre-computed query pairs
//    for (auto i = 0; i != num_queries; i++){
//        int s = rand() % surface_mesh.num_vertices(),
//                t = rand() % surface_mesh.num_vertices();
//        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
//        query_pairs.emplace_back(s, t);
//    }
//
//    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
//    vector<double> face_max_length = {}; // maximum edge length for each face.
//
//    map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
//    map<int, int> point_face_map;
//    map<int, Base::Point> point_location_map;
//
//    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsJACM(surface_mesh, eps, gama,
//                                                                       edge_bisector_map, bisector_point_map,
//                                                                       point_face_map, point_location_map, face_point_map);
//    int number_placed_points_jacm = ret_place_points.second;
//
//    Highway::HighwayGraph g;
//    g.init(number_placed_points_jacm);
//    Highway::constructHighwayGraph(g, surface_mesh, face_weight, number_placed_points_jacm, edge_bisector_map,
//                                   bisector_point_map, point_face_map, point_location_map);
//
//    for (auto q: query_pairs){
//        auto query_result = g.distanceQuery(q.first, q.second);
//        cout << query_result.first << " " << query_result.second << endl;
//    }
//}
//
//void kSkipTest(string &file_name, double eps, int point_num, int type){
//    K = type;
//    cout << "k in k-skip graph is: " << K << endl;
//    srand((int)time(0));
//    Base::Mesh surface_mesh;
//    ifstream fin(file_name);
//    fin >> surface_mesh;
//    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
//    vector<double> face_max_length = {}; // maximum edge length for each face.
//
//    map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
//    map<int, int> point_face_map;
//    map<int, Base::Point> point_location_map;
//
//    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
//    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsJACM(surface_mesh, eps, gama,
//                                                                            edge_bisector_map, bisector_point_map,
//                                                                            point_face_map, point_location_map, face_point_map);
//    int number_placed_points_jacm = ret_place_points.second;
//    cout << "Place Steiner points finished." << endl;
//
//    kSkip::Graph g;
//    g.init(number_placed_points_jacm);
//    kSkip::constructGraph(g, surface_mesh, face_weight, number_placed_points_jacm, edge_bisector_map,
//                                   bisector_point_map, point_face_map, point_location_map);
//
//    cout << "Construct graph finished." << endl;
//
//    auto k_cover_V = kSkip::adaptiveSampling(g);
//    auto k_cover_vertex_id = kSkip::generateKCoverVertexId(k_cover_V);
//    auto k_skip_graph = kSkip::computeKSkipGraph(g, k_cover_V, k_cover_vertex_id);
//    cout << "V = " << g.num_V << " k-skip V = " << k_skip_graph.num_V << endl;
//    cout << "E = " << g.num_E << " k-skip E = " << k_skip_graph.num_E << endl;
//    cout << "compressed V percent: " << static_cast<double>(k_skip_graph.num_V) / g.num_V << endl;
//    cout << "compressed E percent: " << static_cast<double>(k_skip_graph.num_E) / g.num_E << endl;
////    cout << "k-cover-vertices: ";
////    for (auto vid: k_cover_V){
////        cout << vid << " ";
////    }
////    cout << endl;
//
//    cout << "level 2 start construction:" << endl;
//    auto k2_cover_V = kSkip::adaptiveSampling(k_skip_graph);
//    auto k2_cover_vertex_id = kSkip::generateKCoverVertexId(k2_cover_V);
//    auto k2_skip_graph = kSkip::computeKSkipGraph(k_skip_graph, k2_cover_V, k2_cover_vertex_id);
//    cout << "V = " << g.num_V << " k2-skip V = " << k2_skip_graph.num_V << endl;
//    cout << "E = " << g.num_E << " k2-skip E = " << k2_skip_graph.num_E << endl;
//    cout << "compressed V percent: " << static_cast<double>(k2_skip_graph.num_V) / g.num_V << endl;
//    cout << "compressed E percent: " << static_cast<double>(k2_skip_graph.num_E) / g.num_E << endl;
//
//    int q_num = 10;
//    double tot = 0, tot_dijk = 0;
////    vector<int> vec(k_cover_V.begin(), k_cover_V.end());
//    vector<int> vec(k2_cover_V.begin(), k2_cover_V.end());
//
//    for (auto i = 0; i < q_num; i++){
////        int s = rand() % surface_mesh.num_vertices();
////        int t = rand() % surface_mesh.num_vertices();
////        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
//        int s = rand() % vec.size();
//        int t = rand() % vec.size();
//        if (s == t) t = (t + 1) % vec.size();
//        s = vec[s]; t = vec[t];
//        auto ret = kSkip::queryKSkipGraph(k_skip_graph, k2_skip_graph, k2_cover_V, k2_cover_vertex_id, s, t);
//        auto ret_ori_g = kSkip::dijkstra(k_skip_graph, s, t);
//        if (Base::doubleCmp(ret_ori_g.first - ret.first)){
//            cout << "kSkip dis(" << s << "," << t << "): ";
//            cout << fixed << setprecision(3) << ret.first << ", query time = " << ret.second << " ms" << endl;
//            cout << "exact dis(" << s << "," << t << "): ";
//            cout << fixed << setprecision(3) << ret_ori_g.first << ", query time = " << ret_ori_g.second << " ms" << endl;
//        }
//        tot += ret.second;
//        tot_dijk += ret_ori_g.second;
//    }
//    cout << "AVG k-skip query time = " << tot / 1000 << endl;
//    cout << "AVG dijkstra query time = " << tot_dijk / 1000 << endl;
//}
//
//void kSkipGraphTest(){
//    ifstream fin("../datasets/challenge9/USA-road-d.BAY.gr");
//    kSkip::Graph g;
//    char c;
//    K = 10;
//    int V, E;
//    while (fin >> c){
//        if (c == 'c'){
//            char s[100];
//            fin.getline(s, 100);
//        }
//        else if (c == 'p'){
//            string sp;
//            fin >> sp >> V >> E;
//            g.init(V);
//        }
//        else if (c == 'a'){
//            int u, v; double w;
//            fin >> u >> v >> w;
//            g.addEdge(u - 1, v - 1, w);
//            g.addEdge(v - 1, u - 1, w);
//        }
//    }
//
//    cout << "Construct graph finished." << endl;
//
//    auto k_cover_V = kSkip::adaptiveSampling(g);
//    auto k_cover_vertex_id = kSkip::generateKCoverVertexId(k_cover_V);
//    auto k_skip_graph = kSkip::computeKSkipGraph(g, k_cover_V, k_cover_vertex_id);
//    cout << "V = " << g.num_V << " k-skip V = " << k_skip_graph.num_V << endl;
//    cout << "E = " << g.num_E << " k-skip E = " << k_skip_graph.num_E << endl;
//    cout << "compressed V percent: " << static_cast<double>(k_skip_graph.num_V) / g.num_V << endl;
//    cout << "compressed E percent: " << static_cast<double>(k_skip_graph.num_E) / g.num_E << endl;
//
//
//
//
//    int q_num = 1000;
//    double tot = 0, tot_dijk = 0;
//    vector<int> vec(k_cover_V.begin(), k_cover_V.end());
////    vector<int> vec(k2_cover_V.begin(), k2_cover_V.end());
//
//    for (auto i = 0; i < q_num; i++){
////        int s = rand() % V;
////        int t = rand() % V;
////        if (s == t) t = (t + 1) % V;
//        int s = rand() % vec.size();
//        int t = rand() % vec.size();
//        if (s == t) t = (t + 1) % vec.size();
//        s = vec[s]; t = vec[t];
//
//        auto ret = kSkip::queryKSkipGraph(g, k_skip_graph, k_cover_V, k_cover_vertex_id, s, t);
//        auto ret_ori_g = kSkip::dijkstra(g, s, t);
//        if (Base::doubleCmp(ret_ori_g.first - ret.first)){
//            cout << "kSkip dis(" << s << "," << t << "): ";
//            cout << fixed << setprecision(3) << ret.first << ", query time = " << ret.second << " ms" << endl;
//            cout << "exact dis(" << s << "," << t << "): ";
//            cout << fixed << setprecision(3) << ret_ori_g.first << ", query time = " << ret_ori_g.second << " ms" << endl;
//        }
//        tot += ret.second;
//        tot_dijk += ret_ori_g.second;
//    }
//    cout << "AVG query time = " << tot / 1000 << endl;
//    cout << "AVG dijkstra query time = " << tot_dijk / 1000 << endl;
//}

pair<vector<double>, pair<double, double> > quadTest(string &file_name, double eps, int point_num, int type, int level, int q_num) {
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
    vector<double> face_max_length = {}; // maximum edge length for each face.
    WeightedDistanceOracle::getFaceMaxLength(surface_mesh, face_max_length); // get the maximum edge length for each face.
    double max_len = -1.0;
    for (auto val: face_max_length){
        max_len = max(max_len, val);
    }
//    cout << "max edge length: " << max_len << endl;


    WeightedDistanceOracle::getFaceMaxLength(surface_mesh, face_max_length); // get the maximum edge length for each face.
    map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
    map<int, int> point_face_map;
    map<int, Base::Point> point_location_map;

    auto index_start = chrono::_V2::system_clock::now();

    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
    Quad::quadTree quad_tree(surface_mesh, face_point_map);
    cout << "quad tree level = " << level << endl;
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

    auto ret_base_graph = WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                                     bisector_point_map, point_face_map, point_location_map);


    auto ret_preprocessing_pair = WeightedDistanceOracle::distanceBoundMapPreprocessing(surface_mesh, type,
                                                                                       point_face_map,
                                                                                       face_max_length,
                                                                                       point_location_map,
                                                                                       pid_list);
    WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
    auto ret_build_tree = tree.constructTree(pid_list);

    set<WeightedDistanceOracle::nodePair> node_pairs;
    auto ret_node_pair_generate = tree.generateNodePairSet(eps, type, point_num, node_pairs);
    cout << "Node pair set size: " << node_pairs.size() << endl;
    cout << "Node pair percentage: " << fixed << setprecision(2) << 1.0 * node_pairs.size() / surface_mesh.num_vertices() / surface_mesh.num_vertices() << endl;
//    vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());
    map<int, int> new_id;
    auto spanner = Quad::generateSpanner(pid_list, node_pairs, new_id);
    auto index_end = chrono::_V2::system_clock::now();
    auto index_duration = chrono::duration_cast<chrono::milliseconds>(index_end - index_start);
    double index_time = index_duration.count();

    vector<double> V2V_results = {};
    double query_time = 0.0;
    for (auto i = 0; i < q_num; i++){
        int s = V2V_query[i].first;
        int t = V2V_query[i].second;
//        int s = rand() % surface_mesh.num_vertices();
//        int t = rand() % surface_mesh.num_vertices();
        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
        auto q_start = chrono::_V2::system_clock::now();
        double spanner_distance = Quad::querySpanner(surface_mesh, spanner, s, t, quad_tree, new_id);
        auto q_end = chrono::_V2::system_clock::now();
        auto q_duration = chrono::duration_cast<chrono::milliseconds>(q_end - q_start);
        query_time += static_cast<double>(q_duration.count());

        V2V_results.push_back(spanner_distance);

    }
    return make_pair(V2V_results, make_pair(index_time, query_time));
}

int main(int argc, char* argv[]) {
    string file_name;
    int type, point_num, level, q_num;
    double eps;
    Base::getOpt(argc, argv, file_name, eps, type, point_num, level, q_num);
//    highway_test(file_name, eps, point_num);
    auto res_oracle = run(file_name, eps, type, point_num, q_num);
//    evaluateBaseGraphEffect(file_name, eps, point_num);
//    kSkipTest(file_name, eps, point_num, type);
//    kSkipGraphTest();
    auto res_quad = quadTest(file_name, eps, point_num, type, level, q_num);

    cout << "index time cmp: " << endl;
    cout << fixed << setprecision(3) << "oracle index: " << res_oracle.second.first << " ms." << endl;
    cout << fixed << setprecision(3) << "quad index: " << res_quad.second.first << " ms." << endl;
    cout << "query time cmp: " << endl;
    cout << fixed << setprecision(3) << "oracle query: " << res_oracle.second.second << " ms." << endl;
    cout << fixed << setprecision(3) << "quad query: " << res_quad.second.second << " ms." << endl;
    cout << "totally " << q_num << " queries" << endl;
    double tot_err, min_err, max_err;
    for (auto i = 0; i < q_num; i++){
        double relative_error = fabs(res_oracle.first[i] - res_quad.first[i]) / res_oracle.first[i];
        if (Base::doubleCmp(relative_error - 0.5) > 0){
            cout << fixed << setprecision(2) << "large error for: " << res_oracle.first[i] << " | " << res_quad.first[i] << " | " << relative_error << endl;
        }

        tot_err += relative_error;
        if (Base::doubleCmp(relative_error - min_err) < 0) min_err = relative_error;
        if (Base::doubleCmp(relative_error - max_err) > 0) max_err = relative_error;
    }
    tot_err /= q_num;
    cout << "[QUALITY] min relative = " << fixed << setprecision(3) << min_err << endl;
    cout << "[QUALITY] max relative = " << fixed << setprecision(3) << max_err << endl;
    cout << "[QUALITY] avg relative = " << fixed << setprecision(3) << tot_err << endl;
    return 0;
}
