#include "base.h"
#include "weighted_distance_oracle.h"
#include "highway_well_separated.h"
#include "k_skip.h"

using namespace std;
int K = 3;

void run(string &file_name, double eps, int type, int point_num){
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
//    CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh);
    Base::AABB_tree aabb_tree;
    CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);

    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
    vector<double> face_max_length = {}; // maximum edge length for each face.
    WeightedDistanceOracle::getFaceMaxLength(surface_mesh, face_max_length); // get the maximum edge length for each face.
#ifdef PrintDetails
    for (auto i = 0; i != face_max_length.size(); i++){
        cout << "face " << i << " has maximum edge length " << face_max_length[i] << endl;
    }
#endif
    map<int, vector<int> > edge_bisector_map, bisector_point_map,  face_point_map;
    map<int, int> point_face_map;
    map<int, Base::Point> point_location_map;
    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num, edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
    int num_base_graph_vertices = ret_place_points.second;

    auto ret_base_graph = WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
                                                                     bisector_point_map, point_face_map, point_location_map);

    cout << "[TIME] Base graph construction: " << fixed << setprecision(2) << ret_base_graph << " ms" << endl;
    auto ret_preprocessing_pair = WeightedDistanceOracle::distanceBoundMapPreprocessing(surface_mesh, type,
                                                                                       point_face_map,
                                                                                       face_max_length);
    cout << "[TIME] Distance map preprocessing: " << fixed << setprecision(2) << ret_preprocessing_pair.first << " ms" << endl;
    cout << "[TIME] Bound map preprocessing: " << fixed << setprecision(2) << ret_preprocessing_pair.second << " ms" << endl;
    WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
    auto ret_build_tree = tree.constructTree(surface_mesh.num_vertices());
    cout << "[TIME] Partition tree construction: " << fixed << setprecision(2) << ret_build_tree << " ms" << endl;

    set<WeightedDistanceOracle::nodePair> node_pairs;
    auto ret_node_pair_generate = tree.generateNodePairSet(eps, type, point_num, node_pairs);
    cout << "Node pair set size: " << node_pairs.size() << endl;
    cout << "Node pair percentage: " << fixed << setprecision(2) << 1.0 * node_pairs.size() / surface_mesh.num_vertices() / surface_mesh.num_vertices() << endl;
    cout << "[TIME] Node pair generation: " << fixed << setprecision(2) << ret_node_pair_generate << " ms" << endl;

    vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());

    int cnt = 0;
    double tot_err = 0.0, min_err = 1e20, max_err = -1.0;

    for (auto i = 0; i < 10000; i++){
        auto p1_pair = Base::generateArbitrarySurfacePoint(surface_mesh, aabb_tree);
        auto p2_pair = Base::generateArbitrarySurfacePoint(surface_mesh, aabb_tree);
        auto p1 = p1_pair.first, p2 = p2_pair.first;
        double query_distance = WeightedDistanceOracle::distanceQueryA2A(p1, p2, surface_mesh, tree, node_pairs, aabb_tree);
        if (Base::doubleCmp(query_distance) >= 0){
            cnt++;
//            cout << "p1 = " << p1 << " ; p2 = " << p2 << endl;
//            cout << "query_distance = " << query_distance << endl;
            cout << "oracle distance: " << query_distance << endl;

        }
        else{
            cout << "[INFO] node found in oracles!" << endl;
        }
//        else{
            int s = WeightedDistanceOracle::addVertexEdgeToBaseGraph(surface_mesh, p1, p1_pair.second, face_point_map, point_location_map);
            int t = WeightedDistanceOracle::addVertexEdgeToBaseGraph(surface_mesh, p2, p2_pair.second, face_point_map, point_location_map);
            vector<double> tmp_dall;
            vector<int> tmp_fa;
            Base::dijkstra_SSAD(WeightedDistanceOracle::base_graph, s, tmp_dall, tmp_fa);
            double real_distance = tmp_dall[t];
            cout << "real distance: " << real_distance << endl;

//        }
        if (Base::doubleCmp(query_distance) < 0) continue;
        double relative_error = fabs(query_distance - real_distance) / real_distance;
        tot_err += relative_error;
        if (Base::doubleCmp(relative_error - min_err) < 0) min_err = relative_error;
        if (Base::doubleCmp(relative_error - max_err) > 0) max_err = relative_error;
        if (relative_error > eps){
            cout << "[WARNING] large error found in (" << s << "," << t << "): " << relative_error << " > " << eps << endl;
//            cout << "oracle distance of (" << s << "," << t << "): " << fixed << setprecision(2) << query_distance << endl;
//            cout << "real distance of (" << s << "," << t << "): " << fixed << setprecision(2) << real_distance << endl;
        }
    }
    cout << fixed << setprecision(2) << (double)cnt / 10000 << " percent queries in the node_pair_set." << endl;
    cout << "[QUALITY] min relative = " << fixed << setprecision(2) << min_err << endl;
    cout << "[QUALITY] max relative = " << fixed << setprecision(2) << max_err << endl;
    cout << "[QUALITY] avg relative = " << fixed << setprecision(2) << tot_err / cnt << endl;
    return;

    tot_err = 0.0, min_err = 1e20, max_err = -1.0;
    Base::Surface_mesh_shortest_path shortest_paths(surface_mesh);
    for (auto i = 0; i < 1000; i++){
        int s = rand() % surface_mesh.num_vertices(),
            t = rand() % surface_mesh.num_vertices();
        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
        assert(leaf_nodes[s]->center_idx == s);
        assert(leaf_nodes[t]->center_idx == t);
        vector<WeightedDistanceOracle::PartitionTreeNode*> As, At;
        tree.getPathToRoot(leaf_nodes[s], As);
        tree.getPathToRoot(leaf_nodes[t], At);
        double oracle_distance = WeightedDistanceOracle::distanceQueryBf(node_pairs, As, At);

        auto svd = *(surface_mesh.vertices().begin() + s),
             tvd = *(surface_mesh.vertices().begin() + t);
        shortest_paths.add_source_point(svd);
        auto dis_pair = shortest_paths.shortest_distance_to_source_points(tvd);
        double real_distance = dis_pair.first;
//        double real_distance = 0.0;

        shortest_paths.remove_all_source_points();

        double relative_error = fabs(oracle_distance - real_distance) / real_distance;
        tot_err += relative_error;
        if (Base::doubleCmp(relative_error - min_err) < 0) min_err = relative_error;
        if (Base::doubleCmp(relative_error - max_err) > 0) max_err = relative_error;
        if (relative_error > eps){
            cout << "[WARNING] large error found in (" << s << "," << t << "): " << relative_error << " > " << eps << endl;
            cout << "oracle distance of (" << s << "," << t << "): " << fixed << setprecision(2) << oracle_distance << endl;
            cout << "real distance of (" << s << "," << t << "): " << fixed << setprecision(2) << real_distance << endl;
        }
    }
    cout << "[QUALITY] min relative = " << fixed << setprecision(2) << min_err << endl;
    cout << "[QUALITY] max relative = " << fixed << setprecision(2) << max_err << endl;
    cout << "[QUALITY] avg relative = " << fixed << setprecision(2) << tot_err / 1000 << endl;

#ifdef PrintDetails
    for (auto it = edge_bisector_map.begin(); it != edge_bisector_map.end(); it++){
        cout << "edge id = " << it->first << endl;
        for (auto bid: it->second){
            cout << bid << " ";
        }
        cout << endl;
    }
    for (auto it = bisector_point_map.begin(); it != bisector_point_map.end(); it++){
        cout << "bisector id = " << it->first << endl;
        for (auto pid: it->second){
            cout << pid << " ";
        }
        cout << endl;
    }
    int u = rand() % surface_mesh.num_vertices();
    cout << "u = " << u << endl;
    vector<double> d;
    vector<int> fa;
    Base::dijkstra_SSAD(WeightedDistanceOracle::base_graph, u, d, fa);
    for (auto i = 0; i < 100; i++){
        int v = rand() % surface_mesh.num_vertices();
        if (u == v) v = (v + 1) % surface_mesh.num_vertices();
        cout << "dis " << u << " " << v << " = " << d[v] << " —— ";
        cout << WeightedDistanceOracle::distance_map[u][v] << endl;
//        cout << "fa[" << v << "] = " << fa[v] << endl;
//        cout << "u = " << u << " v = " << v << endl;
//        cout << "distance of " << point_location_map[u] << " and " << point_location_map[v] << " is " << d[v] << endl;
    }
#endif
}

void evaluateBaseGraphEffect(string &file_name, double eps, int point_num){
    srand((int)time(0));

    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
    int num_queries = 1000;
    vector<pair<int, int> > query_pairs = {};   //  pre-computed query pairs
    for (auto i = 0; i != num_queries; i++){
        int s = rand() % surface_mesh.num_vertices(),
                t = rand() % surface_mesh.num_vertices();
        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
        query_pairs.emplace_back(s, t);
    }

    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
    vector<double> face_max_length = {}; // maximum edge length for each face.

    map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
    map<int, int> point_face_map;
    map<int, Base::Point> point_location_map;
    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num,
                                                                             edge_bisector_map, bisector_point_map, point_face_map, point_location_map, face_point_map);
    double time_place_points_fixed = ret_place_points.first;
    int number_placed_points_fixed = ret_place_points.second;

    auto ret_base_graph = WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, number_placed_points_fixed, edge_bisector_map,
                                                                     bisector_point_map, point_face_map,
                                                                     point_location_map);
    double time_base_graph_fixed = ret_base_graph;

    vector<double> fixed_results = {};
    for (auto query: query_pairs){
        vector<double> dall;
        vector<int> fa;
        Base::dijkstra_SSAD(WeightedDistanceOracle::base_graph, query.first, dall, fa);
        fixed_results.push_back(dall[query.second]);
    }

    edge_bisector_map.clear(); bisector_point_map.clear();
    point_face_map.clear(); point_location_map.clear(); face_point_map.clear();
    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
    ret_place_points = WeightedDistanceOracle::placeBisectorPointsJACM(surface_mesh, eps, gama,
                                                                            edge_bisector_map, bisector_point_map,
                                                                            point_face_map, point_location_map, face_point_map);
    double time_place_points_jacm = ret_place_points.first;
    int number_placed_points_jacm = ret_place_points.second;
    ret_base_graph = WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, number_placed_points_jacm, edge_bisector_map,
                                                                     bisector_point_map, point_face_map, point_location_map);
    double time_base_graph_jacm = ret_base_graph;

    cout << "[TIME] Place Steiner points: " << fixed << setprecision(2) << time_place_points_fixed << " vs. " << time_place_points_jacm << " ms" << endl;
    cout << "Number of Steiner points: " << number_placed_points_fixed << " vs. " << number_placed_points_jacm  << " nodes" << endl;
    cout << "[TIME] Base graph construction: " << fixed << setprecision(2) << time_base_graph_fixed << " vs. " << time_base_graph_jacm << " ms" << endl;

    vector<double> jacm_results = {};
    for (auto query: query_pairs){
        vector<double> dall;
        vector<int> fa;
        Base::dijkstra_SSAD(WeightedDistanceOracle::base_graph, query.first, dall, fa);
        jacm_results.push_back(dall[query.second]);
    }

    double tot_err = 0.0;
    double max_err = -1.0;
    double min_err = 1e100;
    for (auto i = 0; i != num_queries; i++){
        double relative_err = fabs(fixed_results[i] - jacm_results[i]) / jacm_results[i];
        cout << fixed_results[i] << " vs. " << jacm_results[i] << " error = " << relative_err << endl;
        tot_err += relative_err;
        if (Base::doubleCmp(relative_err - max_err) > 0) max_err = relative_err;
        if (Base::doubleCmp(relative_err - min_err) < 0) min_err = relative_err;
    }
    cout << "min error = " << min_err << endl;
    cout << "max error = " << max_err << endl;
    cout << "average error = " << tot_err / num_queries << endl;
}

void highway_test(string &file_name, double eps, int point_num){

    srand((int)time(0));

    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
    int num_queries = 1000;
    vector<pair<int, int> > query_pairs = {};   //  pre-computed query pairs
    for (auto i = 0; i != num_queries; i++){
        int s = rand() % surface_mesh.num_vertices(),
                t = rand() % surface_mesh.num_vertices();
        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
        query_pairs.emplace_back(s, t);
    }

    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
    vector<double> face_max_length = {}; // maximum edge length for each face.

    map<int, vector<int> > edge_bisector_map, bisector_point_map, face_point_map;
    map<int, int> point_face_map;
    map<int, Base::Point> point_location_map;

    vector<double> gama = WeightedDistanceOracle::getVertexGamma(surface_mesh, face_weight);
    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsJACM(surface_mesh, eps, gama,
                                                                       edge_bisector_map, bisector_point_map,
                                                                       point_face_map, point_location_map, face_point_map);
    int number_placed_points_jacm = ret_place_points.second;

    Highway::HighwayGraph g;
    g.init(number_placed_points_jacm);
    Highway::constructHighwayGraph(g, surface_mesh, face_weight, number_placed_points_jacm, edge_bisector_map,
                                   bisector_point_map, point_face_map, point_location_map);

    for (auto q: query_pairs){
        auto query_result = g.distanceQuery(q.first, q.second);
        cout << query_result.first << " " << query_result.second << endl;
    }
}

void kSkipTest(){
    ifstream fin("../datasets/grid.txt");
    int V, E;
    fin >> V >> E;
    kSkip::Graph g;
    g.init(V);
    int u, v; double w;
    for (auto i = 0; i < E; i++){
        fin >> u >> v >> w;
        g.addEdge(u, v, w);
    }
    auto ret = kSkip::adaptiveSampling(g);
    for (auto &vid: ret){
        cout << vid << " ";
    }
    cout << endl;
}

int main(int argc, char* argv[]) {
    string file_name;
    int type, point_num;
    double eps;
    Base::getOpt(argc, argv, file_name, eps, type, point_num);
//    highway_test(file_name, eps, point_num);
//    run(file_name, eps, type, point_num);
//    evaluateBaseGraphEffect(file_name, eps, point_num);
    kSkipTest();
    return 0;
}
