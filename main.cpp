#include "base.h"
#include "weighted_distance_oracle.h"

using namespace std;

void run(string &file_name, double eps, int type, int point_num){
    srand((int)time(0));
    Base::Mesh surface_mesh;
    ifstream fin(file_name);
    fin >> surface_mesh;
    vector<double> face_weight(surface_mesh.num_faces(), 1.0); // face weight for each face.
    vector<double> face_max_length = {}; // maximum edge length for each face.
    WeightedDistanceOracle::getFaceMaxLength(surface_mesh, face_max_length); // get the maximum edge length for each face.
#ifdef PrintDetails
    for (auto i = 0; i != face_max_length.size(); i++){
        cout << "face " << i << " has maximum edge length " << face_max_length[i] << endl;
    }
#endif
    map<int, vector<int> > edge_bisector_map, bisector_point_map;
    map<int, int> point_face_map;
    map<int, Base::Point> point_location_map;
    auto ret_place_points = WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num, edge_bisector_map, bisector_point_map, point_face_map, point_location_map);
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
    auto ret_build_tree = tree.construct_tree(surface_mesh.num_vertices());
    cout << "[TIME] Partition tree construction: " << fixed << setprecision(2) << ret_build_tree << " ms" << endl;

    set<WeightedDistanceOracle::nodePair> node_pairs;
    auto ret_node_pair_generate = tree.generate_node_pair_set(eps, type, point_num, node_pairs);
    cout << "Node pair set size: " << node_pairs.size() << endl;
    cout << "Node pair percentage: " << fixed << setprecision(2) << 1.0 * node_pairs.size() / surface_mesh.num_vertices() / surface_mesh.num_vertices() << endl;
    cout << "[TIME] Node pair generation: " << fixed << setprecision(2) << ret_node_pair_generate << " ms" << endl;

    vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());
    double tot_err = 0.0, min_err = 1e20, max_err = -1.0;
    Base::Surface_mesh_shortest_path shortest_paths(surface_mesh);
    for (auto i = 0; i < 1000; i++){
        int s = rand() % surface_mesh.num_vertices(),
            t = rand() % surface_mesh.num_vertices();
        if (s == t) t = (t + 1) % surface_mesh.num_vertices();
        assert(leaf_nodes[s]->center_idx == s);
        assert(leaf_nodes[t]->center_idx == t);
        vector<WeightedDistanceOracle::PartitionTreeNode*> As, At;
        tree.get_path_to_root(leaf_nodes[s], As);
        tree.get_path_to_root(leaf_nodes[t], At);
        double oracle_distance = WeightedDistanceOracle::distance_query_bf(node_pairs, As, At);

        auto svd = *(surface_mesh.vertices().begin() + s),
             tvd = *(surface_mesh.vertices().begin() + t);
        shortest_paths.add_source_point(svd);
        auto dis_pair = shortest_paths.shortest_distance_to_source_points(tvd);
        double real_distance = dis_pair.first;
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

int main(int argc, char* argv[]) {
    string file_name;
    int type, point_num;
    double eps;
    Base::getOpt(argc, argv, file_name, eps, type, point_num);
    run(file_name, eps, type, point_num);
    return 0;
}
