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
    auto ret_preprocessing = WeightedDistanceOracle::distanceBoundMapPreprocessing(surface_mesh, 0, point_face_map, face_max_length);
//    auto ret_preprocessing = WeightedDistanceOracle::distanceBoundMapPreprocessing(surface_mesh, 1, point_face_map, face_max_length);

    WeightedDistanceOracle::PartitionTree tree = WeightedDistanceOracle::PartitionTree();
    auto ret_build_tree = tree.construct_tree(surface_mesh.num_vertices());

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
