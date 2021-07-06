#include "base.h"
#include "weighted_distance_oracle.h"

using namespace std;

void run(string &file_name, double eps, int type, int point_num){
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
    WeightedDistanceOracle::placeBisectorPointsFixed(surface_mesh, point_num, edge_bisector_map, bisector_point_map, point_face_map, point_location_map);
}

int main(int argc, char* argv[]) {
    string file_name;
    int type, point_num;
    double eps;
    Base::getOpt(argc, argv, file_name, eps, type, point_num);
    run(file_name, eps, type, point_num);
    return 0;
}
