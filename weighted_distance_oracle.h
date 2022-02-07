//
// Created by huang on 2021/7/6.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H
#define WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H

#include "base.h"
#include "k_skip.h"

namespace WeightedDistanceOracle {

    using namespace std;

    Base::Graph base_graph;
    map<int, vector<double> > distance_map, bound_map;

    // return the vertex id of inserted point
    int addVertexEdgeToBaseGraph(Base::Mesh &mesh, Base::Point &p, int fid, map<int, vector<int> > &face_point_map,
                                  map<int, Base::Point> &point_location_map){
        auto fd = *(mesh.faces().begin() + fid);
        int V = boost::num_vertices(base_graph);
//        cout << "old V size = " << V << endl;
        for (auto pid: face_point_map[fid]){
            double w = sqrt(CGAL::squared_distance(p, point_location_map[pid]));
//            cout << "u,v,w = " << V << " " << pid << " " << w << endl;
            boost::add_edge(V, pid, w, base_graph);
            boost::add_edge(pid, V, w, base_graph);
        }
        return V;
//        cout << "new V size = " << boost::num_vertices(base_graph) << endl;
    }

    double getFaceMaxLength(Base::Mesh &mesh, vector<double> &face_max_length) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        face_max_length.resize(mesh.num_faces(), 0.0);
        double max_length = -1.0, e_length;
        for (auto fd: mesh.faces()) {
            max_length = -1.0;
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))) {
                e_length = sqrt(
                        CGAL::squared_distance(mesh.points()[mesh.source(hed)], mesh.points()[mesh.target(hed)]));
                max_length = Base::doubleCmp(e_length - max_length) > 0 ? e_length : max_length;
            }
            face_max_length[fd.idx()] = max_length;
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    double getFaceMaxMinLength(Base::Mesh &mesh, vector<double> &face_max_length, vector<double> &face_min_length) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        face_max_length.resize(mesh.num_faces(), 0.0);
        face_min_length.resize(mesh.num_faces(), 0.0);

        double max_length = -1.0, min_length = 1e60, e_length;
        for (auto fd: mesh.faces()) {
            max_length = -1.0;
            min_length = 1e60;
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))) {
                e_length = sqrt(
                        CGAL::squared_distance(mesh.points()[mesh.source(hed)], mesh.points()[mesh.target(hed)]));
                max_length = Base::doubleCmp(e_length - max_length) > 0 ? e_length : max_length;
                min_length = Base::doubleCmp(e_length - min_length) < 0 ? e_length : min_length;
            }
            face_max_length[fd.idx()] = max_length;
            face_min_length[fd.idx()] = min_length;
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    vector<double> getVertexGamma(Base::Mesh &mesh, vector<double> &face_weights){
        map<int, vector<int> > vertex_face_map = {};
        vector<double> gama(mesh.num_vertices());
        double dis, w_min, w_max;
        for (auto fd: mesh.faces()) {
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))) {
                int source_idx = mesh.source(hed);
                assert(mesh.face(hed) != mesh.null_face());
                vertex_face_map[source_idx].push_back(mesh.face(hed).idx());
            }
        }

        for (auto vd: mesh.vertices()){
            int vid = vd.idx();
            Base::Point v(mesh.points()[vd]);
            dis = 1e100;
            w_max = -1.0;
            w_min = 1e100;
            for (auto fid: vertex_face_map[vid]){
                if (Base::doubleCmp(face_weights[fid] - w_max) > 0) w_max = face_weights[fid];
                if (Base::doubleCmp(face_weights[fid] - w_min) < 0) w_min = face_weights[fid];
                double tmp_face_dis = -1.0;
                for (auto &hed: mesh.halfedges_around_face(mesh.halfedge(*(mesh.faces_begin() + fid)))){
                    auto source = mesh.source(hed), target = mesh.target(hed);
                    Base::Point p(mesh.points()[source]), q(mesh.points()[target]);
                    double t = CGAL::squared_distance(v, Base::Segment(p, q));
                    if (Base::doubleCmp(t - tmp_face_dis) > 0){
                        tmp_face_dis = t;
                    }
                }
                if (Base::doubleCmp(tmp_face_dis - dis) < 0){
                    dis = tmp_face_dis;
                }
            }
            gama[vid] = w_min * sqrt(dis) / w_max / 7;
        }
        return gama;
    }

    pair<double, int> //    int num_base_graph_vertices = ret_place_points.second;
//
//    auto ret_base_graph = WeightedDistanceOracle::constructBaseGraph(surface_mesh, face_weight, num_base_graph_vertices, edge_bisector_map,
//                                                                     bisector_point_map, point_face_map, point_location_map);
//
//    cout << "[TIME] Base graph construction: " << fixed << setprecision(2) << ret_base_graph << " ms" << endl;
    placeSteinerPointsJACM(Base::Mesh &mesh, const double &eps,
                           vector<double> &gama,
                           map<int, vector<int>> &edge_bisector_map,
                           map<int, vector<int>> &bisector_point_map,
                           map<int, int> &point_face_map,
                           map<int, Base::Point> &point_location_map,
                           map<int, vector<int>> &face_point_map){
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        int num_vertices = static_cast<int>(mesh.num_vertices());
        int bisector_num = 0;
        for (auto  fd: mesh.faces()){
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))){
                int eid = static_cast<int>(mesh.edge(hed).idx()),
                    eid2 = static_cast<int>(mesh.edge(mesh.prev(hed)).idx());
                vector<Base::Point> p(3);
                vector<int> pid(3);
                for (auto i = 0; i != 3; hed = mesh.next(hed), i++){
                    p[i] = mesh.points()[mesh.source(hed)];
                    pid[i] = static_cast<int>(mesh.source(hed).idx());
                    if (point_location_map.find(pid[i]) == point_location_map.end()){
                        point_location_map[pid[i]] = p[i];
                    }
                    face_point_map[static_cast<int>(fd.idx())].push_back(pid[i]);
                }
                vector<int> cur_bisector = {};
                cur_bisector.push_back(pid[0]);

                double len1 = sqrt(CGAL::squared_distance(p[1], p[0])),
                        len2 = sqrt(CGAL::squared_distance(p[2], p[0]));
                Base::Point p_end(p[1] + Base::Vector(p[2] - p[1]) * len1 / (len1 + len2));
                Base::Vector vec_bisector(p_end - p[0]);
                Base::Point aux1 = p[0] + Base::Vector(p[1] - p[0]) * gama[pid[0]] / len1;
                Base::Point aux2 = p[0] + Base::Vector(p[2] - p[0]) * gama[pid[0]] / len2;
                Base::Point bisector_p0 = aux1 + 0.5 * Base::Vector(aux2 - aux1);
                double angle = Base::PI * CGAL::approximate_angle(p[1], p[0], p[2]) / 180;
                double sin_val = sin(angle * 0.5);
                double limit_distance = sqrt(CGAL::squared_distance(p[0], p_end));
                double cur_distance = sqrt(CGAL::squared_distance(p[0], bisector_p0));

                while (Base::doubleCmp(cur_distance - limit_distance) < 0) {
                    Base::Point bisector_p = p[0] + cur_distance / limit_distance * vec_bisector;
                    point_location_map[num_vertices] = bisector_p;
                    point_face_map[num_vertices] = static_cast<int>(fd.idx());
                    face_point_map[static_cast<int>(fd.idx())].push_back(num_vertices);
                    cur_bisector.push_back(num_vertices++);
                    double distance_delta = sin_val * sqrt(0.5 * eps) * cur_distance;
                    cur_distance += distance_delta;
                }
                bisector_point_map[bisector_num] = cur_bisector;
                edge_bisector_map[eid].push_back(bisector_num);
                edge_bisector_map[eid2].push_back(bisector_num++);
            }
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return make_pair(static_cast<double>(duration.count()), num_vertices);
    }

    pair<double, int> placeSteinerPointsFixed(Base::Mesh &mesh, const int &point_num,
                                              map<int, vector<int>> &edge_bisector_map,
                                              map<int, vector<int>> &bisector_point_map,
                                              map<int, int> &point_face_map,
                                              map<int, Base::Point> &point_location_map,
                                              map<int, vector<int>> &face_point_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        int num_vertices = static_cast<int>(mesh.num_vertices());
        int bisector_num = 0;
        for (auto fd: mesh.faces()) {
            for (auto hed:mesh.halfedges_around_face(mesh.halfedge(fd))) {
                int eid = static_cast<int>(mesh.edge(hed).idx()),
                    eid2 = static_cast<int>(mesh.edge(mesh.prev(hed)).idx());
                vector<Base::Point> p(3);
                vector<int> pid(3);
                for (auto i = 0; i != 3; hed = mesh.next(hed), i++) {
                    p[i] = mesh.points()[mesh.source(hed)];
                    pid[i] = static_cast<int>(mesh.source(hed).idx());
                    if (point_location_map.find(pid[i]) == point_location_map.end()) {
                        point_location_map[pid[i]] = p[i];
                    }
                    face_point_map[static_cast<int>(fd.idx())].push_back(pid[i]);
                }

                vector<int> cur_bisector = {};
                cur_bisector.push_back(pid[0]);

                double len1 = sqrt(CGAL::squared_distance(p[1], p[0])),
                        len2 = sqrt(CGAL::squared_distance(p[2], p[0]));
                Base::Point p_end(p[1] + Base::Vector(p[2] - p[1]) * len1 / (len1 + len2));
                Base::Vector vec_bisector(p_end - p[0]);
                double limit_distance = sqrt(CGAL::squared_distance(p[0], p_end));
                double cur_distance = 0;
                int k = 0;

                while (Base::doubleCmp(cur_distance - limit_distance) < 0 && Base::doubleCmp(limit_distance / point_num) > 0) {
                    Base::Point bisector_p = p[0] + static_cast<double>(++k) / point_num * vec_bisector;
                    point_location_map[num_vertices] = bisector_p;
                    point_face_map[num_vertices] = static_cast<int>(fd.idx());
                    face_point_map[static_cast<int>(fd.idx())].push_back(num_vertices);
                    cur_bisector.push_back(num_vertices++);
                    cur_distance += limit_distance / point_num;
                }
                bisector_point_map[bisector_num] = cur_bisector;
                edge_bisector_map[eid].push_back(bisector_num);
                edge_bisector_map[eid2].push_back(bisector_num++);
            }
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return make_pair(static_cast<double>(duration.count()), num_vertices);
    }

    pair<double, int> placeSteinerPointsKAlgo(Base::Mesh &mesh, const int &K, const double &l_min,
                                               map<int, vector<int>> &edge_point_map,
                                               map<int, int> &point_face_map,
                                               map<int, Base::Point> &point_location_map,
                                               map<int, vector<int>> &face_point_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        int num_vertices = static_cast<int>(mesh.num_vertices());
        for (auto fd: mesh.faces()) {

            vector<Base::Point> p(3);
            vector<int> pid(3);
            int i = 0;
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))) {
                p[i] = mesh.points()[mesh.source(hed)];
                pid[i] = static_cast<int>(mesh.source(hed).idx());
                if (point_location_map.find(pid[i]) == point_location_map.end()) {
                    point_location_map[pid[i]] = p[i];
                }
                face_point_map[static_cast<int>(fd.idx())].push_back(pid[i]);
                i++;
            }

//            cout << "i = " << i << endl;

            double theta_m = 1e60;
            theta_m = min(theta_m, CGAL::approximate_angle(p[0], p[1], p[2]));
            theta_m = min(theta_m, CGAL::approximate_angle(p[1], p[2], p[0]));
            theta_m = min(theta_m, CGAL::approximate_angle(p[2], p[0], p[1]));
            double delta_I = l_min / 2 / sqrt(2 * (1 - cos(theta_m))) / K;

            for (auto hed:mesh.halfedges_around_face(mesh.halfedge(fd))) {
                auto p0 = mesh.points()[mesh.source(hed)],
                    p1 = mesh.points()[mesh.target(hed)];
                int eid = static_cast<int>(mesh.edge(hed).idx());
                if (edge_point_map.find(eid) == edge_point_map.end()){
                    edge_point_map[eid] = {};
                }
                edge_point_map[eid].push_back(mesh.source(hed).idx());

                Base::Vector vec01(p1 - p0);
                double limit_distance = sqrt(CGAL::squared_distance(p0, p1));
                int num_point_this_edge = floor(limit_distance / delta_I);
//                cout << "limit_distance = " << limit_distance << " # of points this edge = " << num_point_this_edge << endl;

                for (auto i = 1; i < num_point_this_edge; i++){
                    Base::Point next_p = p0 + 1.0 * i / num_point_this_edge * vec01;
                    point_location_map[num_vertices] = next_p;
                    point_face_map[num_vertices] = static_cast<int>(fd.idx());
                    face_point_map[static_cast<int>(fd.idx())].push_back(num_vertices);
                    edge_point_map[eid].push_back(num_vertices++);
                }
            }
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return make_pair(static_cast<double>(duration.count()), num_vertices);
    }

    double constructKAlgoGraph(Base::Mesh &mesh,
                              vector<double> &face_weight,
                              int num_vertices,
                              map<int, vector<int>> &edge_point_map,
                              map<int, int> &point_face_map,
                              map<int, Base::Point> &point_location_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        bool g_flag = true;
        if (kSkip::my_base_graph.num_V != num_vertices){
            g_flag = false;
            kSkip::my_base_graph.init(num_vertices);
        }
        vector<pair<int, int> > base_graph_edges = {};
        vector<double> base_graph_weights = {};

        for (auto fd: mesh.faces()){
            vector<int> eids;
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))){
                eids.push_back(mesh.edge(hed).idx());
            }
            for (auto i = 0; i < 3; i++){
                for (auto pid_1: edge_point_map[eids[i]]){
                    // connect the other two edges
                    for (auto pid_2: edge_point_map[eids[(i + 1) % 3]]){
                        double dis = sqrt(CGAL::squared_distance(point_location_map[pid_1], point_location_map[pid_2]));
                        base_graph_edges.emplace_back(pid_1, pid_2);
                        base_graph_weights.push_back(dis);
                        if (!g_flag) kSkip::my_base_graph.addEdge(pid_1, pid_2, dis);
                    }

                    for (auto pid_2: edge_point_map[eids[(i + 2) % 3]]){
                        double dis = sqrt(CGAL::squared_distance(point_location_map[pid_1], point_location_map[pid_2]));
                        base_graph_edges.emplace_back(pid_1, pid_2);
                        base_graph_weights.push_back(dis);
                        if (!g_flag) kSkip::my_base_graph.addEdge(pid_1, pid_2, dis);
                    }
                }
            }
        }

        base_graph = Base::Graph(begin(base_graph_edges), end(base_graph_edges), begin(base_graph_weights),
                                 num_vertices);


        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }


    double constructBaseGraph(Base::Mesh &mesh,
                              vector<double> &face_weight,
                              int num_vertices,
                              map<int, vector<int>> &edge_bisector_map,
                              map<int, vector<int>> &bisector_point_map,
                              map<int, int> &point_face_map,
                              map<int, Base::Point> &point_location_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        bool g_flag = true;
        if (kSkip::my_base_graph.num_V != num_vertices){
            g_flag = false;
            kSkip::my_base_graph.init(num_vertices);
        }
        vector<pair<int, int> > base_graph_edges = {};
        vector<double> base_graph_weights = {};
        for (auto it = edge_bisector_map.begin(); it != edge_bisector_map.end(); it++) {
            auto neighbor_bisectors = it->second;
            auto common_eid = it->first;
            for (auto &bisector_1: neighbor_bisectors) {
                for (auto &bisector_2: neighbor_bisectors) {
                    for (auto i = 0; i != bisector_point_map[bisector_1].size(); i++) {
                        for (auto j = 0; j != bisector_point_map[bisector_2].size(); j++) {
                            auto pid_1 = bisector_point_map[bisector_1][i];
                            auto pid_2 = bisector_point_map[bisector_2][j];
                            if (pid_1 == pid_2) continue;
                            auto fid_1 = i ? pid_1 : bisector_point_map[bisector_1][i + 1],
                                    fid_2 = j ? pid_2 : bisector_point_map[bisector_2][j + 1];
                            fid_1 = point_face_map[fid_1];
                            fid_2 = point_face_map[fid_2];
//                            auto point_1 = point_location_map[pid_1];
//                            auto point_2 = point_location_map[pid_2];
                            double dis = Base::distanceSnell(mesh, face_weight,
                                                             point_location_map[pid_1], fid_1,
                                                             point_location_map[pid_2], fid_2,
                                                             common_eid);
//                            if (pid_1 == 7) {
//                                cout << "dbg: " << pid_1 << " " << pid_2 << " " << dis << endl;
//                            }
                            base_graph_edges.emplace_back(pid_1, pid_2);
                            base_graph_weights.push_back(dis);
                            if (!g_flag) kSkip::my_base_graph.addEdge(pid_1, pid_2, dis);
                        }
                    }
                }
            }
        }
        base_graph = Base::Graph(begin(base_graph_edges), end(base_graph_edges), begin(base_graph_weights),
                                 num_vertices);


        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    // method type{0: enlarged surface disk; 1: surface disk} # Do not use method 1!!!
    pair<double, double> distanceBoundMapPreprocessing(Base::Mesh &mesh, int method_type,
                                                       map<int, int> &point_face_map,
                                                       vector<double> &face_max_length,
                                                       map<int, Base::Point> &point_location_map,
                                                       vector<int> &pid_list) {
        double bound_pre_time = 0.0;
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        distance_map.clear();
        bound_map.clear();

        auto geodesic_pre_start = chrono::_V2::system_clock::now();
        Base::Surface_mesh_shortest_path shortest_paths(mesh);
        auto geodesic_pre_end = chrono::_V2::system_clock::now();
        auto geodesic_pre_duration = chrono::duration_cast<chrono::milliseconds>(geodesic_pre_end - geodesic_pre_start);
        auto geodesic_pre_time = static_cast<double>(geodesic_pre_duration.count());

        auto V = pid_list.size();
        int cnt = 1;
        vector<double> dall(V, -1);
        vector<int> fa(V, -1);
        vector<double> b(V, 0);
        for (auto i = 0; i < V; i++) {
            auto s = pid_list[i];
            if (i > cnt * V / 10) {
                cout << cnt++ << "0% nodes are preprocessed." << endl;
            }
            dall.resize(num_vertices(base_graph), Base::unreachable);
            fa.resize(num_vertices(base_graph), -1);
            b.resize(num_vertices(base_graph), 0);
            if (!method_type) {
                Base::dijkstra_SSAD(base_graph, s, dall, fa);
            }
            else {
                auto vd = *(mesh.vertices().begin() + s);
                shortest_paths.add_source_point(vd);
            }
            vector<double> d(dall.begin(), dall.begin() + num_vertices(base_graph));
            for (auto j = 0; j < V; j++) {
                auto t = pid_list[j];
                if (!method_type) {
                    auto bound_pre_start = chrono::_V2::system_clock::now();

//                    int id = t;
//                    while (id != s) {
//                        if (id < V && fa[id] < V) {
//                            b[t] += sqrt(CGAL::squared_distance(point_location_map[id], point_location_map[fa[id]]));
//                        } else {
//                            b[t] += face_max_length[point_face_map[id]];
//                        }
//                        id = fa[id];
//                    }

                    auto bound_pre_end = chrono::_V2::system_clock::now();
                    auto bound_pre_duration = chrono::duration_cast<chrono::milliseconds>(bound_pre_end - bound_pre_start);
                    bound_pre_time += static_cast<double>(bound_pre_duration.count());
                }
                else {
                    auto target_vd = *(mesh.vertices().begin() + t);
                    auto dis_pair = shortest_paths.shortest_distance_to_source_points(target_vd);
                    d[t] = dis_pair.first;
                }
            }
            if (method_type == 1) {
                shortest_paths.remove_all_source_points();
            }
            distance_map[s] = d;
            bound_map[s] = d;
//            bound_map[s] = b;

//            cout << "s = " << s;
//            for (auto x = 0; x < V; x++){
//                if (!Base::doubleCmp(d[pid_list[x]])){
//                    cout << " id = " << pid_list[x] << " dis = " << d[pid_list[x]] << endl;
//                    cout << point_location_map[s] << " on face " << point_face_map[s] << endl;
//                    cout << point_location_map[pid_list[x]] << " on face " << point_face_map[pid_list[x]] << endl;
//                }
//            }
        }
        cout << "Distance/Bound map preprocessing finished." << endl;

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        if (!method_type) {
//            cout << "[Debug Info] Bound pre time = " << bound_pre_time << " ms" << endl;
            return make_pair(static_cast<double>(duration.count()) - geodesic_pre_time - bound_pre_time, bound_pre_time);
        } else {
            return make_pair(static_cast<double>(duration.count()), 0.0);
        }
    }

    class PartitionTreeNode {
    public:
        int node_id{};
        int center_idx;
        double radius;
        int level;
        PartitionTreeNode *parent;
        map<int, PartitionTreeNode *> children;
        map<int, bool> cover_set;

        bool operator()(PartitionTreeNode *a, PartitionTreeNode *b) {
            return a->center_idx < b->center_idx;
        }

        PartitionTreeNode();

        PartitionTreeNode(int center_idx_value, double radius_value, int node_id_value);

        ~PartitionTreeNode();

    private:
    };

    inline PartitionTreeNode::PartitionTreeNode() {
        parent = nullptr;
        children.clear();
        center_idx = -1;
        level = -1;
        radius = -1.0;
    }

    inline PartitionTreeNode::PartitionTreeNode(int center_idx_value, double radius_value, int node_id_value) {
        parent = nullptr;
        children.clear();
        center_idx = center_idx_value;
        radius = radius_value;
        node_id = node_id_value;
    }

    inline PartitionTreeNode::~PartitionTreeNode() {
    }

    typedef pair<double, pair<PartitionTreeNode *, PartitionTreeNode *> > nodePairWithValue;
    typedef pair<PartitionTreeNode*, PartitionTreeNode*> nodePair;
    map<pair<int, int>, double> enhanced_edges; // store the enhanced edges

    class PartitionTree {
    public:
        PartitionTreeNode *root;
        int max_level;
        int num_tree_nodes;
        vector<vector<PartitionTreeNode *> > level_nodes;
        bool stop_flag;

        PartitionTree();
        ~PartitionTree();

        double constructTree(int num_vertices);

        double generateNodePairSet(double eps, int type, int bisector_point_num, set<nodePair> &node_pair_set);

        static void getPathToRoot(PartitionTreeNode* node, vector<PartitionTreeNode*> &path);

        void rootNodeSelect(vector<int> &pid_list);

        int buildLevel(vector<int> &pid_list, double l);

        double constructTree(ofstream &fout, vector<int> &pid_list, double l);

        double findEnhancedDistance(int uid, int vid);
    };

    PartitionTree::PartitionTree() { stop_flag = 0; num_tree_nodes = 0; }

    PartitionTree::~PartitionTree() {}

    inline void PartitionTree::getPathToRoot(PartitionTreeNode *node, vector<PartitionTreeNode *> &path) {
        path.clear();
        while (node != nullptr){
            path.push_back(node);
            if (node->parent != nullptr){
                node = node->parent;
            }
            else{
                node = nullptr;
            }
        }
    }

    inline double wellSeparatedValue(double center_distance, double eps, int type,
                                     PartitionTreeNode* node1, PartitionTreeNode* node2){
        if (!Base::doubleCmp(node1->radius) && !Base::doubleCmp(node2->radius)) return 0;
        assert(Base::doubleCmp(eps));
        double rhs = (2 / eps + 2) * max(node1->radius, node2->radius);
        if (!type){
            rhs = (8 / eps + 4) * max(node1->radius, node2->radius);
        }
        return rhs;
    }

    double PartitionTree::findEnhancedDistance(int uid, int vid){
        vector<PartitionTreeNode*> u_path;
        vector<PartitionTreeNode*> v_path;
        auto u_node = level_nodes[max_level][uid];
        auto v_node = level_nodes[max_level][vid];
        assert(u_node->center_idx == uid);
        assert(v_node->center_idx == vid);
        getPathToRoot(u_node, u_path);
        getPathToRoot(v_node, v_path);
        for (auto i = 0; i < u_path.size(); i++){
            if (enhanced_edges.find(make_pair(u_path[i]->center_idx, v_path[i]->center_idx)) != enhanced_edges.end()){
                return enhanced_edges[make_pair(u_path[i]->center_idx, v_path[i]->center_idx)];
            }
        }
        cout << "<ERROR> No Enhanced Edges found." << endl;
        return -1;
    }

    double PartitionTree::generateNodePairSet(double eps, int type, int bisector_point_num, set<nodePair> &node_pair_set) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        priority_queue<nodePairWithValue , vector<nodePairWithValue>, greater<nodePairWithValue> > q;
        double center_distance = 0.0;
        double ws_value = -1.0;

        ws_value = wellSeparatedValue(center_distance, eps, type, root, root);

        q.push(make_pair(0 - ws_value, make_pair(root, root)));
        while (Base::doubleCmp(q.top().first) < 0){
            nodePairWithValue cur_node_pair;
            cur_node_pair = q.top(); q.pop();
            auto l_node = cur_node_pair.second.first;
            auto r_node = cur_node_pair.second.second;
            //  Split the left node
            if (Base::doubleCmp(l_node->radius - r_node->radius) > 0 ||
                !Base::doubleCmp(l_node->radius - r_node->radius)
                && l_node->center_idx <= r_node->center_idx){
                for (auto &child: l_node->children){
                    int child_idx = child.first;
                    auto child_node = child.second;
//                    center_distance = distance_map[child_idx][r_node->center_idx];
                    center_distance = enhanced_edges[make_pair(child_idx, r_node->center_idx)];

                    ws_value = wellSeparatedValue(center_distance, eps, type, child_node, r_node);
                    q.push(make_pair(center_distance - ws_value, make_pair(child_node, r_node)));
                }
            }
            //  Split the right node
            else{
                for (auto &child: r_node->children){
                    int child_idx = child.first;
                    auto child_node = child.second;
//                    center_distance = distance_map[l_node->center_idx][child_idx];
                    center_distance = enhanced_edges[make_pair(l_node->center_idx, child_idx)];

                    ws_value = wellSeparatedValue(center_distance, eps, type, l_node, child_node);

                    q.push(make_pair(center_distance - ws_value, make_pair(l_node, child_node)));
                }
            }
        }
        node_pair_set.clear();
        while (!q.empty()){
            node_pair_set.insert(q.top().second);
            q.pop();
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    double PartitionTree::constructTree(ofstream &fout, vector<int> &pid_list, double l) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer
        auto num_vertices = pid_list.size();
        rootNodeSelect(pid_list);
        while (!stop_flag){
            int cur_num_level_nodes = buildLevel(pid_list, l);
            int max_cover_set_size = -1, min_cover_set_size = num_vertices + 1;
            for (auto &node: level_nodes[max_level]){
                max_cover_set_size = max(max_cover_set_size, static_cast<int>(node->cover_set.size()));
                min_cover_set_size = min(min_cover_set_size, static_cast<int>(node->cover_set.size()));
            }
            fout << "level " << max_level << " has " << cur_num_level_nodes << " nodes" << endl;
            fout << "max cover set size: " << max_cover_set_size << " min cover set size: " << min_cover_set_size << endl;
        }
        fout << "[Partition Tree] totally " << num_tree_nodes << " nodes " << endl;

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    void PartitionTree::rootNodeSelect(vector<int> &pid_list) {
        auto num_vertices = pid_list.size();
        int root_index = pid_list[rand() % num_vertices];   // random select one point
//        vector<double> d(distance_map[root_index]);
        vector<double> d; kSkip::bounded_dijkstra(kSkip::my_base_graph, root_index, Base::unreachable, d);
        double radius = -1.0;
        for (auto i = 0; i < num_vertices; i++){
            if (pid_list[i] == root_index) continue;
            double distance = d[pid_list[i]];
            if (Base::doubleCmp(radius - distance) < 0){
                radius = distance;
            }
        }
        root = new PartitionTreeNode(root_index, radius, num_tree_nodes++);
        max_level = 0;
        root->level = max_level;
        vector<PartitionTreeNode*> cur_level_nodes(1, root);
        level_nodes.push_back(cur_level_nodes);
        for (auto i = 0; i < num_vertices; i++){
            root->cover_set[pid_list[i]] = true;
        }
//        cout << "[Partition Tree] Root = (" << root_index << "," << radius << ")" << endl;
    }

    int PartitionTree::buildLevel(vector<int> &pid_list, double l){
        auto num_vertices = pid_list.size();
        if (stop_flag){
            cout << "[Partition Tree] Already finish the building of partition tree." << endl;
            return -1;
        }
//        max_level++;
        vector<PartitionTreeNode*> cur_level_nodes = {};
        map<int, bool> father_candidates;
        set<int> vertex_to_be_covered = {};
        for (auto i = 0; i < num_vertices; i++){
            vertex_to_be_covered.insert(pid_list[i]);
        }
        double cur_level_radius = level_nodes[max_level][0]->radius * 0.5;
//        cout << "[Partition Tree] Current level radius is: " << cur_level_radius << endl;
        for (auto &last_level_node: level_nodes[max_level]){
            int vertex_idx = last_level_node->center_idx;
            PartitionTreeNode* cur_node = new PartitionTreeNode(vertex_idx, cur_level_radius, num_tree_nodes++);
//            vector<double> d = distance_map[vertex_idx];
            vector<double> d; kSkip::bounded_dijkstra(kSkip::my_base_graph, vertex_idx, l * cur_level_radius, d);
            for (auto it = vertex_to_be_covered.begin(); it != vertex_to_be_covered.end(); ){
                int id = *it;
                assert(Base::doubleCmp(d[id]) >= 0);
                if (Base::doubleCmp(d[id] - cur_level_radius) <= 0){
//                    cout << "s = " << vertex_idx << " id = " << id << " dis = " << d[id] << endl;
                    cur_node->cover_set[id] = true;
                    it = vertex_to_be_covered.erase(it);
                }
                else{
                    it++;
                }
            }
            //            for (auto i = 0; i < d.size(); i++){
            for (auto i = 0; i < pid_list.size(); i++){
                    int id = pid_list[i];
                    if (Base::doubleCmp(d[id] - l * cur_level_radius) < 0){
                    if (!(enhanced_edges.find(make_pair(vertex_idx, id)) == enhanced_edges.end()))continue;
                    enhanced_edges[make_pair(vertex_idx, id)] = d[id];
                }
            }
            assert(!cur_node->cover_set.empty());
            cur_node->parent = last_level_node;
            last_level_node->children[vertex_idx] = cur_node;
            cur_node->level = max_level + 1;
            cur_level_nodes.push_back(cur_node);
        }
        //  Still has nodes to be covered
        if (!vertex_to_be_covered.empty()){
            vector<int> vertex_sequence(vertex_to_be_covered.begin(), vertex_to_be_covered.end());
#ifdef PrintDetails
            cout << "to be covered vertices: "; for (auto val: vertex_sequence){ cout << val << " "; } cout << endl;
#endif
            unsigned seed = chrono::system_clock::now().time_since_epoch().count();
            shuffle(vertex_sequence.begin(), vertex_sequence.end(), default_random_engine(seed));
            for (auto &vertex_idx: vertex_sequence){
                //  Vertex already covered by some other vertex.
                if (vertex_to_be_covered.find(vertex_idx) == vertex_to_be_covered.end()) continue;
                PartitionTreeNode* cur_node = new PartitionTreeNode(vertex_idx, cur_level_radius, num_tree_nodes++);
//                vector<double> d = distance_map[vertex_idx];
                vector<double> d; kSkip::bounded_dijkstra(kSkip::my_base_graph, vertex_idx, l * cur_level_radius, d);
                for (auto it = vertex_to_be_covered.begin(); it != vertex_to_be_covered.end(); ){
                    int id = *it;
                    if (Base::doubleCmp(d[id] - cur_level_radius) <= 0){
                        cur_node->cover_set[id] = true;
                        it = vertex_to_be_covered.erase(it);
                    }
                    else{
                        it++;
                    }
                }
//                for (auto i = 0; i < d.size(); i++){
                for (auto i = 0; i < pid_list.size(); i++){
                    int id = pid_list[i];
                    if (Base::doubleCmp(d[id] - l * cur_level_radius) < 0){
                        if (!(enhanced_edges.find(make_pair(vertex_idx, id)) == enhanced_edges.end()))continue;
                        enhanced_edges[make_pair(vertex_idx, id)] = d[id];
                    }
                }
                assert(!cur_node->cover_set.empty());  //  Each node cover itself at least
                int father_idx = -1;
                PartitionTreeNode* father;
                for (auto &last_level_node: level_nodes[max_level]){
                    int fa_id = last_level_node->center_idx;
                    assert(Base::doubleCmp(d[father_idx]) >= 0);
                    assert(Base::doubleCmp(d[fa_id]) >= 0);
                    if (father_idx < 0 || Base::doubleCmp(d[fa_id] - d[father_idx]) < 0){
                        father_idx = fa_id;
                        father = last_level_node;
                    }
                }
                assert(Base::doubleCmp(d[father_idx] - cur_level_radius * 2.0) <= 0);
                cur_node->parent = father;
                father->children[vertex_idx] = cur_node;
                cur_node->level = max_level + 1;
                cur_level_nodes.push_back(cur_node);
            }
        }
        assert(vertex_to_be_covered.empty());   //  All vertices will be covered.
        max_level++;
        if (cur_level_nodes.size() == num_vertices){
            stop_flag = true;
            sort(cur_level_nodes.begin(), cur_level_nodes.end(), PartitionTreeNode());
            for (auto &leaf_node: cur_level_nodes){
                leaf_node->radius = 0.0;
            }
        }
        level_nodes.push_back(cur_level_nodes);
        return static_cast<int>(cur_level_nodes.size());
    }

    double distanceQueryBf(set<nodePair> &node_pairs,
                           vector<PartitionTreeNode*> &As,
                           vector<PartitionTreeNode*> &At){
        double approximate_distance = 0.0;
        int cnt = 0;
        for (auto i = 0; i != As.size(); i++){
            for (auto j = 0; j != At.size(); j++){
                nodePair query_pair = make_pair(As[i], At[j]);
                if (node_pairs.find(query_pair) != node_pairs.end()){
                    cnt++;
//                    approximate_distance = distance_map[As[i]->center_idx][At[j]->center_idx];
                    approximate_distance = enhanced_edges[make_pair(As[i]->center_idx, At[j]->center_idx)];
                }
            }
        }
        assert(cnt == 1);
        if (cnt != 1){
            cout << "[ERROR] multiple node pairs found!!! cnt = " << cnt << endl;
        }
        return approximate_distance;
    }

    double distanceQueryEfficient(set<nodePair> &node_pairs,
                                  vector<PartitionTreeNode*> &As,
                                  vector<PartitionTreeNode*> &At){
        double approximate_distance = 0.0;
        int cnt = 0;
        // same-layer
        for (auto i = 0; i < As.size(); i++){
            nodePair query_pair = make_pair(As[i], At[i]);
            if (node_pairs.find(query_pair) != node_pairs.end()){
                cnt++;
                approximate_distance = enhanced_edges[make_pair(As[i]->center_idx, At[i]->center_idx)];
            }
        }

        //first-higher
        for (auto i = 1; i < As.size(); i++){
            nodePair query_pair = make_pair(As[i - 1], At[i]);
            if (node_pairs.find(query_pair) != node_pairs.end()){
                cnt++;
                approximate_distance = enhanced_edges[make_pair(As[i - 1]->center_idx, At[i]->center_idx)];
            }
        }

        //first-lower
        for (auto i = 1; i < As.size(); i++){
            nodePair query_pair = make_pair(As[i], At[i - 1]);
            if (node_pairs.find(query_pair) != node_pairs.end()){
                cnt++;
                approximate_distance = enhanced_edges[make_pair(As[i]->center_idx, At[i - 1]->center_idx)];
            }
        }
        assert(cnt == 1);
        if (cnt != 1){
            cout << "Multiple distance found! cnt = " << cnt << endl;
        }
        return approximate_distance;
    }

    double distanceQueryBf(double source_dis,
                           double target_dis,
                           PartitionTree &tree,
                           set<nodePair> &node_pairs,
                           vector<PartitionTreeNode*> &As,
                           vector<PartitionTreeNode*> &At,
                           map<int, int> &new_id){
        double approximate_distance = -1.0;
        int nearest_sid = As[0]->center_idx, nearest_tid = At[0]->center_idx, cnt = 0;

        for (auto i = 0; i != As.size(); i++){
            for (auto j = 0; j != At.size(); j++){
                nodePair query_pair = make_pair(As[i], At[j]);
//                cout << "As.radius = " << As[i]->radius << " source_dis = " << source_dis << endl;
//                cout << "At.radius = " << At[i]->radius << " target_dis = " << target_dis << endl;

                vector<PartitionTreeNode*> tAs, tAt;

                tree.getPathToRoot(tree.level_nodes[tree.max_level][new_id[As[i]->center_idx]], tAs);
                tree.getPathToRoot(tree.level_nodes[tree.max_level][new_id[nearest_sid]], tAt);

                bool is_cover_source = false, is_cover_target = false;
//                if (Base::doubleCmp(As[i]->radius - (distance_map[As[i]->center_idx][nearest_sid] + source_dis)) >= 0){
                double t_dis = distanceQueryBf(node_pairs, tAs, tAt);
//                if (Base::doubleCmp(As[i]->radius - (enhanced_edges[make_pair(As[i]->center_idx, nearest_sid)] + source_dis)) >= 0){
                if (Base::doubleCmp(As[i]->radius - (t_dis + source_dis)) >= 0){
                        is_cover_source = true;
                }
//                if (Base::doubleCmp(At[j]->radius - (distance_map[At[j]->center_idx][nearest_tid] + target_dis)) >= 0){
//                if (Base::doubleCmp(At[j]->radius - (enhanced_edges[make_pair(At[j]->center_idx, nearest_tid)] + target_dis)) >= 0){
                tree.getPathToRoot(tree.level_nodes[tree.max_level][new_id[At[j]->center_idx]], tAs);
                tree.getPathToRoot(tree.level_nodes[tree.max_level][new_id[nearest_tid]], tAt);
                t_dis = distanceQueryBf(node_pairs, tAs, tAt);

                if (Base::doubleCmp(At[j]->radius - (t_dis + target_dis)) >= 0){
                    is_cover_target = true;
                }

                if (is_cover_source && is_cover_target && node_pairs.find(query_pair) != node_pairs.end()){
                    cnt++;
//                    approximate_distance = distance_map[As[i]->center_idx][At[j]->center_idx];
                    approximate_distance = enhanced_edges[make_pair(As[i]->center_idx, At[j]->center_idx)];

                }
            }
        }
        if (cnt > 0){
            assert(cnt == 1);
            if (cnt != 1){
                cout << "[ERROR] Multiple node pairs found!" << endl;
            }
        }
        return approximate_distance;
    }

    pair<int, double> findNearestVertex(Base::Point &p,
                          Base::Mesh &mesh,
                          Base::AABB_tree &aabb_tree){
        auto location = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(p, aabb_tree, mesh);
        auto fd = location.first;
        double nearest_dis = 1e100;
        int nearest_id = -1;
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd))){
            double dis = CGAL::squared_distance(mesh.points()[vd], p);
//            cout << "dis = " << dis << endl;
            if (Base::doubleCmp(dis - nearest_dis) < 0){
                nearest_dis = dis;
                nearest_id = vd.idx();
            }
        }
        return make_pair(nearest_id, sqrt(nearest_dis));
    }

    double distanceQueryA2A(Base::Point &query_source,
                          int fid_s,
                          Base::Point &query_target,
                          int fid_t,
                          PartitionTree &tree,
                          map<int, vector<int> > &face_point_map,
                          map<int, Base::Point> &point_location_map,
                          set<nodePair> &node_pairs
                          ){
        double ret_dis = Base::unreachable;
        auto leaf_nodes = tree.level_nodes[tree.max_level];
        for (auto pid1: face_point_map[fid_s]){
            for (auto pid2: face_point_map[fid_t]){
                vector<PartitionTreeNode*> As, At;
                tree.getPathToRoot(leaf_nodes[pid1], As);
                tree.getPathToRoot(leaf_nodes[pid2], At);
                double t_dis = sqrt(CGAL::squared_distance(query_source, point_location_map[pid1])) +
                        distanceQueryEfficient(node_pairs, As, At) +
                        sqrt(CGAL::squared_distance(point_location_map[pid2], query_target));
                if (Base::doubleCmp(t_dis - ret_dis) < 0) ret_dis = t_dis;
            }
        }
        if (!Base::doubleCmp(ret_dis - Base::unreachable)){
            cout << "ERROR: A2A distance not found!" << endl;
        }
        return ret_dis;
    }

}

#endif //WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H
