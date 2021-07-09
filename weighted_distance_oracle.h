//
// Created by huang on 2021/7/6.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H
#define WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H

#include "base.h"

namespace WeightedDistanceOracle {

    using namespace std;

    Base::Graph base_graph;
    map<int, vector<double> > distance_map, bound_map;

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

    pair<double, int> placeBisectorPointsFixed(Base::Mesh &mesh, const int &point_num,
                                               map<int, vector<int>> &edge_bisector_map,
                                               map<int, vector<int>> &bisector_point_map,
                                               map<int, int> &point_face_map,
                                               map<int, Base::Point> &point_location_map) {
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
#ifdef PrintDetails
                cout << "current V = " << pid[0] << endl;
#endif
                while (Base::doubleCmp(cur_distance - limit_distance) < 0) {
                    Base::Point bisector_p = p[0] + static_cast<double>(++k) / point_num * vec_bisector;
#ifdef PrintDetails
                    cout << "k = " << k << ", p = " << bisector_p << endl;
#endif
                    point_location_map[num_vertices] = bisector_p;
                    point_face_map[num_vertices] = static_cast<int>(fd.idx());
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

    double constructBaseGraph(Base::Mesh &mesh,
                              vector<double> &face_weight,
                              int num_vertices,
                              map<int, vector<int>> &edge_bisector_map,
                              map<int, vector<int>> &bisector_point_map,
                              map<int, int> &point_face_map,
                              map<int, Base::Point> &point_location_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

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
                            auto point_1 = point_location_map[pid_1];
                            auto point_2 = point_location_map[pid_2];
                            double dis = Base::distanceSnell(mesh, face_weight,
                                                             point_location_map[pid_1], fid_1,
                                                             point_location_map[pid_2], fid_2,
                                                             common_eid);
                            base_graph_edges.emplace_back(pid_1, pid_2);
                            base_graph_weights.push_back(dis);
#ifdef PrintDetails
                            cout << "add edge: " << pid_1 << " " << pid_2 << " " << dis << endl;
#endif
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

    // method type{0: enlarged surface disk; 1: surface disk}
    pair<double, double> distanceBoundMapPreprocessing(Base::Mesh &mesh, int method_type,
                                         map<int, int> &point_face_map,
                                         vector<double> face_max_length) {
        double bound_pre_time;
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        distance_map.clear();
        bound_map.clear();

        auto geodesic_pre_start = chrono::_V2::system_clock::now();
        Base::Surface_mesh_shortest_path shortest_paths(mesh);
        auto geodesic_pre_end = chrono::_V2::system_clock::now();
        auto geodesic_pre_duration = chrono::duration_cast<chrono::milliseconds>(geodesic_pre_end - geodesic_pre_start);
        auto geodesic_pre_time = static_cast<double>(geodesic_pre_duration.count());

        int V = mesh.num_vertices();
        int cnt = 1;
        for (auto i = 0; i < V; i++) {
            if (i > cnt * V / 10) {
                cout << cnt++ << "0% nodes are preprocessed." << endl;
            }
            vector<double> dall(V, -1);
            vector<int> fa(V, -1);
            vector<double> b(V, 0);
            if (!method_type) {
                Base::dijkstra_SSAD(base_graph, i, dall, fa);
            }
            else {
                auto vd = *(mesh.vertices().begin() + i);
                shortest_paths.add_source_point(vd);
            }
            vector<double> d(dall.begin(), dall.begin() + V);
            for (auto j = 0; j < V; j++) {
                if (!method_type) {
                    auto bound_pre_start = chrono::_V2::system_clock::now();

                    int id = j;
                    while (id != i) {
                        if (id < V && fa[id] < V) {
                            auto vds = *(mesh.vertices().begin() + id),
                                    vdt = *(mesh.vertices().begin() + fa[id]);
                            b[j] += sqrt(CGAL::squared_distance(mesh.points()[vds], mesh.points()[vdt]));
                        } else {
                            b[j] += face_max_length[point_face_map[id]];
                        }
                        id = fa[id];
                    }

                    auto bound_pre_end = chrono::_V2::system_clock::now();
                    auto bound_pre_duration = chrono::duration_cast<chrono::milliseconds>(bound_pre_end - bound_pre_start);
                    bound_pre_time = static_cast<double>(bound_pre_duration.count());
                }
                else {
                    auto target_vd = *(mesh.vertices().begin() + j);
                    auto dis_pair = shortest_paths.shortest_distance_to_source_points(target_vd);
                    d[j] = dis_pair.first;
                }
            }
            if (method_type == 1) {
                shortest_paths.remove_all_source_points();
            }
            distance_map[i] = d;
            bound_map[i] = b;
        }
        cout << "Distance/Bound map preprocessing finished." << endl;

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        if (!method_type) {
            return make_pair(static_cast<double>(duration.count()) - geodesic_pre_time, bound_pre_time);
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

    class PartitionTree {
    public:
        PartitionTreeNode *root;
        int max_level;
        int num_tree_nodes;
        vector<vector<PartitionTreeNode *> > level_nodes;
        bool stop_flag;

        PartitionTree();
        ~PartitionTree();

        void root_node_select(int num_vertices);

        int build_level(int num_vertices);

        double construct_tree(int num_vertices);

        double generate_node_pair_set(double eps, int type, int bisector_point_num, set<nodePair> &node_pair_set);

        void get_path_to_root(PartitionTreeNode* node, vector<PartitionTreeNode*> &path);
    };

    PartitionTree::PartitionTree() { stop_flag = 0; num_tree_nodes = 0; }

    PartitionTree::~PartitionTree() {}

    inline void PartitionTree::get_path_to_root(PartitionTreeNode *node, vector<PartitionTreeNode *> &path) {
        path.clear();
        while (node != nullptr){
            path.push_back(node);
            node = node->parent;
        }
    }

    inline double well_separated_value(double center_distance, double eps, int type,
                                       PartitionTreeNode* node1, PartitionTreeNode* node2,
                                       int bisector_point_num = 6, double bound_dis = 1e100){
        if (!Base::doubleCmp(node1->radius) && !Base::doubleCmp(node2->radius)) return 0;
        assert(Base::doubleCmp(eps));
        double rhs = (2 / eps + 2) * max(node1->radius, node2->radius);
        if (!type){
            double C = sqrt(3) / 2 / bisector_point_num * bound_dis;
            rhs += C;
        }
        return rhs;
    }

    double PartitionTree::generate_node_pair_set(double eps, int type, int bisector_point_num, set<nodePair> &node_pair_set) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        priority_queue<nodePairWithValue , vector<nodePairWithValue>, greater<nodePairWithValue> > q;
        double center_distance = 0.0;
        double ws_value = -1.0;
        if (!type){
            ws_value = well_separated_value(center_distance, eps, type,
                                            root, root,
                                            bisector_point_num, 0);
        }
        else{
            ws_value = well_separated_value(center_distance, eps, type, root, root);
        }

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
                    center_distance = distance_map[child_idx][r_node->center_idx];
                    if (!type){
                        ws_value = well_separated_value(center_distance, eps, type,
                                                        child_node, r_node,
                                                        bisector_point_num, bound_map[child_idx][r_node->center_idx]);
                    }
                    else{
                        ws_value = well_separated_value(center_distance, eps, type,
                                                        child_node, r_node);
                    }
                    q.push(make_pair(center_distance - ws_value, make_pair(child_node, r_node)));
                }
            }
            //  Split the right node
            else{
                for (auto &child: r_node->children){
                    int child_idx = child.first;
                    auto child_node = child.second;
                    center_distance = distance_map[l_node->center_idx][child_idx];
                    if (!type){
                        ws_value = well_separated_value(center_distance, eps, type,
                                                        l_node, child_node,
                                                        bisector_point_num, bound_map[l_node->center_idx][child_idx]);
                    }
                    else{
                        ws_value = well_separated_value(center_distance, eps, type,
                                                        l_node, child_node);
                    }
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

    double PartitionTree::construct_tree(int num_vertices) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        root_node_select(num_vertices);
        while (!stop_flag){
            int cur_num_level_nodes = build_level(num_vertices);
            int max_cover_set_size = -1, min_cover_set_size = num_vertices + 1;
            for (auto &node: level_nodes[max_level]){
                max_cover_set_size = max(max_cover_set_size, static_cast<int>(node->cover_set.size()));
                min_cover_set_size = min(min_cover_set_size, static_cast<int>(node->cover_set.size()));
            }
            cout << "level " << max_level << " has " << cur_num_level_nodes << " nodes" << endl;
            cout << "max cover set size: " << max_cover_set_size << " min cover set size: " << min_cover_set_size << endl;
        }
        cout << "[Partition Tree] totally " << num_tree_nodes << " nodes " << endl;

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    void PartitionTree::root_node_select(int num_vertices) {
        srand((int)time(0));
        int root_index = rand() % num_vertices;
        vector<double> d(distance_map[root_index]);
        double radius = -1.0;
        for (auto i = 0; i < num_vertices; i++){
            if (i == root_index) continue;
            double distance = d[i];
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
            root->cover_set[i] = true;
        }
        cout << "[Partition Tree] Root = (" << root_index << "," << radius << ")" << endl;
    }

    int PartitionTree::build_level(int num_vertices){
        if (stop_flag){
            cout << "[Partition Tree] Already finish the building of partition tree." << endl;
            return -1;
        }
//        max_level++;
        vector<PartitionTreeNode*> cur_level_nodes = {};
        map<int, bool> father_candidates;
        vector<int> mesh_vertices(num_vertices);
        set<int> vertex_to_be_covered = {};
        for (auto i = 0; i < num_vertices; i++){
            vertex_to_be_covered.insert(i);
        }
        double cur_level_radius = level_nodes[max_level][0]->radius * 0.5;
        cout << "[Partition Tree] Current level radius is: " << cur_level_radius << endl;
        for (auto &last_level_node: level_nodes[max_level]){
            int vertex_idx = last_level_node->center_idx;
            PartitionTreeNode* cur_node = new PartitionTreeNode(vertex_idx, cur_level_radius, num_tree_nodes++);
            vector<double> d = distance_map[vertex_idx];
            for (auto it = vertex_to_be_covered.begin(); it != vertex_to_be_covered.end(); ){
                int id = *it;
                assert(Base::doubleCmp(d[id]) >= 0);
                if (Base::doubleCmp(d[id] - cur_level_radius) <= 0){
                    cur_node->cover_set[id] = true;
                    it = vertex_to_be_covered.erase(it);
                }
                else{
                    it++;
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
                vector<double> d = distance_map[vertex_idx];
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

    double distance_query_bf(set<nodePair> &node_pairs,
                             vector<PartitionTreeNode*> &As,
                             vector<PartitionTreeNode*> &At){
        double approximate_distance = 0.0;
        int cnt = 0;
        for (auto i = 0; i != As.size(); i++){
            for (auto j = 0; j != At.size(); j++){
                auto u = As[i], v = At[j];
                nodePair query_pair = make_pair(As[i], At[j]);
                if (node_pairs.find(query_pair) != node_pairs.end()){
                    cnt++;
                    approximate_distance = distance_map[As[i]->center_idx][At[j]->center_idx];
                }
            }
        }
        assert(cnt == 1);
        if (cnt != 1){
            cout << "[ERROR] multiple node pairs found!!!" << endl;
        }
        return approximate_distance;
    }

}

#endif //WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H
