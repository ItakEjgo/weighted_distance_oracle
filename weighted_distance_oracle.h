#ifndef WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H
#define WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H

#include "base.h"
#include "k_skip.h"

namespace WeightedDistanceOracle {

    using namespace std;

    Base::Graph base_graph;
    map<int, vector<float> > distance_map, bound_map;

    //get the value of Gamme for each vertex of the terrain. Ref. JACM'2005 bisector
    vector<float> getVertexGamma(const Base::Mesh &mesh, const vector<float> &face_weights){
        map<unsigned, vector<unsigned> > vertex_face_map = {};
        vector<float> gama(mesh.num_vertices());
        float dis, w_min, w_max;
        for (auto fd: mesh.faces()) {
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))) {
                unsigned source_idx = mesh.source(hed);
                assert(mesh.face(hed) != mesh.null_face());
                vertex_face_map[source_idx].push_back(mesh.face(hed).idx());
            }
        }

        for (auto vd: mesh.vertices()){
            unsigned vid = vd.idx();
            Base::Point v(mesh.points()[vd]);
            dis = 1e100;
            w_max = -1.0;
            w_min = 1e100;
            for (auto fid: vertex_face_map[vid]){
                w_max = Base::floatCmp(face_weights[fid] - w_max) > 0 ? face_weights[fid] : w_max;
                w_min = Base::floatCmp(face_weights[fid] - w_min) < 0 ? face_weights[fid] : w_min;
                float tmp_face_dis = -1.0;
                for (auto &hed: mesh.halfedges_around_face(mesh.halfedge(*(mesh.faces_begin() + fid)))){
                    auto source = mesh.source(hed), target = mesh.target(hed);
                    Base::Point p(mesh.points()[source]), q(mesh.points()[target]);
                    float t = CGAL::squared_distance(v, Base::Segment(p, q));
                    if (Base::floatCmp(t - tmp_face_dis) > 0){
                        tmp_face_dis = t;
                    }
                }
                if (Base::floatCmp(tmp_face_dis - dis) < 0){
                    dis = tmp_face_dis;
                }
            }
            gama[vid] = w_min * sqrt(dis) / w_max / 7;
        }
        return gama;
    }

    pair<float, int>
    placeSteinerPointsJACM(const Base::Mesh &mesh, const float &eps,
                           const vector<float> &gama,
                           map<unsigned, vector<unsigned> > &edge_bisector_map,
                           map<unsigned, vector<unsigned> > &bisector_point_map,
                           map<unsigned, unsigned> &point_face_map,
                           map<unsigned, Base::Point> &point_location_map,
                           map<unsigned, vector<unsigned>> &face_point_map){
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        unsigned num_vertices = mesh.num_vertices();
        unsigned bisector_num = 0;
        for (auto  fd: mesh.faces()){
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))){
                unsigned eid = mesh.edge(hed).idx(),
                        eid2 = mesh.edge(mesh.prev(hed)).idx();
                vector<Base::Point> p(3);
                vector<unsigned> pid(3);
                for (auto i = 0; i != 3; hed = mesh.next(hed), i++){
                    p[i] = mesh.points()[mesh.source(hed)];
                    pid[i] = mesh.source(hed).idx();
                    if (point_location_map.find(pid[i]) == point_location_map.end()){
                        point_location_map[pid[i]] = p[i];
                    }
                    face_point_map[fd.idx()].push_back(pid[i]);
                }
                vector<unsigned> cur_bisector = {};
                cur_bisector.push_back(pid[0]);

                float len1 = sqrt(CGAL::squared_distance(p[1], p[0])),
                        len2 = sqrt(CGAL::squared_distance(p[2], p[0]));
                Base::Point p_end(p[1] + Base::Vector(p[2] - p[1]) * len1 / (len1 + len2));
                Base::Vector vec_bisector(p_end - p[0]);
                Base::Point aux1 = p[0] + Base::Vector(p[1] - p[0]) * gama[pid[0]] / len1;
                Base::Point aux2 = p[0] + Base::Vector(p[2] - p[0]) * gama[pid[0]] / len2;
                Base::Point bisector_p0 = aux1 + 0.5 * Base::Vector(aux2 - aux1);
                float angle = Base::PI * CGAL::approximate_angle(p[1], p[0], p[2]) / 180;
                float sin_val = sin(angle * 0.5);
                float limit_distance = sqrt(CGAL::squared_distance(p[0], p_end));
                float cur_distance = sqrt(CGAL::squared_distance(p[0], bisector_p0));

                float num_Steiner_points = 1.61 / sin(angle) * log(2 * limit_distance / gama[pid[0]]);
                num_Steiner_points *= 1 / sqrt(eps) * log(2 / eps);
                bool upper_flag = 0;
                if (Base::floatCmp(num_Steiner_points) <= 0 || Base::floatCmp(num_Steiner_points - 25) > 0){
                    cur_distance = limit_distance / 25;
                    upper_flag = 1;
                }

                while (Base::floatCmp(cur_distance - limit_distance) < 0) {
                    Base::Point bisector_p = p[0] + cur_distance / limit_distance * vec_bisector;
                    point_location_map[num_vertices] = bisector_p;
                    point_face_map[num_vertices] = fd.idx();
                    face_point_map[fd.idx()].push_back(num_vertices);
                    cur_bisector.push_back(num_vertices++);
                    float distance_delta = sin_val * sqrt(0.5 * eps) * cur_distance;
                    if (upper_flag){
                        distance_delta = limit_distance / 50;
                    }
                    cur_distance += distance_delta;
                }
                bisector_point_map[bisector_num] = cur_bisector;
                edge_bisector_map[eid].push_back(bisector_num);
                edge_bisector_map[eid2].push_back(bisector_num++);
            }
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return make_pair(static_cast<float>(duration.count()), num_vertices);
    }

    pair<float, unsigned> placeSteinerPointsFixed(const Base::Mesh &mesh, const int &sp_num,
                                                  map<unsigned, vector<unsigned>> &edge_bisector_map,
                                                  map<unsigned, vector<unsigned>> &bisector_point_map,
                                                  map<unsigned, unsigned> &point_face_map,
                                                  map<unsigned, Base::Point> &point_location_map,
                                                  map<unsigned, vector<unsigned>> &face_point_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        unsigned num_vertices = mesh.num_vertices();
        unsigned bisector_num = 0;
        for (auto fd: mesh.faces()) {
            for (auto hed:mesh.halfedges_around_face(mesh.halfedge(fd))) {
                unsigned eid = mesh.edge(hed).idx(),
                        eid2 = mesh.edge(mesh.prev(hed)).idx();
                vector<Base::Point> p(3);
                vector<unsigned> pid(3);
                for (auto i = 0; i != 3; hed = mesh.next(hed), i++) {
                    p[i] = mesh.points()[mesh.source(hed)];
                    pid[i] = mesh.source(hed).idx();
                    if (point_location_map.find(pid[i]) == point_location_map.end()) {
                        point_location_map[pid[i]] = p[i];
                    }
                    face_point_map[fd.idx()].push_back(pid[i]);
                }

                vector<unsigned> cur_bisector_points = {};
                cur_bisector_points.push_back(pid[0]);

                float len1 = sqrt(CGAL::squared_distance(p[1], p[0])),
                        len2 = sqrt(CGAL::squared_distance(p[2], p[0]));
                Base::Point p_end(p[1] + Base::Vector(p[2] - p[1]) * len1 / (len1 + len2));
                Base::Vector vec_bisector(p_end - p[0]);
                float limit_distance = sqrt(CGAL::squared_distance(p[0], p_end));
                float cur_distance = 0;
                unsigned k = 0;

                while (Base::floatCmp(cur_distance - limit_distance) < 0 && Base::floatCmp(limit_distance / sp_num) > 0) {
                    Base::Point bisector_p = p[0] + static_cast<float>(++k) / sp_num * vec_bisector;
                    point_location_map[num_vertices] = bisector_p;
                    point_face_map[num_vertices] = fd.idx();
                    face_point_map[fd.idx()].push_back(num_vertices);
                    cur_bisector_points.push_back(num_vertices++);
                    cur_distance += limit_distance / sp_num;
                }
                bisector_point_map[bisector_num] = cur_bisector_points;
                edge_bisector_map[eid].push_back(bisector_num);
                edge_bisector_map[eid2].push_back(bisector_num++);
            }
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return make_pair(static_cast<float>(duration.count()), num_vertices);
    }

    float constructBaseGraph(const Base::Mesh &mesh,
                             const vector<float> &face_weight,
                             const unsigned &num_vertices,
                             map<unsigned, vector<unsigned>> &edge_bisector_map,
                             map<unsigned, vector<unsigned>> &bisector_point_map,
                             map<unsigned, unsigned> &point_face_map,
                             map<unsigned, Base::Point> &point_location_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        bool g_flag = true;
        if (kSkip::my_base_graph.num_V != num_vertices){
            g_flag = false;
            kSkip::my_base_graph.init(num_vertices);
        }
        vector<pair<unsigned, unsigned> > base_graph_edges = {};
        vector<float> base_graph_weights = {};
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
                            float dis = Base::distanceSnell(mesh, face_weight,
                                                            point_location_map[pid_1], fid_1,
                                                            point_location_map[pid_2], fid_2,
                                                            common_eid);
                            base_graph_edges.emplace_back(pid_1, pid_2);
                            base_graph_weights.push_back(dis);
                            if (!g_flag) kSkip::my_base_graph.addEdge(pid_1, pid_2, dis);
                        }
                    }
                }
            }
        }
        //boost graph.
        base_graph = Base::Graph(begin(base_graph_edges), end(base_graph_edges), begin(base_graph_weights),
                                 num_vertices);


        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<float>(duration.count());
    }

    class PartitionTreeNode {
    public:
        int node_id{};
        int center_idx;
        float radius;
        int level;
        PartitionTreeNode *parent;
        map<int, PartitionTreeNode *> children;
        map<int, bool> cover_set;

        bool operator()(PartitionTreeNode *a, PartitionTreeNode *b) {
            return a->center_idx < b->center_idx;
        }

        PartitionTreeNode();

        PartitionTreeNode(int center_idx_value, float radius_value, int node_id_value);

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

    inline PartitionTreeNode::PartitionTreeNode(int center_idx_value, float radius_value, int node_id_value) {
        parent = nullptr;
        children.clear();
        center_idx = center_idx_value;
        radius = radius_value;
        node_id = node_id_value;
    }

    inline PartitionTreeNode::~PartitionTreeNode() {
    }

    typedef pair<float, pair<PartitionTreeNode *, PartitionTreeNode *> > nodePairWithValue;
    typedef pair<PartitionTreeNode*, PartitionTreeNode*> nodePair;
    map<pair<int, int>, float> enhanced_edges; // store the enhanced edges

    class PartitionTree {
    public:
        PartitionTreeNode *root;
        int max_level;
        int num_tree_nodes;
        vector<vector<PartitionTreeNode *> > level_nodes;
        bool stop_flag;

        PartitionTree();
        ~PartitionTree();

        float constructTree(int num_vertices);

        float generateNodePairSet(const float &eps, const int &type, const unsigned &sp_num, set<nodePair> &node_pair_set);

        static void getPathToRoot(PartitionTreeNode* node, vector<PartitionTreeNode*> &path);

        void rootNodeSelect(const vector<unsigned> &pid_list);

        int buildLevel(const vector<unsigned> &pid_list, const float &l);

        float constructTree(ofstream &fout, const vector<unsigned> &pid_list, const float &l);

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

    inline float wellSeparatedValue(float center_distance, float eps, int type,
                                    PartitionTreeNode* node1, PartitionTreeNode* node2){
        if (!Base::floatCmp(node1->radius) && !Base::floatCmp(node2->radius)) return 0;
        assert(Base::floatCmp(eps));
        float rhs = (2 / eps + 2) * max(node1->radius, node2->radius);
        if (!type){
            rhs = (8 / eps + 4) * max(node1->radius, node2->radius);
        }
        return rhs;
    }

    float PartitionTree::generateNodePairSet(const float &eps, const int &type, const unsigned &sp_num, set<nodePair> &node_pair_set) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        priority_queue<nodePairWithValue , vector<nodePairWithValue>, greater<nodePairWithValue> > q;
        float center_distance = 0.0;
        float ws_value = -1.0;

        ws_value = wellSeparatedValue(center_distance, eps, type, root, root);

        q.push(make_pair(0 - ws_value, make_pair(root, root)));
        while (Base::floatCmp(q.top().first) < 0){
            nodePairWithValue cur_node_pair;
            cur_node_pair = q.top(); q.pop();
            auto l_node = cur_node_pair.second.first;
            auto r_node = cur_node_pair.second.second;
            //  Split the left node
            if (Base::floatCmp(l_node->radius - r_node->radius) > 0 ||
                !Base::floatCmp(l_node->radius - r_node->radius)
                && l_node->center_idx <= r_node->center_idx){
                for (auto &child: l_node->children){
                    auto child_idx = child.first;
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
                    auto child_idx = child.first;
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
        return static_cast<float>(duration.count());
    }

    // Construct partition tree.
    float PartitionTree::constructTree(ofstream &fout, const vector<unsigned> &pid_list, const float &l) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        auto num_vertices = pid_list.size();
        rootNodeSelect(pid_list);
        while (!stop_flag){
            unsigned cur_num_level_nodes = buildLevel(pid_list, l);
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

        return static_cast<float>(duration.count());
    }

    void PartitionTree::rootNodeSelect(const vector<unsigned> &pid_list) {
        auto num_vertices = pid_list.size();
        int root_index = pid_list[rand() % num_vertices];   // random select one point
//        vector<float> d(distance_map[root_index]);
        vector<float> d;
        kSkip::bounded_dijkstra(kSkip::my_base_graph, root_index, Base::unreachable, d);
        float radius = -1.0;
//        for (auto &dis: d){
//            if (Base::floatCmp(radius - distance) < 0){
//                radius = distance;
//            }
//            if (!Base::floatCmp(radius - Base::unreachable)){
//                cout << "ERROR: Unconnected Steiner Points found!" << endl;
//            }
//        }
//        cout << "num vertices = " << num_vertices << " d.size() = " << d.size() << endl;
        for (auto i = 0; i < num_vertices; i++){
            if (pid_list[i] == root_index) continue;
            float distance = d[pid_list[i]];
//            cout << "[dis] = " << distance << endl;
            if (Base::floatCmp(radius - distance) < 0){
                radius = distance;
            }
            if (!Base::floatCmp(radius - Base::unreachable)){
                int x; cin >> x;
                cout << "ERROR: Unconnected Steiner Points found!" << endl;
                cout << "ERROR Steiner Point ID = " << i << endl;
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

    int PartitionTree::buildLevel(const vector<unsigned> &pid_list, const float &l){
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
        float cur_level_radius = level_nodes[max_level][0]->radius * 0.5;
        cout << "[Partition Tree] Current level radius is: " << cur_level_radius << endl;
        for (auto &last_level_node: level_nodes[max_level]){
            int vertex_idx = last_level_node->center_idx;
            PartitionTreeNode* cur_node = new PartitionTreeNode(vertex_idx, cur_level_radius, num_tree_nodes++);
//            vector<float> d = distance_map[vertex_idx];
            vector<float> d; kSkip::bounded_dijkstra(kSkip::my_base_graph, vertex_idx, l * cur_level_radius, d);
            for (auto it = vertex_to_be_covered.begin(); it != vertex_to_be_covered.end(); ){
                int id = *it;
                assert(Base::floatCmp(d[id]) >= 0);
                if (Base::floatCmp(d[id] - cur_level_radius) <= 0){
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
                if (Base::floatCmp(d[id] - l * cur_level_radius) < 0){
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
//                vector<float> d = distance_map[vertex_idx];
                vector<float> d; kSkip::bounded_dijkstra(kSkip::my_base_graph, vertex_idx, l * cur_level_radius, d);
                for (auto it = vertex_to_be_covered.begin(); it != vertex_to_be_covered.end(); ){
                    int id = *it;
                    if (Base::floatCmp(d[id] - cur_level_radius) <= 0){
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
                    if (Base::floatCmp(d[id] - l * cur_level_radius) < 0){
                        if (!(enhanced_edges.find(make_pair(vertex_idx, id)) == enhanced_edges.end()))continue;
                        enhanced_edges[make_pair(vertex_idx, id)] = d[id];
                    }
                }
                assert(!cur_node->cover_set.empty());  //  Each node cover itself at least
                int father_idx = -1;
                PartitionTreeNode* father;
                for (auto &last_level_node: level_nodes[max_level]){
                    int fa_id = last_level_node->center_idx;
                    assert(Base::floatCmp(d[father_idx]) >= 0);
                    assert(Base::floatCmp(d[fa_id]) >= 0);
                    if (father_idx < 0 || Base::floatCmp(d[fa_id] - d[father_idx]) < 0){
                        father_idx = fa_id;
                        father = last_level_node;
                    }
                }
                assert(Base::floatCmp(d[father_idx] - cur_level_radius * 2.0) <= 0);
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

    float distanceQueryEfficient(set<nodePair> &node_pairs,
                                 vector<PartitionTreeNode*> &As,
                                 vector<PartitionTreeNode*> &At){
        float approximate_distance = 0.0;
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

    float distanceQueryA2A(const Base::Point &query_source, const unsigned &fid_s,
                           const Base::Point &query_target, const unsigned &fid_t,
                           const PartitionTree &tree,
                           map<unsigned, vector<unsigned> > &face_point_map,
                           map<unsigned, Base::Point> &point_location_map,
                           set<nodePair> &node_pairs
    ){
        float ret_dis = Base::unreachable;
        auto leaf_nodes = tree.level_nodes[tree.max_level];
        for (auto pid1: face_point_map[fid_s]){
            for (auto pid2: face_point_map[fid_t]){
                vector<PartitionTreeNode*> As, At;
                tree.getPathToRoot(leaf_nodes[pid1], As);
                tree.getPathToRoot(leaf_nodes[pid2], At);
                float t_dis = sqrt(CGAL::squared_distance(query_source, point_location_map[pid1])) +
                              distanceQueryEfficient(node_pairs, As, At) +
                              sqrt(CGAL::squared_distance(point_location_map[pid2], query_target));
                if (Base::floatCmp(t_dis - ret_dis) < 0) ret_dis = t_dis;
            }
        }
        if (!Base::floatCmp(ret_dis - Base::unreachable)){
            cout << "ERROR: A2A distance not found!" << endl;
        }
        return ret_dis;
    }

}

#endif //WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H