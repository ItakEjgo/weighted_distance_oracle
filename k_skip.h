//
// Created by huang on 2021/8/3.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_K_SKIP_H
#define WEIGHTED_DISTANCE_ORACLE_K_SKIP_H
#include "base.h"

extern int K;

namespace kSkip{
    using namespace std;

    struct Edge{
        double w;
        int from, to, next;
        Edge(int u, int v, double len, int nxt): from(u), to(v),  w(len), next(nxt){}
    };

    struct QNode{
        int p;
        int hop_cnt;
        double dis;
        bool operator < (const QNode& x) const{
            return Base::doubleCmp(dis - x.dis) > 0 || !Base::doubleCmp(dis - x.dis) && p < x.p;
        }
        QNode(int pid, double d_val, int hops = 0): p(pid), dis(d_val), hop_cnt(hops){}
    };

    struct Graph {
        int num_V;
        int num_E;
        int corner_vertex_flag;

        vector<int> head;
        vector<Edge> edges;

        void init(int v) {
            num_V = v;
            corner_vertex_flag = v;
            head.resize(v, 0);
            edges.clear();
            num_E = 1;
            edges.emplace_back(-1, -1, -1.0, -1); // stop edge.
        }

        void addEdge(int u, int v, double w) {
            edges.emplace_back(u, v, w, head[u]);
            head[u] = num_E;
            num_E++;
        }

        int addVertex(){
            head.emplace_back(0);
            return num_V++;
        }

        bool isCorner(int vid) {
            return vid < corner_vertex_flag;
        }
    };

    Graph my_base_graph;

    // k-hop-vector
    vector<int> generateKHopVector(Graph &g, int s, vector<double> &d, vector<int> &fa, int max_hops = K){
        vector<int> k_hop_vector = {};
        vector<bool> vis(g.num_V, false);
        d.resize(g.num_V, Base::unreachable);
        fa.resize(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        while (!q.empty()){
            QNode f = q.top(); q.pop();
            if (vis[f.p]) continue;
            if (f.hop_cnt < max_hops - 1){
                for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                    int v = g.edges[eid].to;
                    double w = g.edges[eid].w;
                    if (Base::doubleCmp(d[f.p] + w - d[v]) < 0){
                        d[v] = d[f.p] + w;
                        fa[v] = f.p;
                        q.push(QNode(v, d[v], f.hop_cnt + 1));
                    }
                }
            }
            vis[f.p] = true;
            k_hop_vector.emplace_back(f.p);
        }
        return k_hop_vector;
    }

    //generate vLambda
    set<int> generateVLambda(Graph &g, int s, map<int, int> &v_id, map<int, int> &ori_v_id, int max_hops = K){
        set<int> v_lambda = {};
        queue<pair<int, int> > q;
        vector<bool> vis(g.num_V, false);
        vis[s] = true;
        q.emplace(s, 0);
        while (!q.empty()){
            pair<int, int> f = q.front(); q.pop();
            if (f.second > max_hops) break;
            for (int eid = g.head[f.first]; eid; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                if (!vis[v]) {
                    vis[v] = true;
                    q.emplace(v, f.second + 1);
                }
            }
            v_lambda.insert(f.first);
        }
        int cnt = 0;
        for (auto &val: v_lambda){
            v_id[val] = cnt;
            ori_v_id[cnt++] = val;
        }

        return v_lambda;
    }

    vector<int> generateShortestPathTree(Graph &g, set<int> &v_lambda, int s, vector<double> &d, vector<int> &fa){
        vector<int> k_hop_vector = {};
        vector<bool> vis(g.num_V, false);
        d.resize(g.num_V, Base::unreachable);
        fa.resize(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        while (!q.empty()){
            QNode f = q.top(); q.pop();
            if (vis[f.p]) continue;
            for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                if (v_lambda.find(v) == v_lambda.end()) continue;
                double w = g.edges[eid].w;
                if (Base::doubleCmp(d[f.p] + w - d[v]) < 0){
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
            k_hop_vector.emplace_back(f.p);
        }
        return k_hop_vector;
    }

    set<int> adaptiveSampling(Graph &g){
        set<int> k_cover_V = {};
        vector<int> vid(g.num_V);
        for (auto i = 0; i < g.num_V; i++) vid[i] = i;
        auto seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(vid.begin(), vid.end(), default_random_engine(seed));
        for (auto s: vid){
            map<int, int> ori_v_id, v_id;
            auto v_lambda = generateVLambda(g, s, v_id, ori_v_id, K - 1);
//            cout << "s = " << s << " size = " << v_lambda.size() << endl;
            vector<double> d;
            vector<int> fa;
            auto k_hop_vector = generateShortestPathTree(g, v_lambda, s, d, fa);
//            cout << "v-lambda size = " << v_lambda.size() << "  k_hop_vector size = " << k_hop_vector.size() << endl;

            bool is_not_covered = false;
            for (auto t: k_hop_vector){
                if (is_not_covered) break;
                bool appear_flag = false;
                stack<int> path = {};
                while (t != -1){
                    path.emplace(t);
                    t = fa[t];
                }
                if (path.size() >= K - 1){
                    while (!path.empty()){
                        int top = path.top(); path.pop();
                        if (k_cover_V.find(top) != k_cover_V.end()){
                            appear_flag = true;
                            break;
                        }
                    }
                    if (!appear_flag) is_not_covered = true;
                }
            }
            if (is_not_covered){
                k_cover_V.insert(s);
            }
        }
        return k_cover_V;
    }

    map<int, int> generateKCoverVertexId(set<int> &k_cover_V){
        map<int, int> ret = {};
        int pos = 0;
        for (auto s: k_cover_V){
            ret[s] = pos++;
        }
        return ret;
    }

    Graph computeKSkipGraph(Graph &g, set<int> &k_cover_V, map<int, int> &k_cover_vertex_id){
        Graph k_skip_graph;
        k_skip_graph.init(static_cast<int>(k_cover_V.size()));
        map<pair<int, int>, int> edge_count;
        for (auto s: k_cover_V){
            vector<double> d; vector<int> fa;
//            auto k_hop_vector = generateKHopVector(g, s, d, fa, K + 1);
            map<int, int> ori_v_id, v_id;
            auto v_lambda = generateVLambda(g, s, v_id, ori_v_id, K);

            auto k_hop_vector = generateShortestPathTree(g, v_lambda, s, d, fa);
//            cout << "v-lambda size = " << v_lambda.size() << "  k_hop_vector size = " << k_hop_vector.size() << endl;
            for (auto t: k_hop_vector){
                if (k_cover_V.find(t) != k_cover_V.end()){
                    bool first_appear = true;
                    int cur = t;
                    stack<int> path = {};
                    while (cur != -1){
                        path.emplace(cur);
                        cur = fa[cur];
                    }
                    while (!path.empty()){
                        int top = path.top(); path.pop();
                        if (top != s && top != t && k_cover_V.find(top) != k_cover_V.end()){
                            first_appear = false;
                            break;
                        }
                    }
                    if (!first_appear) continue;
//                    cout << "Edge: " << s + 1 << " " << t + 1 << " " <<  d[t] << endl;
                    k_skip_graph.addEdge(k_cover_vertex_id[s], k_cover_vertex_id[t], d[t]);
//                    k_skip_graph.addEdge(k_cover_vertex_id[t], k_cover_vertex_id[s], d[t]);

                    edge_count[make_pair(k_cover_vertex_id[s], k_cover_vertex_id[t])]++;
//                    k_skip_graph.addEdge(k_cover_vertex_id[t], k_cover_vertex_id[s], d[t]);
                }
            }
        }
        for (auto & it : edge_count){
            if (it.second > 1){
                cout << "edge: " << it.first.first << "," << it.first.second << " appears " << it.second << "times" << endl;
            }
        }
        return k_skip_graph;
    }

    // s,t dijkstra query
    pair<double, double> dijkstra(Graph &g, int s, int t){
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        vector<bool> vis(g.num_V, false);
        vector<double> d(g.num_V, Base::unreachable);
        vector<int> fa(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        while (!q.empty()){
            QNode f = q.top(); q.pop();
            if (f.p == t) break;
            if (vis[f.p]) continue;
            for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                double w = g.edges[eid].w;
//                if (Base::doubleCmp(w) < 0){
//                    cout << "w < 0: " << f.p << " " << v << " " << w << endl;
//                }
                if (Base::doubleCmp(d[f.p] + w - d[v]) < 0){
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
        }
        assert(Base::doubleCmp(d[t] - Base::unreachable) < 0);

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return make_pair(d[t], static_cast<double>(duration.count()));
    }

    double bounded_dijkstra(Graph &g, int s, double bound, vector<double> &d){
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        vector<bool> vis(g.num_V, false);
        d.resize(g.num_V, Base::unreachable);
        vector<int> fa(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        while (!q.empty()){
            QNode f = q.top(); q.pop();
            if (Base::doubleCmp(d[f.p] - bound) > 0) break;
            if (vis[f.p]) continue;
            for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                double w = g.edges[eid].w;
//                if (Base::doubleCmp(w) < 0){
//                    cout << "w < 0: " << f.p << " " << v << " " << w << endl;
//                }
                if (Base::doubleCmp(d[f.p] + w - d[v]) < 0){
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    double covered_dijkstra(Graph &g, int s, set<int> &cover_id, vector<double> &d){
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        vector<bool> vis(g.num_V, false);
        d.resize(g.num_V, Base::unreachable);
        vector<int> fa(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        int cnt = cover_id.size();
        while (!q.empty()){
            QNode f = q.top(); q.pop();
            if (!cnt) break;
            if (vis[f.p]) continue;
            for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                double w = g.edges[eid].w;
//                if (Base::doubleCmp(w) < 0){
//                    cout << "w < 0: " << f.p << " " << v << " " << w << endl;
//                }
                if (Base::doubleCmp(d[f.p] + w - d[v]) < 0){
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
            if (cover_id.find(f.p) != cover_id.end()) cnt--;
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    double queryGraphA2A(Graph g, Base::Point s, int fid_s, Base::Point t, int fid_t,
                         map<int, vector<int> > &face_point_map, map<int, Base::Point> &point_location_map){
        //  add edges from s to its neighbor Steiner points
        int sid = g.addVertex();
        for (auto pid: face_point_map[fid_s]){
            auto pt = point_location_map[pid];
            double dis = sqrt(CGAL::squared_distance(s, pt));
            g.addEdge(sid, pid, dis);
        }
        //  add edges from t's neighbor Steiner points to t
        int tid = g.addVertex();
        for (auto pid: face_point_map[fid_t]){
            auto pt = point_location_map[pid];
            double dis = sqrt(CGAL::squared_distance(t, pt));
            g.addEdge(pid, tid, dis);
        }
        return dijkstra(g, sid, tid).first;
    }

    pair<double, double> queryKSkipGraph(Graph &g, Graph k_skip_graph, set<int> &k_cover_V, map<int, int> &k_cover_vertex_id, int query_s, int query_t){
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        // add s to k_skip_graph
        int sid, tid;
        if (k_cover_V.find(query_s) == k_cover_V.end()){
            sid = k_skip_graph.addVertex();
            vector<double> d; vector<int> fa;
            map<int, int> ori_v_id, v_id;
            auto v_lambda = generateVLambda(g, query_s, v_id, ori_v_id, K);
            auto k_hop_vector = generateShortestPathTree(g, v_lambda, query_s, d, fa);

            for (auto t: k_hop_vector){
                if (k_cover_V.find(t) != k_cover_V.end()){
                    bool first_appear = true;
                    int cur = t;
                    stack<int> path = {};
                    while (cur != -1){
                        path.emplace(cur);
                        cur = fa[cur];
                    }
                    while (!path.empty()){
                        int top = path.top(); path.pop();
                        if (top != query_s && top != t && k_cover_V.find(top) != k_cover_V.end()){
                            first_appear = false;
                            break;
                        }
                    }
                    if (!first_appear) continue;
                    k_skip_graph.addEdge(sid, k_cover_vertex_id[t], d[t]);
                }
            }
        }
        else{
//            cout << "s in the k-skip graph." << endl;
            sid = k_cover_vertex_id[query_s];
        }
        if (k_cover_V.find(query_t) == k_cover_V.end()){
            tid = k_skip_graph.addVertex();
            vector<double> d; vector<int> fa;
            map<int, int> ori_v_id, v_id;
            auto v_lambda = generateVLambda(g, query_t, v_id, ori_v_id, K);
            auto k_hop_vector = generateShortestPathTree(g, v_lambda, query_t, d, fa);

            for (auto t: k_hop_vector){
                if (k_cover_V.find(t) != k_cover_V.end()){
                    bool first_appear = true;
                    int cur = t;
                    stack<int> path = {};
                    while (cur != -1){
                        path.emplace(cur);
                        cur = fa[cur];
                    }
                    while (!path.empty()){
                        int top = path.top(); path.pop();
                        if (top != query_t && top != t && k_cover_V.find(top) != k_cover_V.end()){
                            first_appear = false;
                            break;
                        }
                    }
                    if (!first_appear) continue;
                    k_skip_graph.addEdge(k_cover_vertex_id[t], tid, d[t]);
//                    k_skip_graph.addEdge(k_cover_vertex_id[query_t], k_cover_vertex_id[t], d[t]);

                }
            }
        }
        else{
//            cout << "t in the k-skip graph." << endl;
            tid = k_cover_vertex_id[query_t];
        }

        double ret_dis = dijkstra(k_skip_graph, sid, tid).first;

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

        return make_pair(ret_dis, static_cast<double>(duration.count()));
    }

    double constructGraph(Graph &g, Base::Mesh &mesh,
                                 vector<double> &face_weight,
                                 int num_vertices,
                                 map<int, vector<int>> &edge_bisector_map,
                                 map<int, vector<int>> &bisector_point_map,
                                 map<int, int> &point_face_map,
                                 map<int, Base::Point> &point_location_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

//        vector<pair<int, int> > base_graph_edges = {};
//        vector<double> base_graph_weights = {};
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
//                            base_graph_edges.emplace_back(pid_1, pid_2);
                            g.addEdge(pid_1, pid_2, dis);
//                            base_graph_weights.push_back(dis);
                        }
                    }
                }
            }
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    void constructMeshGraph(Base::Mesh &mesh, Graph &g){
        g.init(mesh.num_vertices());
        g.head.resize(mesh.num_vertices(), 0);
        for (auto &fd: mesh.faces()){
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))){
                auto sid = mesh.source(hed);
                auto tid = mesh.target(hed);
                double w = sqrt(CGAL::squared_distance(mesh.points()[sid], mesh.points()[tid]));
                g.addEdge(sid.idx(), tid.idx(), w);
                g.addEdge(tid.idx(), sid.idx(), w);
            }
        }
    }

    void getAdjacentFaces(Base::Mesh &m, Graph &g,
                          map<int, int> &cut_vertex_halfedge, int u, vector<int> &adj_face_ids){
        adj_face_ids.clear();
        if (g.isCorner(u)){
            CGAL::Face_around_target_circulator<Base::Mesh> vbegin(m.halfedge(CGAL::SM_Vertex_index(u)), m), done(vbegin);
            do{
                auto fd = *vbegin++;
                adj_face_ids.push_back(fd.idx());
            }while (vbegin != done);
        }
        else{
            int halfedge_id = cut_vertex_halfedge[u];
            // cout << "halfedge id = " << halfedge_id << endl;
            auto hed = CGAL::SM_Halfedge_index(halfedge_id);
            adj_face_ids.push_back(m.face(hed).idx());
            adj_face_ids.push_back(m.face(m.opposite(hed)).idx());
        }
        // cout << "adjacent_faces: ";
        // for (auto fid: adj_face_ids){
        //     cout << fid << " ";
        // }
        // cout << endl;
    }

    // VLDB'2015 implementation
    double computeDistanceBound(
            Base::Mesh m, Graph g,
            int K, Base::Point point_s, int fid_s, Base::Point point_t, int fid_t,
            double l_min
//            vector<double> &D,
//            pair<double, double> &bounds,
//            map<int, vector<int> > &edge_cut_vertex
    ){
        auto vds = m.add_vertex(point_s);
        int s = vds.idx();
        g.addVertex();
        auto fds = *(m.faces_begin() + fid_s);
        for (auto vd: m.vertices_around_face(m.halfedge(fds))){
            m.add_edge(vds, vd);
            double dis = sqrt(CGAL::squared_distance(point_s, m.points()[vd]));
            g.addEdge(s, vd.idx(), dis);
        }
        auto vdt = m.add_vertex(point_t);
        int t = vdt.idx();
        g.addVertex();
        auto fdt = *(m.faces_begin() + fid_t);
        for (auto vd: m.vertices_around_face(m.halfedge(fdt))){
            m.add_edge(vdt, vd);
            double dis = sqrt(CGAL::squared_distance(point_t, m.points()[vd]));
            g.addEdge(vd.idx(), t, dis);
        }


        map<int, vector<int> > edge_cut_vertex = {};
        map<int, int> cut_vertex_halfedge = {};
        int vertices_num = m.num_vertices();
        vector<double> D;
        D.resize(vertices_num, -1.0);   // -1.0 means unreachable
        D[s] = 0.0;
        vector<bool> vis;
        vis.resize(vertices_num, 0);
        priority_queue<QNode> q = {};
        q.push(QNode(s, D[s]));

        int cnt = 0;
        while (!q.empty()){
            QNode f = q.top();
            q.pop();
//             cout << f.u << " " << f.d_u << endl;
//             cout << "cnt = " << cnt++ << endl;
            // int w; cin >> w;
            if (t == f.p){
//                double lambda = 1.0 - 1.0 / K;
//                bounds.first = lambda * D[f.u];
//                bounds.second = D[f.u];
                return D[f.p];
            }
            if (vis[f.p]) continue;
            // S = S U {u};
            double d_ut = sqrt(CGAL::squared_distance(m.points()[CGAL::SM_Vertex_index(f.p)], m.points()[CGAL::SM_Vertex_index(t)]));

            auto u = CGAL::SM_Vertex_index(f.p);
            vector<int> adj_face_ids;
            getAdjacentFaces(m, g, cut_vertex_halfedge, f.p, adj_face_ids);
            for (auto fid: adj_face_ids){
                auto fd = CGAL::SM_Face_index(fid);
                // auto fd = *vbegin++;
                if (fd != Base::Mesh::null_face()){
                    vector<int> vids = {};
                    pair<int, int> opposite_edge = {};
                    int opposite_edge_id = -1;
                    // vector<pair<int, int> > adjacent_edges = {};
                    vector<int> L_id = {};
                    auto hed = m.halfedge(fd);
                    vector<pair<int, int> > L = {};
                    for (int i = 0; i != 3; i++){
                        auto t_u = m.points()[m.source(hed)];
                        auto t_v = m.points()[m.target(hed)];
                        int uid = m.source(hed).idx();
                        int vid = m.target(hed).idx();
                        Base::Segment seg(t_u, t_v);
                        vids.push_back(uid);
                        if (!seg.has_on(m.points()[u])){
                            L.push_back(make_pair(uid, vid));
                            L_id.push_back(hed);
                        }
                        hed = m.next(hed);
                    }
                    double theta_m = -1.0;
                    for (int i = 0; i != 3; i++){
                        auto t_u = CGAL::SM_Vertex_index(vids[i]),
                                t_v = CGAL::SM_Vertex_index(vids[(i + 1) % 3]),
                                t_w = CGAL::SM_Vertex_index(vids[(i + 2) % 3]);
                        double angle = CGAL::approximate_angle(m.points()[t_u], m.points()[t_v], m.points()[t_w]);
                        if (Base::doubleCmp(theta_m) < 0 || Base::doubleCmp(angle - theta_m) < 0){
                            theta_m = angle;
                        }
                    }
                    // cout << "theta_m = " << theta_m << endl;
                    // compute |O_iO_i+1|min and delta_I.
                    double OO_min = 0.5 * l_min * sqrt(2 * (1 - cos(theta_m / 180 * Base::PI)));
                    double delta_I = OO_min / K;
                    // cout << "O_min = " << OO_min << endl;
                    // cout << "delta_I = " << delta_I << endl;

                    for (auto x = 0; x != L.size(); x++){
                        auto e = L[x];
                        // cout << "e = " << e.first << " " << e.second << endl;
                        auto v_a = m.points()[CGAL::SM_Vertex_index(e.first)];
                        auto v_b = m.points()[CGAL::SM_Vertex_index(e.second)];
                        double d_et = min(
                                sqrt(CGAL::squared_distance(v_a, m.points()[CGAL::SM_Vertex_index(t)])),
                                sqrt(CGAL::squared_distance(v_b, m.points()[CGAL::SM_Vertex_index(t)]))
                        );
                        if (Base::doubleCmp(d_ut - d_et) >= 0){
                            int j = 1;
                            double e_len = sqrt(CGAL::squared_distance(v_a, v_b));
//                             cout << "num cut vertices = " << floor(e_len / delta_I) << endl;
                            while (j <= floor(e_len / delta_I)){
                                int v_c_id = -1;
                                decltype(v_a) v_c;
                                if (edge_cut_vertex.find(L_id[x]) == edge_cut_vertex.end() || edge_cut_vertex[L_id[x]].size() < j){
                                    v_c = v_a + (v_b - v_a) / e_len * j * delta_I; // cut-vertex coordinate
                                    m.add_vertex(v_c);  //  add the vertex to mesh
                                    v_c_id = g.addVertex(); //   vertex id
                                    assert(v_c != m.points()[u]);
                                    cut_vertex_halfedge[v_c_id] = L_id[x];
                                    D.push_back(-1.0);
                                    vis.push_back(0);
                                    assert(v_c_id == g.num_V - 1);
                                    if (!edge_cut_vertex[L_id[x]].size()){
                                        edge_cut_vertex[L_id[x]] = {};
                                    }
                                    edge_cut_vertex[L_id[x]].push_back(v_c_id);
                                }
                                else{
                                    v_c_id = edge_cut_vertex[L_id[x]][j - 1];
                                    assert(edge_cut_vertex[L_id[x]].size() >= j);
                                    v_c = m.points()[CGAL::SM_Vertex_index(v_c_id)];
                                }
                                double d_u_vc = sqrt(CGAL::squared_distance(m.points()[u], v_c));
                                if (Base::doubleCmp(d_u_vc - delta_I) >= 0){
                                    g.addEdge(f.p, v_c_id, d_u_vc);
                                }
                                // D[v_c_id] = D[f.u] + d_u_vc;
                                // q.insert(Qnode(v_c_id, D[v_c_id]));
                                j++;
                            }
                        }
                    }
                    if (!g.isCorner(f.p)){
                        for (auto id: vids){
                            auto v = CGAL::SM_Vertex_index(id);
                            double w = sqrt(CGAL::squared_distance(m.points()[u], m.points()[v]));
                            g.addEdge(f.p, v.idx(), w);
                        }
                    }
                }
            }
            for (auto eid = g.head[f.p]; eid > 0; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                double w = g.edges[eid].w;
                if (Base::doubleCmp(D[v]) < 0 || Base::doubleCmp(D[v] - D[f.p] - w) > 0){
                    D[v] = D[f.p] + w;
                    q.push(QNode(v, D[v]));
                }
            }
            vis[f.p] = 1;
        }

    }
}

#endif //WEIGHTED_DISTANCE_ORACLE_K_SKIP_H
