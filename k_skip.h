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

        vector<int> head;
        vector<Edge> edges;

        void init(int v) {
            num_V = v;
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
//                    k_skip_graph.addEdge(k_cover_vertex_id[t], k_cover_vertex_id[query_s], d[t]);

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

}

#endif //WEIGHTED_DISTANCE_ORACLE_K_SKIP_H
