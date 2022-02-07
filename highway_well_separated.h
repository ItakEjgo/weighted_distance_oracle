//
// Created by huang on 2021/7/26.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_HIGHWAY_WELL_SEPARATED_H
#define WEIGHTED_DISTANCE_ORACLE_HIGHWAY_WELL_SEPARATED_H
#include "base.h"

namespace Highway{

    using namespace std;

    struct Edge{
        double w;
        int from, to, next, label;
        Edge(int u, int v, double len, int nxt, int lb = 0): from(u), to(v), next(nxt), w(len), label(lb){}
    };

    struct QNode{
        int p;
        double dis;
        bool operator < (const QNode& x) const{
//            return Base::doubleCmp(dis - x.dis) < 0 || !Base::doubleCmp(dis - x.dis) && p < x.p;
            return Base::doubleCmp(dis - x.dis) > 0 || !Base::doubleCmp(dis - x.dis) && p < x.p;
        }
        QNode(int pid, double d_val): p(pid), dis(d_val){}
    };

    struct HighwayGraph{
        int num_V;
        int num_E;
        int num_H;
        vector<int> head;
        vector<Edge> edges;

        vector<int> head_highway;
        vector<Edge> edges_highway;

        void init(int v){
            num_V = v;
            num_E = 0;
            num_H = 0;
            head.resize(v, 0);
            head_highway.resize(v, 0);
            edges.clear();
            edges_highway.clear();
            edges.emplace_back(-1, -1, -1.0, -1); // stop edge.
            edges_highway.emplace_back(-1, -1, -1.0, -1); // stop edge.
            num_E++;
            num_H++;
        }

        void addEdge(int u, int v, double w){
            edges.emplace_back(u, v, w, head[u]);
            head[u] = num_E;
            num_E++;
        }

        void addHighway(int u, int v, double w, int lb){
            edges_highway.emplace_back(u, v, w, head[u], lb);
            head_highway[u] = num_H;
            num_H++;
        }

        static inline double wellSeparatedRadius(const double dis, const double eps){
            return dis / (2 / eps + 2);
        }

        void highwayPropagate(int s, const double eps){ //  generate first-well-separated-highway
            vector<double> d(num_V, Base::unreachable);
            vector<bool> vis(num_V, false), vis_cir(num_V, false);
            d[s] = 0;
            priority_queue<QNode> q = {};
            priority_queue<QNode> circle_heap = {};
            q.push(QNode(s, d[s]));
            int cir_label = 1;
            bool cir_exist_flag = false;
            while (!q.empty()){
                QNode f = q.top(); q.pop();
                if (vis[f.p]) continue;

                // get the distance of q.top(), generate highways.
                if (f.p != s) circle_heap.push(f);
                double max_radius = wellSeparatedRadius(d[f.p], eps);
                bool new_cir_level_flag = false;
                while (!circle_heap.empty()){
                    QNode f_cir = circle_heap.top();
                    if (vis_cir[f_cir.p]) continue;
                    if (Base::doubleCmp(d[f_cir.p] - max_radius) <= 0){
                        if (!new_cir_level_flag) new_cir_level_flag = true;
                        vis_cir[f_cir.p] = true;
                        addHighway(s, f_cir.p, d[f_cir.p], cir_label);
                        addHighway(f_cir.p, s, d[f_cir.p], cir_label);
                        circle_heap.pop();
                        if (!cir_exist_flag) cir_exist_flag = true;
                    }
                    else{
                        break;
                    }
                }
                if (cir_exist_flag){
                    if (new_cir_level_flag){
                        addHighway(s, f.p, d[f.p], cir_label);
                        addHighway(f.p, s, d[f.p], cir_label);
                        cir_label++;
                    }
                    else{
                        addHighway(s, f.p, d[f.p], cir_label - 1);
                        addHighway(f.p, s, d[f.p], cir_label - 1);
                    }
                    vis[f.p] = true;
                    continue;
                }

                for (int eid = head[f.p]; eid; eid = edges[eid].next){
                    int v = edges[eid].to;
                    double w = edges[eid].w;
                    if (Base::doubleCmp(d[f.p] + w - d[v]) < 0){
                        d[v] = d[f.p] + w;
                        q.push(QNode(v, d[v]));
                    }
                }
                vis[f.p] = true;
            }
        }

        double queryHighway(int s, int t){
            vector<double> d(num_V, Base::unreachable);
            vector<bool> vis(num_V, false);
            d[s] = 0;
            priority_queue<QNode> q = {};
            q.push(QNode(s, d[s]));
            while (!q.empty()){
                QNode f = q.top(); q.pop();
                if (f.p == t) break;
                if (vis[f.p]) continue;

                for (int eid = head_highway[f.p]; eid; eid = edges_highway[eid].next){
                    int v = edges_highway[eid].to;
                    double w = edges_highway[eid].w;
                    if (Base::doubleCmp(d[f.p] + w - d[v]) < 0){
                        d[v] = d[f.p] + w;
                        q.push(QNode(v, d[v]));
                    }
                }

                vis[f.p] = true;
            }
            return d[t];
        }

        double queryOriginGraph(int s, int t){
            vector<double> d(num_V, Base::unreachable);
            vector<bool> vis(num_V, false);
            d[s] = 0;
            priority_queue<QNode> q = {};
            q.push(QNode(s, d[s]));
            while (!q.empty()){
                QNode f = q.top(); q.pop();
                if (f.p == t) break;
                if (vis[f.p]) continue;

                for (int eid = head[f.p]; eid; eid = edges[eid].next){
                    int v = edges[eid].to;
                    double w = edges[eid].w;
                    if (Base::doubleCmp(d[f.p] + w - d[v]) < 0){
                        d[v] = d[f.p] + w;
                        q.push(QNode(v, d[v]));
                    }
                }

                vis[f.p] = true;
            }
            return d[t];
        }

        // return the distance and whether is hit on Highway Graph
        pair<double, bool> distanceQuery(int s, int t){
            double ret = queryHighway(s, t);
            bool hit = false;
            if (!Base::doubleCmp(ret - Base::unreachable)){
                ret = queryOriginGraph(s, t);
            }
            else{
                hit = true;
            }
            if (!Base::doubleCmp(ret - Base::unreachable)) {
                cout << "ERROR! could not found distance on both Highway & OriginGraph" << endl;
            }
            return make_pair(ret, hit);
        }
    };

    double constructHighwayGraph(HighwayGraph &g, Base::Mesh &mesh,
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


#endif //WEIGHTED_DISTANCE_ORACLE_HIGHWAY_WELL_SEPARATED_H
