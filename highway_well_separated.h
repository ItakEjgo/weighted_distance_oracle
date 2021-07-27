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
        int from, to, next;
        Edge(int u, int v, double len, int nxt): from(u), to(v), next(nxt), w(len){}
    };

    struct QNode{
        int p;
        double dis;
        bool operator < (const QNode& x) const{
            return Base::doubleCmp(dis - x.dis) < 0 || !Base::doubleCmp(dis - x.dis) && p < x.p;
        }
        QNode(int pid, double d_val): p(pid), dis(d_val){}
    };

    struct HighwayGraph{
        int num_V;
        int num_E;
        vector<int> head;
        vector<Edge> edges;

        void init(int v){
            num_V = v;
            num_E = 0;
            head.resize(v, 0);
            edges.clear();
            edges.emplace_back(-1, -1, -1.0, -1); // stop edge.
            num_E++;
        }

        void addEdge(int u, int v, double w){
            edges.emplace_back(u, v, w, head[u]);
            head[u] = num_E;
            num_E++;
        }

        static inline double wellSeparatedRadius(const double dis, const double eps){
            return dis / (2 / eps + 2);
        }

        void highwayPropagate(int s, const double eps){
            vector<double> d(num_V, Base::unreachable);
            vector<bool> vis(num_V, false), vis_cir(num_V, false);
            d[s] = 0;
            priority_queue<QNode> q = {};
            priority_queue<QNode> circle_heap = {};
            q.push(QNode(s, d[s]));
            int cir_label = 0;
            while (!q.empty()){
                QNode f = q.top(); q.pop();
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
                if (f.p != s) circle_heap.push(f);
                double max_radius = wellSeparatedRadius(d[f.p], eps);
                bool new_cir_level_flag = false;
                while (!circle_heap.empty()){
                    QNode f_cir = circle_heap.top();
                    if (vis_cir[f_cir.p]) continue;
                    if (Base::doubleCmp(d[f_cir.p] - max_radius) <= 0){
                        if (!new_cir_level_flag) new_cir_level_flag = true;
                        addEdge(s, f_cir.p, d[f_cir.p]);
                        addEdge(f_cir.p, s, d[f_cir.p]);
                        circle_heap.pop();
                    }
                }
            }
        }
    };
}


#endif //WEIGHTED_DISTANCE_ORACLE_HIGHWAY_WELL_SEPARATED_H
