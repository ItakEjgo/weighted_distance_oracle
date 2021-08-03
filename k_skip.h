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

        void addVertex(){
            num_V++;
            head.emplace_back(0);
        }
    };

    // k-hop-vector (distance from s, father id)
    vector<int> generateKHopVector(Graph &g, int s, vector<double> &d, vector<int> &fa){
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
            if (f.hop_cnt < K - 1){
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

    set<int> adaptiveSampling(Graph &g){
        set<int> k_cover_V = {};
        vector<int> vid(g.num_V);
        for (auto i = 0; i < g.num_V; i++) vid[i] = i;
        auto seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(vid.begin(), vid.end(), default_random_engine(seed));
        for (auto s: vid){
            vector<double> d;
            vector<int> fa;
            auto k_hop_vector = generateKHopVector(g, s, d, fa);
            bool add_s_flag = false;
            for (auto t: k_hop_vector){
                if (add_s_flag) break;
                int hop_cnt = 0;
                bool appear_flag = false;
                while (t != -1){
                    if (!appear_flag && k_cover_V.find(t) != k_cover_V.end()){
                        appear_flag = true;
                    }
                    if (!appear_flag && hop_cnt >= K - 1){
                        add_s_flag = true;
                    }
                    t = fa[t];
                    hop_cnt++;
                }
            }
            if (add_s_flag){
                k_cover_V.insert(s);
            }
        }
        return k_cover_V;
    }

}

#endif //WEIGHTED_DISTANCE_ORACLE_K_SKIP_H
