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
        int from, to, next, prev;
        Edge(int u, int v, double len, int nxt, int prv): from(u), to(v),  w(len), next(nxt), prev(prv){}
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
        unsigned num_V;
        unsigned num_E;
        unsigned corner_vertex_flag;

        vector<unsigned> head;
        vector<Edge> edges;

        void init(unsigned v) {
            num_V = v;
            corner_vertex_flag = v;
            head.resize(v, 0);
            edges.clear();
            num_E = 1;
            edges.emplace_back(-1, -1, -1.0, -1, -1); // stop edge.
        }

        void addEdge(unsigned u, unsigned v, double w) {
            if (edges.size() == num_E){
                edges.emplace_back(u, v, w, head[u], 0);
            }
            else{
                edges[num_E] = Edge(u, v, w, head[u], 0);
            }
            edges[head[u]].prev = num_E;
            head[u] = num_E;
            num_E++;
        }

        void removeEdge(unsigned eid){
            int prev = edges[eid].prev, next = edges[eid].next;
            if (!prev){
                head[edges[eid].from] = next;
            }
            else{
                edges[prev].next = next;
            }
            if (next > 0){
                edges[next].prev = prev;
            }
            num_E--;
        }

        int addVertex(){
            if (head.size() == num_V){
                head.emplace_back(0);
            }
            else{
                head[num_V] = 0;
            }
            return num_V++;
        }

        void removeVertex(unsigned vid){
            head[vid] = 0;
            num_V--;
        }

        bool isCorner(unsigned vid) {
            return vid < corner_vertex_flag;
        }
    };

    Graph my_base_graph;

    // s,t dijkstra query
    pair<double, double> dijkstra(const Graph &g, const unsigned &s, const unsigned &t){
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
            for (auto eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                auto v = g.edges[eid].to;
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

    int collect_covered_vertices(Graph &g, int s, double bound, vector<int> &covered_vertices, int &boundary_vertex_flag){
        vector<bool> vis(g.num_V, false);
        vector<double> d;
        d.resize(g.num_V, Base::unreachable);
        covered_vertices.clear();
        vector<int> fa(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        while (!q.empty()) {
            QNode f = q.top();
            q.pop();
            if (Base::doubleCmp(d[f.p] - bound) > 0) break;
            if (f.p < boundary_vertex_flag) covered_vertices.emplace_back(f.p);
            if (vis[f.p]) continue;
            for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next) {
                int v = g.edges[eid].to;
                double w = g.edges[eid].w;
                if (Base::doubleCmp(d[f.p] + w - d[v]) < 0) {
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
        }
        return covered_vertices.size();
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

    double queryGraphA2A(Graph &g, const Base::Point &s, const unsigned &fid_s, const Base::Point &t, const unsigned &fid_t,
                         const map<unsigned, vector<unsigned> > &face_point_map, const map<unsigned, Base::Point> &point_location_map){
        unsigned V_flag = g.num_V, E_flag = g.num_E;
        //  add edges from s to its neighbor Steiner points
        unsigned sid = g.addVertex();
        for (auto pid: face_point_map[fid_s]){
            auto pt = point_location_map[pid];
            double dis = sqrt(CGAL::squared_distance(s, pt));
            g.addEdge(sid, pid, dis);
        }
        //  add edges from t's neighbor Steiner points to t
        unsigned tid = g.addVertex();
        for (auto pid: face_point_map[fid_t]){
            auto pt = point_location_map[pid];
            double dis = sqrt(CGAL::squared_distance(t, pt));
            g.addEdge(pid, tid, dis);
        }
        double res = dijkstra(g, sid, tid).first;

        while (g.num_E > E_flag){
            int eid = g.num_E - 1;
            g.removeEdge(eid);
        }
        while (g.num_V > V_flag){
            int vid = g.num_V - 1;
            g.removeVertex(vid);
        }
        return res;
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
            auto hed = CGAL::SM_Halfedge_index(halfedge_id);
            adj_face_ids.push_back(m.face(hed).idx());
            adj_face_ids.push_back(m.face(m.opposite(hed)).idx());
        }

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
            if (t == f.p){
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
//                    cout << "theta_m = " << theta_m << endl;
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
                            double t_delta_I = max(delta_I, e_len / 12);
                            while (j <= floor(e_len / t_delta_I)){
                                int v_c_id = -1;
                                decltype(v_a) v_c;
                                if (edge_cut_vertex.find(L_id[x]) == edge_cut_vertex.end() || edge_cut_vertex[L_id[x]].size() < j){
                                    v_c = v_a + (v_b - v_a) / e_len * j * t_delta_I; // cut-vertex coordinate
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
