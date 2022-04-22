#ifndef WEIGHTED_DISTANCE_ORACLE_K_SKIP_H
#define WEIGHTED_DISTANCE_ORACLE_K_SKIP_H
#include "base.h"

// parameter for K-Algo
extern int K;

namespace kSkip{
    using namespace std;

    struct Edge{
        float w;
        int from, to, next, prev;
        Edge(int u, int v, float len, int nxt, int prv): from(u), to(v),  w(len), next(nxt), prev(prv){}
    };

    struct QNode{
        int p;
        int hop_cnt;
        float dis;
        bool operator < (const QNode& x) const{
            return Base::floatCmp(dis - x.dis) > 0 || !Base::floatCmp(dis - x.dis) && p < x.p;
        }
        QNode(int pid, float d_val, int hops = 0): p(pid), dis(d_val), hop_cnt(hops){}
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

        void addEdge(unsigned u, unsigned v, float w) {
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
    pair<float, float> dijkstra(const Graph &g, const unsigned &s, const unsigned &t){
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        vector<bool> vis(g.num_V, false);
        vector<float> d(g.num_V, Base::unreachable);
        vector<int> fa(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        while (!q.empty()){
            QNode f = q.top(); q.pop();
            if (f.p == t) {
                break;
            }
            if (vis[f.p]) continue;
            for (auto eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                auto v = g.edges[eid].to;
                float w = g.edges[eid].w;
                if (Base::floatCmp(d[f.p] + w - d[v]) < 0){
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
        }
        assert(Base::floatCmp(d[t] - Base::unreachable) < 0);

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return make_pair(d[t], static_cast<float>(duration.count()));
    }

    //dijkstra for a given hop, return the discovered vertex id
    vector<int> hop_dijkstra(Graph &g, int s, const unsigned &kappa, vector<float> &d){
        vector<bool> vis(g.num_V, false);
        d.resize(g.num_V, Base::unreachable);
        vector<int> fa(g.num_V, -1);
        vector<int> hop(g.num_V, -1);
        d[s] = 0;
        hop[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        vector<int> discovered_point = {};
        while (!q.empty()){
            QNode f = q.top(); q.pop();
            if (vis[f.p]) continue;
            for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                float w = g.edges[eid].w;
                if (hop[f.p] + 1 < kappa && Base::floatCmp(d[f.p] + w - d[v]) < 0){
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    hop[v] = hop[f.p] + 1;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
            discovered_point.emplace_back(f.p);
        }

        return discovered_point;
    }

    //dijkstra for a given distance bound
    float bounded_dijkstra(Graph &g, int s, float bound, vector<float> &d){
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        vector<bool> vis(g.num_V, false);
        d.resize(g.num_V, Base::unreachable);
        vector<int> fa(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        while (!q.empty()){
            QNode f = q.top(); q.pop();
            if (Base::floatCmp(d[f.p] - bound) > 0) break;
            if (vis[f.p]) continue;
            for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                float w = g.edges[eid].w;
                if (Base::floatCmp(d[f.p] + w - d[v]) < 0){
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<float>(duration.count());
    }

    // find covered non-corner vertex, return their id
    int collect_covered_vertices(Graph &g, int s, float bound, vector<int> &covered_vertices, int &boundary_vertex_flag){
        vector<bool> vis(g.num_V, false);
        vector<float> d;
        d.resize(g.num_V, Base::unreachable);
        covered_vertices.clear();
        vector<int> fa(g.num_V, -1);
        d[s] = 0;
        priority_queue<QNode> q = {};
        q.push(QNode(s, d[s]));
        while (!q.empty()) {
            QNode f = q.top();
            q.pop();
            if (Base::floatCmp(d[f.p] - bound) > 0) break;
            if (f.p < boundary_vertex_flag) covered_vertices.emplace_back(f.p);
            if (vis[f.p]) continue;
            for (int eid = g.head[f.p]; eid; eid = g.edges[eid].next) {
                int v = g.edges[eid].to;
                float w = g.edges[eid].w;
                if (Base::floatCmp(d[f.p] + w - d[v]) < 0) {
                    d[v] = d[f.p] + w;
                    fa[v] = f.p;
                    q.push(QNode(v, d[v]));
                }
            }
            vis[f.p] = true;
        }
        return covered_vertices.size();
    }

    float covered_dijkstra(Graph &g, int s, set<int> &cover_id, vector<float> &d){
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
                float w = g.edges[eid].w;
                if (Base::floatCmp(d[f.p] + w - d[v]) < 0){
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
        return static_cast<float>(duration.count());
    }

    float queryGraphA2A(Graph &g, const Base::Point &s, const unsigned &fid_s, const Base::Point &t, const unsigned &fid_t,
                         map<unsigned, vector<unsigned> > &face_point_map, map<unsigned, Base::Point> &point_location_map){
        unsigned V_flag = g.num_V, E_flag = g.num_E;
        //  add edges from s to its neighbor Steiner points
        unsigned sid = g.addVertex();
        for (auto pid: face_point_map[fid_s]){
            auto pt = point_location_map[pid];
            float dis = sqrt(CGAL::squared_distance(s, pt));
            g.addEdge(sid, pid, dis);
        }
        //  add edges from t's neighbor Steiner points to t
        unsigned tid = g.addVertex();
        for (auto pid: face_point_map[fid_t]){
            auto pt = point_location_map[pid];
            float dis = sqrt(CGAL::squared_distance(t, pt));
            g.addEdge(pid, tid, dis);
        }
        float res = dijkstra(g, sid, tid).first;

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

    void constructMeshGraph(const Base::Mesh &mesh, Graph &g){
        g.init(mesh.num_vertices());
        g.head.resize(mesh.num_vertices(), 0);
        for (auto &fd: mesh.faces()){
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))){
                auto sid = mesh.source(hed);
                auto tid = mesh.target(hed);
                float w = sqrt(CGAL::squared_distance(mesh.points()[sid], mesh.points()[tid]));
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

    unsigned initialBisectors(const Base::Mesh &mesh,
                              map<unsigned, vector<unsigned> > &edge_bisector_map,
                              map<unsigned, set<unsigned> > &bisector_edge_map,
                              vector<pair<unsigned, vector<Base::Point> > > &bisector_info,
                              map<unsigned, unsigned> &bisector_face_map,
                              map<unsigned, vector<unsigned> > &face_bisector_map){
        unsigned bisector_id = 0;
        for (auto fd: mesh.faces()){
            auto fid = fd.idx();
            vector<Base::Point> p(3);
            vector<unsigned> pid(3);
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))){
                unsigned eid = mesh.edge(hed).idx(), eid2 = mesh.edge(mesh.prev(hed)).idx();
                for (auto i = 0; i != 3; hed = mesh.next(hed), i++){
                    p[i] = mesh.points()[mesh.source(hed)];
                    pid[i] = mesh.source(hed).idx();
                }
                bisector_info.emplace_back(pid[0], p);
                bisector_edge_map[bisector_id].insert(eid);
                bisector_edge_map[bisector_id].insert(eid2);
                bisector_face_map[bisector_id] = fid;
                face_bisector_map[fid].push_back(bisector_id);
                edge_bisector_map[eid].push_back(bisector_id);
                edge_bisector_map[eid2].push_back(bisector_id++);
            }
        }
        return bisector_id;
    }

    bool placeSteinerPointOnBisector(
            const unsigned &bisector_id,
            const unsigned &bisector_source_point,
            Graph &g,
            const float &gama,
            const float &eps,
            map<unsigned, vector<unsigned> > &bisector_point_map,
            map<unsigned, Base::Point> &point_location_map,
            map<unsigned, unsigned> &point_bisector_map,
            vector<Base::Point> &p){
        if (bisector_point_map[bisector_id].size() > 0) return 0;
//        unsigned vertex_id = g.addVertex();
        unsigned vertex_id_bk = g.num_V;
        vector<Base::Point> cur_bisector_p = {};
        vector<unsigned> cur_bisector_p_id = {};

//        point_bisector_map[vertex_id] = bisector_id;
        cur_bisector_p.emplace_back(p[0]);
        cur_bisector_p_id.emplace_back(vertex_id_bk++);
//        point_location_map[vertex_id] = p[0];
//        bisector_point_map[bisector_id].emplace_back(vertex_id);



        float len1 = sqrt(CGAL::squared_distance(p[1], p[0])),
                len2 = sqrt(CGAL::squared_distance(p[2], p[0]));
        Base::Point p_end(p[1] + Base::Vector(p[2] - p[1]) * len1 / (len1 + len2));
        Base::Vector vec_bisector(p_end - p[0]);
        Base::Point aux1 = p[0] + Base::Vector(p[1] - p[0]) * gama / len1;
        Base::Point aux2 = p[0] + Base::Vector(p[2] - p[0]) * gama / len2;
        Base::Point bisector_p0 = aux1 + 0.5 * Base::Vector(aux2 - aux1);
        float angle = Base::PI * CGAL::approximate_angle(p[1], p[0], p[2]) / 180;
        float sin_val = sin(angle * 0.5);
        float limit_distance = sqrt(CGAL::squared_distance(p[0], p_end));
        float cur_distance = sqrt(CGAL::squared_distance(p[0], bisector_p0));

        double num_Steiner_points = 1.61 / sin(angle) * log(2 * limit_distance / gama);
        num_Steiner_points *= 1 / sqrt(eps) * log(2 / eps);


        bool uniform_flag = 0;
        while (Base::floatCmp(cur_distance - limit_distance) < 0) {
            if (cur_bisector_p.size() > 25) {
                uniform_flag = 1;
                break;
            }
            Base::Point bisector_p = p[0] + cur_distance / limit_distance * vec_bisector;
//            vertex_id = g.addVertex();
            cur_bisector_p_id.emplace_back(vertex_id_bk++);
            cur_bisector_p.emplace_back(bisector_p);
            float distance_delta = sin_val * sqrt(0.5 * eps) * cur_distance;
            cur_distance += distance_delta;
        }

        if (uniform_flag){
            float distance_delta = limit_distance / 25;
            cur_distance = sqrt(CGAL::squared_distance(p[0], bisector_p0));
            vertex_id_bk = g.num_V;

            cur_bisector_p_id.clear();
            cur_bisector_p.clear();
            cur_bisector_p.emplace_back(p[0]);
            cur_bisector_p_id.emplace_back(vertex_id_bk++);

            while (Base::floatCmp(cur_distance - limit_distance) < 0) {
                Base::Point bisector_p = p[0] + cur_distance / limit_distance * vec_bisector;
//                vertex_id = g.addVertex();
                cur_bisector_p_id.emplace_back(vertex_id_bk++);
                cur_bisector_p.emplace_back(bisector_p);
                cur_distance += distance_delta;
            }
        }

        bisector_point_map[bisector_id] = cur_bisector_p_id;
        for (auto i = 0; i < cur_bisector_p_id.size(); i++){
            g.addVertex();
            auto bp_id = cur_bisector_p_id[i];
            point_location_map[bp_id] = cur_bisector_p[i];
            point_bisector_map[bp_id] = bisector_id;
        }
//        cout << "placed points = " << cur_bisector_p_id.size() << endl;

        return 1;
    }

    //Unfixed on-the-fly implementation
    float unfixedOnTheFly(
            Base::Mesh mesh, Graph g,
            const float &eps,
            vector<float> &gama,
            const vector<float> &face_weight,
            const Base::Point &point_s, const unsigned &fid_s,
            const Base::Point &point_t, const unsigned &fid_t) {

        map<unsigned, vector<unsigned> > bisector_point_map;
        map<unsigned, Base::Point> point_location_map;
        map<unsigned, set<unsigned> > bisector_edge_map;
        map<unsigned, unsigned> point_bisector_map;
        map<unsigned, unsigned> bisector_face_map;
        map<unsigned, vector<unsigned> > edge_bisector_map;
        map<unsigned, vector<unsigned> > face_bisector_map;
        vector<pair<unsigned, vector<Base::Point> > > bisector_info;
        map<pair<unsigned, unsigned>, bool> edge_map = {};

        initialBisectors(mesh, edge_bisector_map, bisector_edge_map, bisector_info, bisector_face_map, face_bisector_map);

        unsigned s = g.addVertex();

        auto fds = *(mesh.faces_begin() + fid_s);

        for (auto bisector_id: face_bisector_map[fds.idx()]){
            unsigned bisector_source_point = bisector_info[bisector_id].first;
            placeSteinerPointOnBisector(bisector_id, bisector_source_point, g, gama[bisector_source_point],
                                        eps, bisector_point_map, point_location_map, point_bisector_map, bisector_info[bisector_id].second);
            for (auto sp_id: bisector_point_map[bisector_id]){
                Base::Point sp_p = point_location_map[sp_id];
                float dis = sqrt(CGAL::squared_distance(point_s, sp_p));
                edge_map[make_pair(s, sp_id)] = 1;
                g.addEdge(s, sp_id, dis);
            }
        }

        unsigned t = g.addVertex();

        auto fdt = *(mesh.faces_begin() + fid_t);
        for (auto bisector_id: face_bisector_map[fdt.idx()]){
            unsigned bisector_source_point = bisector_info[bisector_id].first;
            placeSteinerPointOnBisector(bisector_id, bisector_source_point, g, gama[bisector_source_point],
                                        eps, bisector_point_map, point_location_map, point_bisector_map , bisector_info[bisector_id].second);
            for (auto sp_id: bisector_point_map[bisector_id]){
                Base::Point sp_p = point_location_map[sp_id];
                float dis = sqrt(CGAL::squared_distance(point_t, sp_p));
                edge_map[make_pair(t, sp_id)] = 1;
                g.addEdge(sp_id, t, dis);
            }
        }

        unsigned vertices_num = g.num_V;
        vector<float> D;
        D.resize(vertices_num, -1);
        D[s] = 0.0;
        vector<bool> vis;
        vis.resize(vertices_num, 0);
        vis[s] = 1;
        priority_queue<QNode> q = {};
//        q.push(QNode(s, D[s]));


        for (auto eid = g.head[s]; eid > 0; eid = g.edges[eid].next){
            int v = g.edges[eid].to;
            float w = g.edges[eid].w;
            if (Base::floatCmp(D[v]) < 0 || Base::floatCmp(D[v] - D[s] - w) > 0){
                D[v] = D[s] + w;
                q.push(QNode(v, D[v]));
            }
        }


        while (!q.empty()){
            auto f = q.top(); q.pop();
//            cout << fixed << setprecision(6) << "f.dis = " << f.dis << endl;
//            cout << "cur f: " << f.p << " target: " << t << endl;
            auto cur_p = point_location_map[f.p];
//            cout << "cur p = " << cur_p << endl;
            if (f.p == t){
//                int xx; cin >> xx;
                return f.dis;
            }
            if (vis[f.p]) continue;
            auto cur_bisector_id = point_bisector_map[f.p];
//            cout << "cur bisector id = " << cur_bisector_id << endl;

            for (auto eid: bisector_edge_map[cur_bisector_id]){
//                cout << "edge bisector size = " << edge_bisector_map[eid].size() << endl;
                for (auto bisector_id: edge_bisector_map[eid]){
                    if (bisector_point_map[bisector_id].size() == 0){
                        unsigned bisector_source_point = bisector_info[bisector_id].first;
                        placeSteinerPointOnBisector(bisector_id, bisector_source_point, g, gama[bisector_source_point],
                                                    eps, bisector_point_map, point_location_map, point_bisector_map , bisector_info[bisector_id].second);
                        for (auto x = 0; x < bisector_point_map[bisector_id].size(); x++){
                            D.emplace_back(-1);
                            vis.push_back(0);
                        }

                    }
                    for (auto sp_id: bisector_point_map[bisector_id]){
                        Base::Point sp_p = point_location_map[sp_id];
                        float dis = Base::distanceSnell(mesh, face_weight,
                                                        cur_p, bisector_face_map[point_bisector_map[f.p]],
                                                        sp_p, bisector_face_map[bisector_id],
                                                        eid);
//                         cout << "edge len = " << dis << endl;
                        auto cur_edge = make_pair(f.p, sp_id);
                        if (edge_map.find(cur_edge) == edge_map.end()){
                            g.addEdge(f.p, sp_id, dis);
                            edge_map[cur_edge] = 1;
                        }
                    }
                }
            }
            for (auto eid = g.head[f.p]; eid > 0; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                float w = g.edges[eid].w;
                if (Base::floatCmp(D[v]) < 0 || Base::floatCmp(D[v] - D[f.p] - w) > 0){
                    D[v] = D[f.p] + w;
                    q.push(QNode(v, D[v]));
                }
            }
            vis[f.p] = 1;
        }

    }

    // VLDB'2015 implementation
    float computeDistanceBound(
            Base::Mesh mesh, Graph g,
            const unsigned &K, const Base::Point &point_s, const unsigned &fid_s,
            const Base::Point &point_t, const unsigned &fid_t,
            float l_min
    ){
//        cout << "s = " << point_s << " t = " << point_t << endl;
        map<pair<unsigned, unsigned>, bool> edge_map = {};
        auto vds = mesh.add_vertex(point_s);
        int s = vds.idx();
        g.addVertex();
        auto fds = *(mesh.faces_begin() + fid_s);
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fds))){
            mesh.add_edge(vds, vd);
            float dis = sqrt(CGAL::squared_distance(point_s, mesh.points()[vd]));
            edge_map[make_pair(s, vd.idx())] = 1;
            g.addEdge(s, vd.idx(), dis);
        }
        auto vdt = mesh.add_vertex(point_t);
        int t = vdt.idx();
        g.addVertex();
        auto fdt = *(mesh.faces_begin() + fid_t);
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fdt))){
            mesh.add_edge(vdt, vd);
            float dis = sqrt(CGAL::squared_distance(point_t, mesh.points()[vd]));
            edge_map[make_pair(vd.idx(), t)] = 1;
            g.addEdge(vd.idx(), t, dis);
        }


        map<int, vector<int> > edge_cut_vertex = {};
        map<int, int> cut_vertex_halfedge = {};
        int vertices_num = mesh.num_vertices();
        vector<float> D;
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
//            cout << "p = " << f.p << " D = " << f.dis << endl;
            if (t == f.p){
                return D[f.p];
            }
            if (vis[f.p]) continue;
            // S = S U {u};
            float d_ut = sqrt(CGAL::squared_distance(mesh.points()[CGAL::SM_Vertex_index(f.p)], mesh.points()[CGAL::SM_Vertex_index(t)]));

            auto u = CGAL::SM_Vertex_index(f.p);
//            cout << "u = " << u << endl;
            vector<int> adj_face_ids;
            getAdjacentFaces(mesh, g, cut_vertex_halfedge, f.p, adj_face_ids);
            for (auto fid: adj_face_ids){
                auto fd = CGAL::SM_Face_index(fid);
                if (fd != Base::Mesh::null_face()){
                    vector<int> vids = {};
                    pair<int, int> opposite_edge = {};
                    int opposite_edge_id = -1;
                    // vector<pair<int, int> > adjacent_edges = {};
                    vector<int> L_id = {};
                    auto hed = mesh.halfedge(fd);
                    vector<pair<int, int> > L = {};
                    for (int i = 0; i != 3; i++){
                        auto t_u = mesh.points()[mesh.source(hed)];
//                        cout << "vertices " << i << ": " << t_u << endl;
                        auto t_v = mesh.points()[mesh.target(hed)];
                        int uid = mesh.source(hed).idx();
                        int vid = mesh.target(hed).idx();
                        Base::Segment seg(t_u, t_v);
                        vids.push_back(uid);
                        if (!seg.has_on(mesh.points()[u])){
                            L.push_back(make_pair(uid, vid));
                            L_id.push_back(hed);
                        }
                        hed = mesh.next(hed);
                    }
                    float theta_m = -1.0;
                    for (int i = 0; i != 3; i++){
                        auto t_u = CGAL::SM_Vertex_index(vids[i]),
                                t_v = CGAL::SM_Vertex_index(vids[(i + 1) % 3]),
                                t_w = CGAL::SM_Vertex_index(vids[(i + 2) % 3]);
                        float angle = CGAL::approximate_angle(mesh.points()[t_u], mesh.points()[t_v], mesh.points()[t_w]);
                        if (Base::floatCmp(theta_m) < 0 || Base::floatCmp(angle - theta_m) < 0){
                            theta_m = angle;
                        }
                    }
//                    cout << "theta_m = " << theta_m << endl;
                    // compute |O_iO_i+1|min and delta_I.
                    float OO_min = 0.5 * l_min * sqrt(2 * (1 - cos(theta_m / 180 * Base::PI)));
                    float delta_I = OO_min / K;
                    // cout << "O_min = " << OO_min << endl;
//                     cout << "delta_I = " << delta_I << endl;

                    for (auto x = 0; x != L.size(); x++){
                        auto e = L[x];
//                        cout << "e = " << e.first << " " << e.second << endl;
                        auto v_a = mesh.points()[CGAL::SM_Vertex_index(e.first)];
                        auto v_b = mesh.points()[CGAL::SM_Vertex_index(e.second)];
//                        cout << "v_a = " << v_a << " v_b = " << v_b << endl;

                        float d_et = min(
                                sqrt(CGAL::squared_distance(v_a, mesh.points()[CGAL::SM_Vertex_index(t)])),
                                sqrt(CGAL::squared_distance(v_b, mesh.points()[CGAL::SM_Vertex_index(t)]))
                        );
                        if (Base::floatCmp(d_ut - d_et) >= 0){
                            int j = 1;
                            float e_len = sqrt(CGAL::squared_distance(v_a, v_b));
//                             cout << "num cut vertices = " << floor(e_len / delta_I) << endl;
                            float t_delta_I = delta_I;
//                            cout << "placed points = " << e_len / t_delta_I << endl;
                            if (Base::floatCmp(t_delta_I) <= 0 || e_len / t_delta_I > 55){
                                t_delta_I = e_len / 55;
                            }
                            while (j <= floor(e_len / t_delta_I)){
                                int v_c_id = -1;
                                decltype(v_a) v_c;
                                if (edge_cut_vertex.find(L_id[x]) == edge_cut_vertex.end() || edge_cut_vertex[L_id[x]].size() < j){
                                    v_c = v_a + (v_b - v_a) / e_len * j * t_delta_I; // cut-vertex coordinate
                                    mesh.add_vertex(v_c);  //  add the vertex to mesh
                                    v_c_id = g.addVertex(); //   vertex id
                                    assert(v_c != mesh.points()[u]);
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
                                    v_c = mesh.points()[CGAL::SM_Vertex_index(v_c_id)];
                                }
                                float d_u_vc = sqrt(CGAL::squared_distance(mesh.points()[u], v_c));
                                if (Base::floatCmp(d_u_vc - delta_I) >= 0){
                                    auto cur_edge = make_pair(f.p, v_c_id);
                                    if (edge_map.find(cur_edge) == edge_map.end()){
                                        g.addEdge(f.p, v_c_id, d_u_vc);
                                    }
                                }
                                // D[v_c_id] = D[f.u] + d_u_vc;
                                // q.insert(Qnode(v_c_id, D[v_c_id]));
                                j++;
                            }
//                            cout << "placed points = " << j << endl;
                        }
                    }
                    if (!g.isCorner(f.p)){
                        for (auto id: vids){
                            auto v = CGAL::SM_Vertex_index(id);
                            float w = sqrt(CGAL::squared_distance(mesh.points()[u], mesh.points()[v]));
                            auto cur_edge = make_pair(f.p, v.idx());
                            if (edge_map.find(cur_edge) == edge_map.end()){
                                g.addEdge(f.p, v.idx(), w);
                            }
                        }
                    }
                }
            }
            for (auto eid = g.head[f.p]; eid > 0; eid = g.edges[eid].next){
                int v = g.edges[eid].to;
                float w = g.edges[eid].w;
                if (Base::floatCmp(D[v]) < 0 || Base::floatCmp(D[v] - D[f.p] - w) > 0){
                    D[v] = D[f.p] + w;
                    q.push(QNode(v, D[v]));
                }
            }
            vis[f.p] = 1;
        }

    }
}

#endif //WEIGHTED_DISTANCE_ORACLE_K_SKIP_H
