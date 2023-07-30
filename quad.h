#ifndef WEIGHTED_DISTANCE_ORACLE_QUAD_H
#define WEIGHTED_DISTANCE_ORACLE_QUAD_H

#include "base.h"
#include "k_skip.h"
#include "weighted_distance_oracle.h"

namespace Quad{

    using namespace std;

    int WSPD_hit = 0;
    kSkip::Graph my_base_graph;
    map<pair<int, int>, float> EAR_distance_map;
    float query_construction_time = 0.0;
    float dijkstra_running_time = 0.0;

    struct my_point{
        float x, y;
        my_point(float a, float b):x(a), y(b){}
        my_point(){}
        my_point operator - (const my_point& b) const{
            return my_point(x - b.x, y - b.y);
        }
    };

    int cross_product(my_point p1, my_point p2){
        return Base::floatCmp(p1.x * p2.y - p2.x * p1.y);
    }

    bool segment_intersection(my_point a, my_point b, my_point c, my_point d){
        if (Base::floatCmp(std::max(a.x, b.x) - std::min(c.x, d.x)) < 0 || Base::floatCmp(std::max(a.y, b.y) - std::min(c.y, d.y)) < 0 ||
            Base::floatCmp(std::max(c.x, d.x) - std::min(a.x, b.x)) < 0 || Base::floatCmp(std::max(c.y, d.y) - std::min(a.y, b.y)) < 0)
            return false;
        float dir1 = (c.x - a.x) * (b.y - a.y) - (b.x - a.x) * (c.y - a.y);
        float dir2 = (d.x - a.x) * (b.y - a.y) - (b.x - a.x) * (d.y - a.y);
        float dir3 = (a.x - c.x) * (d.y - c.y) - (d.x - c.x) * (a.y - c.y);
        float dir4 = (b.x - c.x) * (d.y - c.y) - (d.x - c.x) * (b.y - c.y);
        return Base::floatCmp(dir1 * dir2) <= 0 && Base::floatCmp(dir3 * dir4) <= 0;
    }

    //  triangle:ABC, segment:PQ
    bool triangle_segment_intersection(my_point a, my_point b, my_point c, my_point p, my_point q){
        int dir1 = cross_product(b - a, p - a);
        int dir2 = cross_product(c - b, p - b);
        int dir3 = cross_product(a - c, p - c);
        bool flag_p = false, flag_q = false;
        if (dir1 == dir2 && dir2 == dir3) flag_p = true;
        dir1 = cross_product(b - a, q - a);
        dir2 = cross_product(c - b, q - b);
        dir3 = cross_product(a - c, q - c);
        if (dir1 == dir2 && dir2 == dir3) flag_q = true;
        if (flag_p && flag_q) return true;
        return segment_intersection(a, b, p, q) || segment_intersection(b, c, p, q) || segment_intersection(c, a, p, q);
    }

    class treeNode{
    public:
        float x_min, y_min, x_max, y_max;
        int node_id;
        set<unsigned> boundary_points_id;
        set<unsigned> covered_faces_id;
        treeNode* parent;
        vector<treeNode*> sons; // 0/1/2/3 -- NE/NW/SW/SE
        treeNode();
        void print();
    };

    treeNode::treeNode() {
        node_id = 0;
        x_min = 1e60; x_max = -1e60;
        y_min = 1e60; y_max = -1e60;
        parent = nullptr;
        boundary_points_id.clear();
        covered_faces_id.clear();
        sons.clear();
    }

    void treeNode::print() {
        cout << "-----------------------------" << endl;
        cout << "(" << x_min << "," << y_min << ")-(" << x_max << "," << y_max << ")" << endl;
        cout << "# of sons " << sons.size() << endl;
        cout << "boundary points: " << endl;
        for (auto pid: boundary_points_id){
            cout << pid << " ";
        }
        cout << endl;
        cout << "-----------------------------" << endl;
    }

    class quadTree{
    public:
        quadTree(const Base::Mesh &mesh, map<unsigned, vector<unsigned> > &face_point_map);
        int level;
        treeNode* root;
        int node_count;
        vector<vector<treeNode* > > level_nodes; // level 0 is the root

        void buildLevel(const Base::Mesh &mesh, map<unsigned, vector<unsigned> > &face_point_map);

    };

    set<unsigned> extractIntersectFaces(const Base::Mesh &mesh, treeNode* tree_node, treeNode* fa_node = nullptr){
        set<unsigned> intersect_face_id = {};
        float x_min = tree_node->x_min, x_max = tree_node->x_max, y_min = tree_node->y_min, y_max = tree_node->y_max;
        if (fa_node == nullptr){
            for (auto &f: mesh.faces()){
                bool in_flag = false;
                vector<my_point> p;
                for (auto v: mesh.vertices_around_face(mesh.halfedge(f))){
                    float x = mesh.points()[v].x();
                    float y = mesh.points()[v].y();
                    p.emplace_back(x, y);
                }
                if (triangle_segment_intersection(p[0], p[1], p[2], my_point(x_min, y_min), my_point(x_max, y_min)) ||
                    triangle_segment_intersection(p[0], p[1], p[2], my_point(x_min, y_min), my_point(x_min, y_max)) ||
                    triangle_segment_intersection(p[0], p[1], p[2], my_point(x_min, y_max), my_point(x_max, y_max)) ||
                    triangle_segment_intersection(p[0], p[1], p[2], my_point(x_max, y_max), my_point(x_max, y_min))
                        ){
                    intersect_face_id.insert(f.idx());
                    tree_node->covered_faces_id.insert(f.idx());
                }
                if (!in_flag){
                    for (auto i = 0; i < 3; i++){
                        if (Base::floatCmp(p[i].x - x_min) >= 0 && Base::floatCmp(p[i].x - x_max) <= 0 &&
                            Base::floatCmp(p[i].y - y_min) >= 0 && Base::floatCmp(p[i].y - y_max) <= 0){
                            in_flag = true;
                            break;
                        }
                    }
                }
                if (in_flag){
                    tree_node->covered_faces_id.insert(f.idx());
                }
            }
        }
        else{
            for (auto &fid: fa_node->covered_faces_id){
                auto f = *(mesh.faces_begin() + fid);
                bool in_flag = false;
                vector<my_point> p;
                for (auto v: mesh.vertices_around_face(mesh.halfedge(f))){
                    float x = mesh.points()[v].x();
                    float y = mesh.points()[v].y();
                    p.emplace_back(x, y);
                }
                if (triangle_segment_intersection(p[0], p[1], p[2], my_point(x_min, y_min), my_point(x_max, y_min)) ||
                    triangle_segment_intersection(p[0], p[1], p[2], my_point(x_min, y_min), my_point(x_min, y_max)) ||
                    triangle_segment_intersection(p[0], p[1], p[2], my_point(x_min, y_max), my_point(x_max, y_max)) ||
                    triangle_segment_intersection(p[0], p[1], p[2], my_point(x_max, y_max), my_point(x_max, y_min))
                        ){
                    intersect_face_id.insert(f.idx());
                    tree_node->covered_faces_id.insert(f.idx());
                }
                if (!in_flag){
                    for (auto i = 0; i < 3; i++){
                        if (Base::floatCmp(p[i].x - x_min) >= 0 && Base::floatCmp(p[i].x - x_max) <= 0 &&
                            Base::floatCmp(p[i].y - y_min) >= 0 && Base::floatCmp(p[i].y - y_max) <= 0){
                            in_flag = true;
                            break;
                        }
                    }
                }
                if (in_flag){
                    tree_node->covered_faces_id.insert(f.idx());
                }
            }
        }
        return intersect_face_id;
    }

    quadTree::quadTree(const Base::Mesh &mesh, map<unsigned, vector<unsigned> > &face_point_map) {
        level = 0;
        root = new treeNode();
        node_count = 0;
        root->node_id = node_count++;
        for (auto &v: mesh.vertices()){
//            cout << mesh.points()[v].x() << " " << mesh.points()[v].y() << " " << mesh.points()[v].z() << endl;
            float x = mesh.points()[v].x(), y = mesh.points()[v].y();
            if (Base::floatCmp(x - root->x_min) < 0) root->x_min = x;
            if (Base::floatCmp(x - root->x_max) > 0) root->x_max = x;
            if (Base::floatCmp(y - root->y_min) < 0) root->y_min = y;
            if (Base::floatCmp(y - root->y_max) > 0) root->y_max = y;
        }
        set<unsigned> intersect_faces = extractIntersectFaces(mesh, root);
        for (auto fid: intersect_faces){
            auto fd = *(mesh.faces_begin() + fid);
            for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd))){
//                cout << "vid = " << vd.idx() << endl;
                root->boundary_points_id.insert(vd.idx());
            }
//            for (auto pid: face_point_map[fid]){
//                root->boundary_points_id.insert(pid);
//            }
        }
//        cout << "Root maximum side length: " << max(root->x_max - root->x_min, root->y_max - root->y_min) << endl;
//        for (auto &val: root.boundary_points_id){
//            cout << val << " ";
//        }
//        cout << endl;
//        for (auto &val: intersect_faces){
//            cout << val << " ";
//        }
//        cout << endl;
        vector<treeNode*> cur_level_nodes;
        cur_level_nodes.push_back(root);
        level_nodes.push_back(cur_level_nodes);
    }

    void quadTree::buildLevel(const Base::Mesh &mesh, map<unsigned, vector<unsigned> > &face_point_map) {
        vector<treeNode*> last_level_nodes = level_nodes[level];
        vector<treeNode*> cur_level_nodes;
        level++;
        for (auto &node: last_level_nodes){
            float pivot_x = 0.5 * (node->x_min + node->x_max),
                   pivot_y = 0.5 * (node->y_min + node->y_max);
            // one node will generate four sons:
            // NW(x_min, pivot_y)-(pivot_x, y_max)
            treeNode* NW_son = new treeNode();
            NW_son->x_min = node->x_min;
            NW_son->y_min = pivot_y;
            NW_son->x_max = pivot_x;
            NW_son->y_max = node->y_max;
            // NE(pivot_x, pivot_y)-(x_max, y_max)
            treeNode* NE_son = new treeNode();
            NE_son->x_min = pivot_x;
            NE_son->y_min = pivot_y;
            NE_son->x_max = node->x_max;
            NE_son->y_max = node->y_max;
            // SW(x_min, y_min)-(pivot_x, pivot_y)
            treeNode* SW_son = new treeNode();
            SW_son->x_min = node->x_min;
            SW_son->y_min = node->y_min;
            SW_son->x_max = pivot_x;
            SW_son->y_max = pivot_y;
            treeNode* SE_son = new treeNode();
            SE_son->x_min = pivot_x;
            SE_son->y_min = node->y_min;
            SE_son->x_max = node->x_max;
            SE_son->y_max = pivot_y;

            set<unsigned> faces = extractIntersectFaces(mesh, NW_son, node);
            for (auto fid: faces){
                auto fd = *(mesh.faces_begin() + fid);
                for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd))){
                    NW_son->boundary_points_id.insert(vd.idx());
                }
            }

            faces = extractIntersectFaces(mesh, NE_son, node);
            for (auto fid: faces){
                auto fd = *(mesh.faces_begin() + fid);
                for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd))){
                    NE_son->boundary_points_id.insert(vd.idx());
                }
            }

            faces = extractIntersectFaces(mesh, SW_son, node);
            for (auto fid: faces){
                auto fd = *(mesh.faces_begin() + fid);
                for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd))){
                    SW_son->boundary_points_id.insert(vd.idx());
                }
            }

            faces = extractIntersectFaces(mesh, SE_son, node);
            for (auto fid: faces){
                auto fd = *(mesh.faces_begin() + fid);
                for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd))){
                    SE_son->boundary_points_id.insert(vd.idx());
                }
            }

            NW_son->parent = node;
            NE_son->parent = node;
            SW_son->parent = node;
            SE_son->parent = node;
            NW_son->node_id = node_count++;
            NE_son->node_id = node_count++;
            SW_son->node_id = node_count++;
            SE_son->node_id = node_count++;
            node->sons.push_back(NW_son);
            node->sons.push_back(NE_son);
            node->sons.push_back(SW_son);
            node->sons.push_back(SE_son);
            cur_level_nodes.push_back(NW_son);
            cur_level_nodes.push_back(NE_son);
            cur_level_nodes.push_back(SW_son);
            cur_level_nodes.push_back(SE_son);
            node_count += 4;
        }
        level_nodes.push_back(cur_level_nodes);
    }

    kSkip::Graph generateSpanner(const vector<unsigned> &pid_list, set<WeightedDistanceOracle::nodePair> &node_pairs,
                                          map<unsigned, unsigned> &new_id){
        new_id.clear();
        for (auto i = 0; i < pid_list.size(); i++){
            new_id[pid_list[i]] = i;
        }
        kSkip::Graph g;
        g.init(static_cast<int>(pid_list.size()));
        for (auto node_pair: node_pairs){
            auto u = node_pair.first->center_idx;
            auto v = node_pair.second->center_idx;
            assert(new_id.find(u) != new_id.end());
            assert(new_id.find(v) != new_id.end());
            float tmp_d = WeightedDistanceOracle::enhanced_edges[make_pair(u, v)];
            g.addEdge(new_id[u], new_id[v], tmp_d);
            g.addEdge(new_id[v], new_id[u], tmp_d);
        }
        return g;
    }

    void distancePreprocessing(quadTree &quad_tree, kSkip::Graph &base_graph, map<unsigned, vector<unsigned> > &face_point_map){
        auto leaf = quad_tree.level_nodes[quad_tree.level];
        for (auto node: leaf){
            set<int> covered_point = {};
            for (auto fid: node->covered_faces_id){
                for (auto pid: face_point_map[fid]){
                    covered_point.insert(pid);
                }
            }
//            cout << "face = " << node->covered_faces_id.size() << endl;
//            cout << "size = " << covered_point.size() << endl;
            for (auto s: node->boundary_points_id){
                vector<float> d;
                kSkip::covered_dijkstra(base_graph, s, covered_point, d);
                for (auto pid: covered_point){
                    EAR_distance_map[make_pair(s, pid)] = d[pid];
                    EAR_distance_map[make_pair(pid, s)] = d[pid];
                }
            }

        }
    }

    // new version of distance query processing
    pair<float, bool> queryA2A(kSkip::Graph &spanner,
                                kSkip::Graph &base_graph, map<unsigned, vector<unsigned> > &face_point_map,
                                WeightedDistanceOracle::PartitionTree &tree, set<WeightedDistanceOracle::nodePair> &node_pairs,
                                map<unsigned, Base::Point> &point_location_map,
                                const Base::Point &s, const unsigned &fid_s,
                                const Base::Point &t, const unsigned &fid_t,
                                quadTree &quad_tree, map<unsigned, unsigned> &new_id, const unsigned &kappa){

        auto box_s = quad_tree.root, box_t = quad_tree.root;
        //find the leaf contains s and t
//        cout << fixed << setprecision(6) << "source: " << "x = " << s.x() << ", y=" << s.y() << endl;
//        cout << fixed << setprecision(6) << "root: " << "xmin = " << box_s->x_min << ", xmax = " << box_s->x_max << ", ymin = " << box_s->y_min << ", ymax = " << box_s->y_max << endl;
        while (box_s->sons.size() > 0){
            for (auto son: box_s->sons){
                if (Base::floatCmp(s.x() - son->x_min) >= 0 && Base::floatCmp(s.x() - son->x_max) <= 0 &&
                    Base::floatCmp(s.y() - son->y_min) >= 0 && Base::floatCmp(s.y() - son->y_max) <= 0){
                    box_s = son;
                    break;
                }
            }
        }
//        cout << "found box s" << endl;
//        cout << fixed << setprecision(6) << "target: " << "x = " << t.x() << ", y=" << t.y() << endl;
//        cout << fixed << setprecision(6) << "root: " << "xmin = " << box_t->x_min << ", xmax = " << box_t->x_max << ", ymin = " << box_t->y_min << ", ymax = " << box_t->y_max << endl;
//        for (auto son: box_t->sons) {
//            cout << fixed << setprecision(6) << "son.xmin = " << son->x_min << ", son.xmax = " << son->x_max << ", son.ymin = " << son->y_min << ", son.ymax = " << son->y_max << endl;
//        }
        while (box_t->sons.size() > 0){
            for (auto son: box_t->sons){
                if (Base::floatCmp(t.x() - son->x_min) >= 0 && Base::floatCmp(t.x() - son->x_max) <= 0 &&
                    Base::floatCmp(t.y() - son->y_min) >= 0 && Base::floatCmp(t.y() - son->y_max) <= 0){
                    box_t = son;
                    break;
                }
            }
        }
//        cout << "found box t" << endl;

        auto basegraph_vflag = base_graph.num_V, basegraph_eflag = base_graph.num_E;
        auto sid = base_graph.addVertex();
        for (auto pid: face_point_map[fid_s]){
            float dis = CGAL::squared_distance(s, point_location_map[pid]);
            base_graph.addEdge(sid, pid, sqrt(dis));
        }
        auto tid = base_graph.addVertex();
        for (auto pid: face_point_map[fid_t]){
            float dis = CGAL::squared_distance(t, point_location_map[pid]);
            base_graph.addEdge(pid, tid, sqrt(dis));
        }

        float res = -1.0;
//        auto basegraph_vflag = base_graph.num_V, basegraph_eflag = base_graph.num_E;

        if (box_s->node_id != box_t->node_id){

            auto spanner_vflag = spanner.num_V, spanner_eflag = spanner.num_E; //  backup
            auto final_s = spanner.addVertex();
            auto final_t = spanner.addVertex();

//             Gs and Gt local search for error bound guarantee on very small terrain surfaces
//            vector<float> d_Gs;
//            auto vertices_Gs = kSkip::hop_dijkstra(base_graph, sid, kappa, d_Gs);
//            for (auto &sp_id: vertices_Gs){
//                auto cur_id = spanner.addVertex();
//                spanner.addEdge(final_s, cur_id, d_Gs[sp_id]);
//            }
//
//            vector<float> d_Gt;
//            auto vertices_Gt = kSkip::hop_dijkstra(base_graph, tid, kappa, d_Gt);
//            for (auto &sp_id: vertices_Gt){
//                auto cur_id = spanner.addVertex();
//                spanner.addEdge(cur_id, final_t, d_Gt[sp_id]);
//            }
//
//            for (auto vid: vertices_Gs){
//                if (Base::floatCmp(d_Gt[vid] - Base::unreachable) < 0){
//                    float t_dis = d_Gs[vid] + d_Gt[vid];
//                    if (Base::floatCmp(res) < 0 || Base::floatCmp(t_dis - res) < 0){
//                        res = t_dis;
//                    }
//                }
//            }
//
//            d_Gs.clear();
//            vertices_Gs.clear();
//            d_Gt.clear();
//            vertices_Gt.clear();

            // result not found, find a highway-network path

            if (Base::floatCmp(res) < 0) {

                for (auto bpid: box_s->boundary_points_id) {
                    float d = Base::unreachable;
                    for (auto pid: face_point_map[fid_s]) {
                        //                    float t_dis = sqrt(CGAL::squared_distance(s, point_location_map[pid])) + EAR_distance_map[make_pair(closest_pid, bpid)];
                        float t_dis = sqrt(CGAL::squared_distance(s, point_location_map[pid])) +
                                      EAR_distance_map[make_pair(pid, bpid)];
                        if (Base::floatCmp(t_dis - d) < 0) d = t_dis;
                    }
                    spanner.addEdge(final_s, new_id[bpid], d);
                }

                for (auto bpid: box_t->boundary_points_id) {
                    float d = Base::unreachable;
                    for (auto pid: face_point_map[fid_t]) {
                        float t_dis = sqrt(CGAL::squared_distance(t, point_location_map[pid])) +
                                      EAR_distance_map[make_pair(pid, bpid)];
                        if (Base::floatCmp(t_dis - d) < 0) d = t_dis;
                    }
                    spanner.addEdge(new_id[bpid], final_t, d);
                }

                auto q_start = chrono::_V2::system_clock::now();    //  break down time evaluation

                res = kSkip::dijkstra(spanner, final_s, final_t).first;

                auto q_end = chrono::_V2::system_clock::now();
                auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);
                dijkstra_running_time += static_cast<float>(q_duration.count());
            }

            // recover spanner
            while (spanner.num_E > spanner_eflag){
                auto eid = spanner.num_E - 1;
                spanner.removeEdge(eid);
            }
            while (spanner.num_V > spanner_vflag){
                auto vid = spanner.num_V - 1;
                spanner.removeVertex(vid);
            }
//            return make_pair(res, box_s->node_id == box_t->node_id);
        }
        else {

            auto sid = base_graph.addVertex();
            for (auto pid: face_point_map[fid_s]){
                float dis = CGAL::squared_distance(s, point_location_map[pid]);
                base_graph.addEdge(sid, pid, sqrt(dis));
            }
            auto tid = base_graph.addVertex();
            for (auto pid: face_point_map[fid_t]){
                float dis = CGAL::squared_distance(t, point_location_map[pid]);
                base_graph.addEdge(pid, tid, sqrt(dis));
            }

            auto q_start = chrono::_V2::system_clock::now();    //  break down time evaluation

            res = kSkip::dijkstra(base_graph, sid, tid).first;

            auto q_end = chrono::_V2::system_clock::now();
            auto q_duration = chrono::duration_cast<chrono::microseconds>(q_end - q_start);
            dijkstra_running_time += static_cast<float>(q_duration.count());
        }

        // recover basegraph
        while (base_graph.num_E > basegraph_eflag) {
            auto eid = base_graph.num_E - 1;
            base_graph.removeEdge(eid);
        }
        while (base_graph.num_V > basegraph_vflag) {
            auto vid = base_graph.num_V - 1;
            base_graph.removeVertex(vid);
        }

        return make_pair(res, box_s->node_id == box_t->node_id);
    }

    // old version
    pair<float, bool> queryA2A(Base::Mesh &m, kSkip::Graph &spanner,
                                kSkip::Graph &base_graph, map<int, vector<int> > &face_point_map,
                                map<int, Base::Point> &point_location_map,
                                Base::Point s, int fid_s,
                                Base::Point t, int fid_t,
                                quadTree &quad_tree, map<int, int> &new_id){

        auto box_s = quad_tree.root, box_t = quad_tree.root;
        //find the leaf contains s and t

        while (box_s->sons.size() > 0){
            for (auto son: box_s->sons){
                if (Base::floatCmp(s.x() - son->x_min) >= 0 && Base::floatCmp(s.x() - son->x_max) <= 0 &&
                    Base::floatCmp(s.y() - son->y_min) >= 0 && Base::floatCmp(s.y() - son->y_max) <= 0){
                    box_s = son;
                    break;
                }
            }
        }
        while (box_t->sons.size() > 0){
            for (auto son: box_t->sons){
                if (Base::floatCmp(t.x() - son->x_min) >= 0 && Base::floatCmp(t.x() - son->x_max) <= 0 &&
                    Base::floatCmp(t.y() - son->y_min) >= 0 && Base::floatCmp(t.y() - son->y_max) <= 0){
                    box_t = son;
                    break;
                }
            }
        }


        if (box_s->node_id != box_t->node_id){

            int V_flag = spanner.num_V, E_flag = spanner.num_E; //  backup

            auto final_s = spanner.addVertex();

            for (auto bpid: box_s->boundary_points_id){
                float d = Base::unreachable;
                for (auto pid: face_point_map[fid_s]){
                    float t_dis = sqrt(CGAL::squared_distance(s, point_location_map[pid])) + EAR_distance_map[make_pair(pid, bpid)];
                    if (Base::floatCmp(t_dis - d) < 0) d = t_dis;
                }
                spanner.addEdge(final_s, new_id[bpid], d);
            }
            auto final_t = spanner.addVertex();

            for (auto bpid: box_t->boundary_points_id){
                float d = Base::unreachable;
                for (auto pid: face_point_map[fid_t]){
                    float t_dis = sqrt(CGAL::squared_distance(t, point_location_map[pid])) + EAR_distance_map[make_pair(pid, bpid)];
                    if (Base::floatCmp(t_dis - d) < 0) d = t_dis;
                }
                spanner.addEdge(new_id[bpid], final_t, d);
            }
            float res = kSkip::dijkstra(spanner, final_s, final_t).first;

            while (spanner.num_E > E_flag){
                int eid = spanner.num_E - 1;
                spanner.removeEdge(eid);
            }
            while (spanner.num_V > V_flag){
                int vid = spanner.num_V - 1;
                spanner.removeVertex(vid);
            }

            return make_pair(res, box_s->node_id == box_t->node_id);
        }
        else{
            int V_flag = base_graph.num_V, E_flag = base_graph.num_E;
            auto sid = base_graph.addVertex();
            for (auto pid: face_point_map[fid_s]){
                float dis = CGAL::squared_distance(s, point_location_map[pid]);
                base_graph.addEdge(sid, pid, sqrt(dis));
            }
            auto tid = base_graph.addVertex();
            for (auto pid: face_point_map[fid_t]){
                float dis = CGAL::squared_distance(t, point_location_map[pid]);
                base_graph.addEdge(pid, tid, sqrt(dis));
            }
            float res = kSkip::dijkstra(base_graph, sid, tid).first;

            while (base_graph.num_E > E_flag){
                int eid = base_graph.num_E - 1;
                base_graph.removeEdge(eid);
            }
            while (base_graph.num_V > V_flag){
                int vid = base_graph.num_V - 1;
                base_graph.removeVertex(vid);
            }

            return make_pair(res, box_s->node_id == box_t->node_id);
        }
    }

}

#endif //WEIGHTED_DISTANCE_ORACLE_QUAD_H
