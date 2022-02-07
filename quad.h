//
// Created by huang on 2021/8/23.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_QUAD_H
#define WEIGHTED_DISTANCE_ORACLE_QUAD_H

#include "base.h"
#include "k_skip.h"
#include "weighted_distance_oracle.h"

namespace Quad{

    using namespace std;

    int WSPD_hit = 0;
    kSkip::Graph my_base_graph;
    map<pair<int, int>, double> LQT_distance_map;

    struct my_point{
        double x, y;
        my_point(double a, double b):x(a), y(b){}
        my_point(){}
        my_point operator - (const my_point& b) const{
            return my_point(x - b.x, y - b.y);
        }
    };

    int cross_product(my_point p1, my_point p2){
        return Base::doubleCmp(p1.x * p2.y - p2.x * p1.y);
    }

    bool segment_intersection(my_point a, my_point b, my_point c, my_point d){
        if (Base::doubleCmp(std::max(a.x, b.x) - std::min(c.x, d.x)) < 0 || Base::doubleCmp(std::max(a.y, b.y) - std::min(c.y, d.y)) < 0 ||
            Base::doubleCmp(std::max(c.x, d.x) - std::min(a.x, b.x)) < 0 || Base::doubleCmp(std::max(c.y, d.y) - std::min(a.y, b.y)) < 0)
            return false;
        double dir1 = (c.x - a.x) * (b.y - a.y) - (b.x - a.x) * (c.y - a.y);
        double dir2 = (d.x - a.x) * (b.y - a.y) - (b.x - a.x) * (d.y - a.y);
        double dir3 = (a.x - c.x) * (d.y - c.y) - (d.x - c.x) * (a.y - c.y);
        double dir4 = (b.x - c.x) * (d.y - c.y) - (d.x - c.x) * (b.y - c.y);
        return Base::doubleCmp(dir1 * dir2) <= 0 && Base::doubleCmp(dir3 * dir4) <= 0;
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
        double x_min, y_min, x_max, y_max;
        int node_id;
        set<int> boundary_points_id;
        set<int> covered_faces_id;
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
        quadTree(Base::Mesh &m, map<int, vector<int> > &face_point_map);
        int level;
        treeNode* root;
        int node_count;
        vector<vector<treeNode* > > level_nodes; // level 0 is the root

        void buildLevel(Base::Mesh &m, map<int, vector<int> > &face_point_map);

    };

    set<int> extractIntersectFaces(Base::Mesh &m, treeNode* tree_node, treeNode* fa_node = nullptr){
        set<int> intersect_face_id = {};
        double x_min = tree_node->x_min, x_max = tree_node->x_max, y_min = tree_node->y_min, y_max = tree_node->y_max;
        if (fa_node == nullptr){
            for (auto &f: m.faces()){
                bool in_flag = false;
                vector<my_point> p;
                for (auto v: m.vertices_around_face(m.halfedge(f))){
                    double x = m.points()[v].x();
                    double y = m.points()[v].y();
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
                        if (Base::doubleCmp(p[i].x - x_min) >= 0 && Base::doubleCmp(p[i].x - x_max) <= 0 &&
                            Base::doubleCmp(p[i].y - y_min) >= 0 && Base::doubleCmp(p[i].y - y_max) <= 0){
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
                auto f = *(m.faces_begin() + fid);
                bool in_flag = false;
                vector<my_point> p;
                for (auto v: m.vertices_around_face(m.halfedge(f))){
                    double x = m.points()[v].x();
                    double y = m.points()[v].y();
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
                        if (Base::doubleCmp(p[i].x - x_min) >= 0 && Base::doubleCmp(p[i].x - x_max) <= 0 &&
                            Base::doubleCmp(p[i].y - y_min) >= 0 && Base::doubleCmp(p[i].y - y_max) <= 0){
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

    quadTree::quadTree(Base::Mesh &m, map<int, vector<int> > &face_point_map) {
        level = 0;
        root = new treeNode();
        node_count = 0;
        root->node_id = node_count++;
        for (auto &v: m.vertices()){
//            cout << m.points()[v].x() << " " << m.points()[v].y() << " " << m.points()[v].z() << endl;
            double x = m.points()[v].x(), y = m.points()[v].y();
            if (Base::doubleCmp(x - root->x_min) < 0) root->x_min = x;
            if (Base::doubleCmp(x - root->x_max) > 0) root->x_max = x;
            if (Base::doubleCmp(y - root->y_min) < 0) root->y_min = y;
            if (Base::doubleCmp(y - root->y_max) > 0) root->y_max = y;
        }
        set<int> intersect_faces = extractIntersectFaces(m, root);
        for (auto fid: intersect_faces){
            auto fd = *(m.faces_begin() + fid);
            for (auto vd: m.vertices_around_face(m.halfedge(fd))){
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

    void quadTree::buildLevel(Base::Mesh &m, map<int, vector<int> > &face_point_map) {
        vector<treeNode*> last_level_nodes = level_nodes[level];
        vector<treeNode*> cur_level_nodes;
        level++;
        for (auto &node: last_level_nodes){
            double pivot_x = 0.5 * (node->x_min + node->x_max),
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

            set<int> faces = extractIntersectFaces(m, NW_son, node);
            for (auto fid: faces){
                auto fd = *(m.faces_begin() + fid);
                for (auto vd: m.vertices_around_face(m.halfedge(fd))){
                    NW_son->boundary_points_id.insert(vd.idx());
                }
            }

            faces = extractIntersectFaces(m, NE_son, node);
            for (auto fid: faces){
                auto fd = *(m.faces_begin() + fid);
                for (auto vd: m.vertices_around_face(m.halfedge(fd))){
                    NE_son->boundary_points_id.insert(vd.idx());
                }
            }

            faces = extractIntersectFaces(m, SW_son, node);
            for (auto fid: faces){
                auto fd = *(m.faces_begin() + fid);
                for (auto vd: m.vertices_around_face(m.halfedge(fd))){
                    SW_son->boundary_points_id.insert(vd.idx());
                }
            }

            faces = extractIntersectFaces(m, SE_son, node);
            for (auto fid: faces){
                auto fd = *(m.faces_begin() + fid);
                for (auto vd: m.vertices_around_face(m.halfedge(fd))){
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

    kSkip::Graph generateSpanner(vector<int> &pid_list, set<WeightedDistanceOracle::nodePair> &node_pairs,
                                          map<int, int> &new_id){
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
            double tmp_d = WeightedDistanceOracle::enhanced_edges[make_pair(u, v)];
            g.addEdge(new_id[u], new_id[v], tmp_d);
            g.addEdge(new_id[v], new_id[u], tmp_d);
        }
        return g;
    }

    pair<double, bool> queryWSPD(Base::Mesh &m, quadTree &tree, int sid, int tid,
                                 set<WeightedDistanceOracle::nodePair> &node_pairs,
                                 WeightedDistanceOracle::PartitionTree &partition_tree,
                                 map<int, int> &new_id){
        auto box_s = tree.root, box_t = tree.root;
        auto point_s = m.points()[*(m.vertices_begin() + sid)];
        auto point_t = m.points()[*(m.vertices_begin() + tid)];
        //find the leaf contains s and t
        while (box_s->sons.size() > 0){
            for (auto son: box_s->sons){
                if (Base::doubleCmp(point_s.x() - son->x_min) >= 0 && Base::doubleCmp(point_s.x() - son->x_max) <= 0 &&
                    Base::doubleCmp(point_s.y() - son->y_min) >= 0 && Base::doubleCmp(point_s.y() - son->y_max) <= 0){
                    box_s = son;
                    break;
                }
            }
        }
        while (box_t->sons.size() > 0){
            for (auto son: box_t->sons){
                if (Base::doubleCmp(point_t.x() - son->x_min) >= 0 && Base::doubleCmp(point_t.x() - son->x_max) <= 0 &&
                    Base::doubleCmp(point_t.y() - son->y_min) >= 0 && Base::doubleCmp(point_t.y() - son->y_max) <= 0){
                    box_t = son;
                    break;
                }
            }
        }

        if (box_s->node_id == box_t->node_id){
            return make_pair(kSkip::dijkstra(kSkip::my_base_graph, sid, tid).first, box_s->node_id == box_t->node_id);
        }
        else{
            vector<double> ds, dt;
            kSkip::covered_dijkstra(kSkip::my_base_graph, sid, box_s->boundary_points_id, ds);
            kSkip::covered_dijkstra(kSkip::my_base_graph, tid, box_t->boundary_points_id, dt);

            double dis_min = Base::unreachable;
            for (auto pid1: box_s->boundary_points_id){
                for (auto pid2: box_t->boundary_points_id){
                    vector<WeightedDistanceOracle::PartitionTreeNode*> As, At;
                    partition_tree.getPathToRoot(partition_tree.level_nodes[partition_tree.max_level][new_id[pid1]], As);
                    partition_tree.getPathToRoot(partition_tree.level_nodes[partition_tree.max_level][new_id[pid2]], At);
                    double cur_dis = ds[pid1] + WeightedDistanceOracle::distanceQueryBf(node_pairs, As, At) + dt[pid2];
                    if (Base::doubleCmp(cur_dis - dis_min) < 0) dis_min = cur_dis;
                }
            }
            return make_pair(dis_min, box_s->node_id == box_t->node_id);
        }
    }

    void distancePreprocessing(quadTree &quad_tree, kSkip::Graph &base_graph, map<int, vector<int> > &face_point_map){
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
                vector<double> d;
                kSkip::covered_dijkstra(base_graph, s, covered_point, d);
                for (auto pid: covered_point){
                    LQT_distance_map[make_pair(s, pid)] = d[pid];
                    LQT_distance_map[make_pair(pid, s)] = d[pid];
                }
            }

        }
    }

    //dealing with V2V queries
    pair<double, bool> querySpanner(Base::Mesh &m, kSkip::Graph spanner, int sid, int tid, quadTree &quad_tree,
                                    map<int, int> &new_id, WeightedDistanceOracle::PartitionTree &tree, Base::AABB_tree &aabb_tree,
                                    set<WeightedDistanceOracle::nodePair> &node_pairs, bool check_first_flag = 0){
        auto box_s = quad_tree.root, box_t = quad_tree.root;
        auto point_s = m.points()[*(m.vertices_begin() + sid)];
        auto point_t = m.points()[*(m.vertices_begin() + tid)];

        //find the leaf contains s and t
        while (box_s->sons.size() > 0){
            for (auto son: box_s->sons){
                if (Base::doubleCmp(point_s.x() - son->x_min) >= 0 && Base::doubleCmp(point_s.x() - son->x_max) <= 0 &&
                    Base::doubleCmp(point_s.y() - son->y_min) >= 0 && Base::doubleCmp(point_s.y() - son->y_max) <= 0){
                    box_s = son;
                    break;
                }
            }
        }
        while (box_t->sons.size() > 0){
            for (auto son: box_t->sons){
                if (Base::doubleCmp(point_t.x() - son->x_min) >= 0 && Base::doubleCmp(point_t.x() - son->x_max) <= 0 &&
                    Base::doubleCmp(point_t.y() - son->y_min) >= 0 && Base::doubleCmp(point_t.y() - son->y_max) <= 0){
                    box_t = son;
                    break;
                }
            }
        }

//        cout << "box_s = " << box_s->node_id << " box_t = " << box_t->node_id << endl;
//        cout << "box_s son = " << box_s->sons.size() << " box_t son = " << box_t->sons.size() << endl;

        double min_ds = 1e60, min_dt = 1e60;
        int closest_pid_s = -1, closest_pid_t = -1;
        if (box_s->node_id != box_t->node_id){
            int final_s, final_t;

            if (new_id.find(sid) == new_id.end()){
                final_s = spanner.addVertex();

                for (auto pid: box_s->boundary_points_id) {
                    if (LQT_distance_map.find(make_pair(sid, pid)) == LQT_distance_map.end()) {
                        cout << "ERROR: distance not found in LQT_distance_map !!!" << endl;
                        cout << sid << "--" << pid << endl;
                    }
                    spanner.addEdge(final_s, new_id[pid], LQT_distance_map[make_pair(sid, pid)]);

                }
            }
            else{
                final_s = new_id[sid];
                min_ds = 0;
                closest_pid_s = sid;
            }
            if (new_id.find(tid) == new_id.end()){
                final_t = spanner.addVertex();

                for (auto pid: box_t->boundary_points_id){
                    if (LQT_distance_map.find(make_pair(tid, pid)) == LQT_distance_map.end()) {
                        cout << "ERROR: distance not found in LQT_distance_map !!!" << endl;
                    }
                    spanner.addEdge(new_id[pid], final_t, LQT_distance_map[make_pair(tid, pid)]);

                }
            }
            else{
                final_t = new_id[tid];
                min_dt = 0;
                closest_pid_t = tid;
            }

            double ret_dis = -1.0;

            if (Base::doubleCmp(ret_dis) < 0){
                ret_dis = kSkip::dijkstra(spanner, final_s, final_t).first;
            }
            else{
                WSPD_hit++;
            }

            return make_pair(ret_dis, box_s->node_id == box_t->node_id);
        }
        else{
            // in the same box, just SE_A2A dijkstra untill find the result.
            return make_pair(kSkip::dijkstra(kSkip::my_base_graph, sid, tid).first, box_s->node_id == box_t->node_id);
        }
    }

    pair<double, bool> queryA2A(Base::Mesh &m, kSkip::Graph &spanner,
                                kSkip::Graph &base_graph, map<int, vector<int> > &face_point_map,
                                WeightedDistanceOracle::PartitionTree &tree, set<WeightedDistanceOracle::nodePair> &node_pairs,
                                map<int, Base::Point> &point_location_map,
                                Base::Point s, int fid_s,
                                Base::Point t, int fid_t,
                                quadTree &quad_tree, map<int, int> &new_id){

        auto box_s = quad_tree.root, box_t = quad_tree.root;
        //find the leaf contains s and t

        while (box_s->sons.size() > 0){
            for (auto son: box_s->sons){
                if (Base::doubleCmp(s.x() - son->x_min) >= 0 && Base::doubleCmp(s.x() - son->x_max) <= 0 &&
                    Base::doubleCmp(s.y() - son->y_min) >= 0 && Base::doubleCmp(s.y() - son->y_max) <= 0){
                    box_s = son;
                    break;
                }
            }
        }
        while (box_t->sons.size() > 0){
            for (auto son: box_t->sons){
                if (Base::doubleCmp(t.x() - son->x_min) >= 0 && Base::doubleCmp(t.x() - son->x_max) <= 0 &&
                    Base::doubleCmp(t.y() - son->y_min) >= 0 && Base::doubleCmp(t.y() - son->y_max) <= 0){
                    box_t = son;
                    break;
                }
            }
        }

//        cout << "box_id = " << box_s->node_id << " " << box_t->node_id << endl;

        if (box_s->node_id != box_t->node_id){

//            if (box_s->boundary_points_id.size() * box_t->boundary_points_id.size() > 0 &&
//                    box_s->boundary_points_id.size() * box_t->boundary_points_id.size() < 2000){
//                double res = Base::unreachable;
//                map<int, double> d_1, d_2;
//
//                for (auto bpid1: box_s->boundary_points_id) {
//                    double dt = Base::unreachable;
//                    for (auto pid1: face_point_map[fid_s]) {
//                        double t_dis = sqrt(CGAL::squared_distance(s, point_location_map[pid1])) +
//                                      LQT_distance_map[make_pair(pid1, bpid1)];
//                        if (Base::doubleCmp(t_dis - dt) < 0) dt = t_dis;
//                    }
//                    d_1[bpid1] = dt;
//                }
//
//                for (auto bpid2: box_t->boundary_points_id) {
//                    double dt = Base::unreachable;
//                    for (auto pid2: face_point_map[fid_t]) {
//                        double t_dis = sqrt(CGAL::squared_distance(t, point_location_map[pid2])) +
//                                      LQT_distance_map[make_pair(pid2, bpid2)];
//                        if (Base::doubleCmp(t_dis - dt) < 0) dt = t_dis;
//                    }
//                    d_2[bpid2] = dt;
//                }
//
//                vector<WeightedDistanceOracle::PartitionTreeNode*> leaf_nodes(tree.level_nodes[tree.max_level].begin(), tree.level_nodes[tree.max_level].end());
//                vector<WeightedDistanceOracle::PartitionTreeNode *> As, At;
//
//                for (auto bpid1: box_s->boundary_points_id){
//                    for (auto bpid2: box_t->boundary_points_id){
//                        As.clear(); At.clear();
//                        tree.getPathToRoot(leaf_nodes[new_id[bpid1]], As);
//                        tree.getPathToRoot(leaf_nodes[new_id[bpid2]], At);
//                        cout << "check: " << leaf_nodes[new_id[bpid1]]->center_idx << " | " << bpid1 << endl;
//                        cout << "check: " << leaf_nodes[new_id[bpid2]]->center_idx << " | " << bpid2 << endl;
//                        cout << "check: " << d_1[bpid1] << " | " << d_2[bpid2] << endl;
//
//                        double t_dis = d_1[bpid1] + WeightedDistanceOracle::distanceQueryEfficient(node_pairs, As, At) + d_2[bpid2];
//                        if (Base::doubleCmp(t_dis - res) < 0) res = t_dis;
//                    }
//                }
//                return make_pair(res, box_s->node_id == box_t->node_id);
//
//            }

            int V_flag = spanner.num_V, E_flag = spanner.num_E; //  backup

            auto final_s = spanner.addVertex();

            for (auto bpid: box_s->boundary_points_id){
                double d = Base::unreachable;
                for (auto pid: face_point_map[fid_s]){
//                    double t_dis = sqrt(CGAL::squared_distance(s, point_location_map[pid])) + LQT_distance_map[make_pair(closest_pid, bpid)];
                    double t_dis = sqrt(CGAL::squared_distance(s, point_location_map[pid])) + LQT_distance_map[make_pair(pid, bpid)];
                    if (Base::doubleCmp(t_dis - d) < 0) d = t_dis;
                }
                spanner.addEdge(final_s, new_id[bpid], d);
            }
            auto final_t = spanner.addVertex();

            for (auto bpid: box_t->boundary_points_id){
                double d = Base::unreachable;
                for (auto pid: face_point_map[fid_t]){
                    double t_dis = sqrt(CGAL::squared_distance(t, point_location_map[pid])) + LQT_distance_map[make_pair(pid, bpid)];
                    if (Base::doubleCmp(t_dis - d) < 0) d = t_dis;
                }
                spanner.addEdge(new_id[bpid], final_t, d);
            }
            double res = kSkip::dijkstra(spanner, final_s, final_t).first;

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
                double dis = CGAL::squared_distance(s, point_location_map[pid]);
                base_graph.addEdge(sid, pid, sqrt(dis));
            }
            auto tid = base_graph.addVertex();
            for (auto pid: face_point_map[fid_t]){
                double dis = CGAL::squared_distance(t, point_location_map[pid]);
                base_graph.addEdge(pid, tid, sqrt(dis));
            }
            double res = kSkip::dijkstra(base_graph, sid, tid).first;

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

    pair<double, bool> queryA2A(Base::Mesh &m, kSkip::Graph &spanner,
                                kSkip::Graph &base_graph, map<int, vector<int> > &face_point_map,
                                map<int, Base::Point> &point_location_map,
                                Base::Point s, int fid_s,
                                Base::Point t, int fid_t,
                                quadTree &quad_tree, map<int, int> &new_id){

        auto box_s = quad_tree.root, box_t = quad_tree.root;
        //find the leaf contains s and t

        while (box_s->sons.size() > 0){
            for (auto son: box_s->sons){
                if (Base::doubleCmp(s.x() - son->x_min) >= 0 && Base::doubleCmp(s.x() - son->x_max) <= 0 &&
                    Base::doubleCmp(s.y() - son->y_min) >= 0 && Base::doubleCmp(s.y() - son->y_max) <= 0){
                    box_s = son;
                    break;
                }
            }
        }
        while (box_t->sons.size() > 0){
            for (auto son: box_t->sons){
                if (Base::doubleCmp(t.x() - son->x_min) >= 0 && Base::doubleCmp(t.x() - son->x_max) <= 0 &&
                    Base::doubleCmp(t.y() - son->y_min) >= 0 && Base::doubleCmp(t.y() - son->y_max) <= 0){
                    box_t = son;
                    break;
                }
            }
        }


        if (box_s->node_id != box_t->node_id){

            int V_flag = spanner.num_V, E_flag = spanner.num_E; //  backup

            auto final_s = spanner.addVertex();

            for (auto bpid: box_s->boundary_points_id){
                double d = Base::unreachable;
                for (auto pid: face_point_map[fid_s]){
                    double t_dis = sqrt(CGAL::squared_distance(s, point_location_map[pid])) + LQT_distance_map[make_pair(pid, bpid)];
                    if (Base::doubleCmp(t_dis - d) < 0) d = t_dis;
                }
                spanner.addEdge(final_s, new_id[bpid], d);
            }
            auto final_t = spanner.addVertex();

            for (auto bpid: box_t->boundary_points_id){
                double d = Base::unreachable;
                for (auto pid: face_point_map[fid_t]){
                    double t_dis = sqrt(CGAL::squared_distance(t, point_location_map[pid])) + LQT_distance_map[make_pair(pid, bpid)];
                    if (Base::doubleCmp(t_dis - d) < 0) d = t_dis;
                }
                spanner.addEdge(new_id[bpid], final_t, d);
            }
            double res = kSkip::dijkstra(spanner, final_s, final_t).first;

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
                double dis = CGAL::squared_distance(s, point_location_map[pid]);
                base_graph.addEdge(sid, pid, sqrt(dis));
            }
            auto tid = base_graph.addVertex();
            for (auto pid: face_point_map[fid_t]){
                double dis = CGAL::squared_distance(t, point_location_map[pid]);
                base_graph.addEdge(pid, tid, sqrt(dis));
            }
            double res = kSkip::dijkstra(base_graph, sid, tid).first;

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
