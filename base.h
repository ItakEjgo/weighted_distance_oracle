//
// Created by huang on 2021/7/6.
//
//#pragma once
#ifndef CODE_BASE_H
#define CODE_BASE_H
//#define PrintDetails

#include <bits/stdc++.h>
//#include "getopt.h"
#include "chrono"
#include "assert.h"
#include <random>

//CGAL headers
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>



//boost headers
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

//Base contains the basic functions
namespace Base {

    using namespace std;
    using namespace chrono;

    using Kernel = CGAL::Simple_cartesian<double>;
//    using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
    using Point = Kernel::Point_3;
    using Vector = Kernel::Vector_3;
    using Triangle = Kernel::Triangle_3;
    using Plane = Kernel::Plane_3;
    using Ray = Kernel::Ray_3;

    using Line = Kernel::Line_3;
    using Segment = Kernel::Segment_3;
    using Mesh = CGAL::Surface_mesh<Point>;
    using Location_property = Mesh::Property_map<Mesh::Vertex_index, Point>;
    using Aff_transformation_3 = CGAL::Aff_transformation_3<Kernel>;

    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, boost::property<boost::edge_weight_t, double> >;
    using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
//    using Traits = CGAL::Surface_mesh_shortest_path_traits<CGAL::Exact_predicates_exact_constructions_kernel, Mesh>;
    using Traits = CGAL::Surface_mesh_shortest_path_traits<Kernel, Mesh>;
    using Surface_mesh_shortest_path = CGAL::Surface_mesh_shortest_path<Traits>;

    using AABB_face_graph_primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
    using AABB_face_graph_traits = CGAL::AABB_traits<Kernel, AABB_face_graph_primitive>;
    using AABB_tree = CGAL::AABB_tree<AABB_face_graph_traits>;

    typedef boost::optional<AABB_tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;


    const double eps = 1e-7;
    const double PI = acos(-1.0);
    const double unreachable = numeric_limits<double>::max();

    vector<pair<Point, Point> > A2A_query;
    vector<pair<unsigned, unsigned> > A2A_fid;


    // memory used in KB.
    size_t physical_memory_used_by_process() {
        FILE *file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != nullptr) {
            if (strncmp(line, "VmRSS:", 6) == 0) {
                int len = strlen(line);

                const char *p = line;
                for (; std::isdigit(*p) == false; ++p) {}

                line[len - 3] = 0;
                result = atoi(p);

                break;
            }
        }

        fclose(file);

        return result;
    }

    int doubleCmp(const double &x) {
        if (fabs(x) < eps) return 0;
        return x < 0 ? -1 : 1;
    }

    // Return the point rotate around a vector with given angle.
    Point pointRotationAroundVector(Point p, Point vec_source, Point vec_target, double theta);

    double distanceSnell(Mesh &mesh, vector<double> &face_weight, Point p1, int fid1, Point p2, int fid2, int eid);

    //get the boundary of a given grid. The order of the returned values are x_min, x_max, y_min, y_max;
    vector<double> fastRetrieveGridBoundary(const int &grid_id, const int &grid_num, const double &mesh_xmin, const double &mesh_xmax,
                                            const double &mesh_ymin, const double &mesh_ymax){
        vector<double> ret = {};
        int side_len = floor(sqrt(grid_num) + eps);
        double delta_x = (mesh_xmax - mesh_xmin) / side_len,
            delta_y = (mesh_ymax - mesh_ymin) / side_len;
        ret.emplace_back(mesh_xmin + (grid_id % side_len) * delta_x);
        ret.emplace_back(mesh_xmin + (grid_id % side_len + 1) * delta_x);
        ret.emplace_back(mesh_ymin + floor(grid_id / side_len + eps) * delta_y);
        ret.emplace_back(mesh_ymin + floor(grid_id / side_len + 1 + eps) * delta_y);
        return ret;
    }

    //get the boundary of the given mesh. The order of the returned values are x_min, x_max, y_min, y_max;
    vector<double> retrieveMeshBoundary(const Mesh &mesh){
        double x_min = 1e60, x_max = -1e60, y_min = 1e60, y_max = -1e60;
        for (auto vd: mesh.vertices()){
            double x = mesh.points()[vd].x(),
                    y = mesh.points()[vd].y();
            x_min = doubleCmp(x - x_min) < 0 ? x : x_min;
            x_max = doubleCmp(x - x_max) > 0 ? x : x_max;
            y_min = doubleCmp(y - y_min) < 0 ? y : y_min;
            y_max = doubleCmp(y - y_max) > 0 ? y : y_max;
        }
        vector<double> ret = {x_min, x_max, y_min, y_max};
        return ret;
    }

    void getOpt2(int argc, char **argv, int &generate_queries, string &file_name, int &weight, int &q_num,
                 double &eps, int &sp_num, int &algo_type, int &lqt_lev, string &output_file) {
//        struct option long_options[] = {
//                {"generate",  required_argument, 0, 'g'},
//                {"mesh",      required_argument, 0, 'm'},
//                {"weight", required_argument, 0, 'w'},
//                {"query_num", required_argument, 0, 'q'},
//                {"algo",      required_argument, 0, 'a'},
//                {"eps",       required_argument, 0, 'e'},
//                {"sp_num",    required_argument, 0, 's'},
//                {"lqt_lev",   required_argument, 0, 'l'},
//                {"output",    required_argument, 0, 'o'},
//                {0, 0,                           0, 0}
//        };
//        int opt;
//        const char *opt_string = "g:m:w:q:a:e:s:l:o:";
//        int option_index = 0;
//        generate_queries = 0;
//        lqt_lev = 0;
//        while ((opt = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
//            char ch = (char) opt;
//            switch (ch) {
//                case 'g':
//                    generate_queries = atoi(optarg);
//                    break;
//                case 'm':
//                    file_name = optarg;
//                    break;
//                case 'w':
//                    weight = atoi(optarg);
//                    break;
//                case 'q':
//                    q_num = atoi(optarg);
//                    break;
//                case 'a':
//                    algo_type = atoi(optarg);
//                    break;
//                case 'e':
//                    eps = atof(optarg);
//                    break;
//                case 's':
//                    sp_num = atoi(optarg);
//                    break;
//                case 'l':
//                    lqt_lev = atoi(optarg);
//                    break;
//                case 'o':
//                    output_file = optarg;
//                    break;
//                default:
//                    break;
//            }
//        }
    }

    double distanceSnell(Mesh &mesh, vector<double> &face_weight, Point p1, int fid1, Point p2, int fid2, int eid) {
        //  same face
        if (!(fid1 - fid2)) {
            return face_weight[fid1] * sqrt(CGAL::squared_distance(p1, p2));
        }
            //  same point
        else if (!doubleCmp(p1.x() - p2.x()) && !doubleCmp(p1.y() - p2.y()) && !doubleCmp(p1.z() - p2.z())) {
            return 0.0;
        }
        //  different points on different faces
        auto fd1 = *(mesh.faces().begin() + fid1), fd2 = *(mesh.faces().begin() + fid2);
        vector<Point> points_f1 = {}, points_f2 = {};
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd1))) {
            points_f1.emplace_back(mesh.points()[vd]);
        }
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd2))) {
            points_f2.emplace_back(mesh.points()[vd]);
        }
        auto norm1 = CGAL::normal(points_f1[0], points_f1[1], points_f1[2]);
        auto norm2 = CGAL::normal(points_f2[0], points_f2[1], points_f2[2]);

        double theta = PI * CGAL::approximate_angle(norm1, norm2) / 180;
        auto common_ed = *(mesh.edges().begin() + eid);
        Point vec_source = mesh.points()[mesh.source(mesh.halfedge(common_ed))];
        Point vec_target = mesh.points()[mesh.target(mesh.halfedge(common_ed))];

        Point p2_rot1 = pointRotationAroundVector(p2, vec_source, vec_target, theta);
        Point p2_rot2 = pointRotationAroundVector(p2, vec_source, vec_target, 2 * PI - theta);

        Point p2_rot; // the further one is the rotated one in the same plane.
        if (doubleCmp(CGAL::squared_distance(p1, p2_rot1) - CGAL::squared_distance(p1, p2_rot2)) > 0) {
            p2_rot = p2_rot1;
        } else {
            p2_rot = p2_rot2;
        }

        Segment common_edge(vec_source, vec_target);
        Line common_edge_line(common_edge.supporting_line());
        Point point_l(common_edge_line.projection(p1)), point_r(common_edge_line.projection(p2_rot));

        double dis1, dis2;
        if (!common_edge.has_on(point_l)) {
            dis1 = CGAL::squared_distance(point_l, common_edge.source());
            dis2 = CGAL::squared_distance(point_l, common_edge.target());
            if (dis1 < dis2) {
                point_l = common_edge.source();
            } else {
                point_l = common_edge.target();
            }
        }
        if (!common_edge.has_on(point_r)) {
            dis1 = CGAL::squared_distance(point_r, common_edge.source());
            dis2 = CGAL::squared_distance(point_r, common_edge.target());
            if (dis1 < dis2) {
                point_r = common_edge.source();
            } else {
                point_r = common_edge.target();
            }
        }
        const double n1 = face_weight[fid1], n2 = face_weight[fid2];
        const double ratio = n2 / n1;
        Point point_mid;

        for (auto i = 0; i != 50; i++) {
            Vector t_vec(point_r - point_l);
            point_mid = point_l + t_vec * 0.5;
            double sin_w1 =
                    sqrt(CGAL::squared_distance(point_l, point_mid)) / sqrt(CGAL::squared_distance(p1, point_mid));
            double sin_w2 =
                    sqrt(CGAL::squared_distance(point_r, point_mid)) / sqrt(CGAL::squared_distance(p2_rot, point_mid));
            if (doubleCmp(sin_w1 - sin_w2 * ratio) >= 0) {
                point_l = point_mid;
            } else {
                point_r = point_mid;
            }
        }
        return n1 * sqrt(CGAL::squared_distance(p1, point_mid)) +
               n2 * sqrt(CGAL::squared_distance(p2_rot, point_mid));
    }

    Point pointRotationAroundVector(Point p, Point vec_source, Point vec_target, double theta) {
        Vector rotation_vec(vec_target - vec_source);
        double rotation_vec_len = sqrt(rotation_vec.squared_length());
        Vector axis(rotation_vec / rotation_vec_len);
        double u = axis.x(), v = axis.y(), w = axis.z();
        double a = vec_source.x(), b = vec_source.y(), c = vec_source.z();
        double M[3][4];
        M[0][0] = u * u + (v * v + w * w) * cos(theta);
        M[0][1] = u * v * (1 - cos(theta)) - w * sin(theta);
        M[0][2] = u * w * (1 - cos(theta)) + v * sin(theta);
        M[0][3] = (a * (v * v + w * w) - u * (b * v + c * w)) * (1 - cos(theta)) + (b * w - c * v) * sin(theta);
        M[1][0] = u * v * (1 - cos(theta)) + w * sin(theta);
        M[1][1] = v * v + (u * u + w * w) * cos(theta);
        M[1][2] = v * w * (1 - cos(theta)) - u * sin(theta);
        M[1][3] = (b * (u * u + w * w) - v * (a * u + c * w)) * (1 - cos(theta)) + (c * u - a * w) * sin(theta);
        M[2][0] = u * w * (1 - cos(theta)) - v * sin(theta);
        M[2][1] = v * w * (1 - cos(theta)) + u * sin(theta);
        M[2][2] = w * w + (u * u + v * v) * cos(theta);
        M[2][3] = (c * (u * u + v * v) - w * (a * u + b * v)) * (1 - cos(theta)) + (a * v - b * u) * sin(theta);
        Aff_transformation_3 trans(
                M[0][0], M[0][1], M[0][2], M[0][3],
                M[1][0], M[1][1], M[1][2], M[1][3],
                M[2][0], M[2][1], M[2][2], M[2][3]
        );
        return trans.transform(p);
    }

    void dijkstra_SSAD(Graph &g, int source_idx, vector<double> &d, vector<int> &fa) {
        d.resize(num_vertices(g));
        vector<vertex_descriptor> p(num_vertices(g));
        auto s = *(vertices(g).first + source_idx);
        boost::dijkstra_shortest_paths(g, s,
                                       predecessor_map(boost::make_iterator_property_map(p.begin(),
                                                                                         get(boost::vertex_index, g)))
                                               .distance_map(boost::make_iterator_property_map(d.begin(),
                                                                                               get(boost::vertex_index,
                                                                                                   g))));
        fa.resize(num_vertices(g));
        for (auto i = 0; i != num_vertices(g); i++) {
            fa[i] = p[i];
        }
    }

    // generate face weights for each face.
    vector<double> generateFaceWeight(const int &num_faces){
        uniform_real_distribution<double> gen(1.000001, 5);
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine e(seed);
        vector<double> ret = {};
        for (auto i = 0; i != num_faces; i++) ret.emplace_back(gen(e));
        return ret;
    }

    vector<double> generateFaceWeight(string &file_name) {
        Mesh surface_mesh;
        ifstream fin(file_name);
        fin >> surface_mesh;
        auto n = surface_mesh.num_faces();
        uniform_real_distribution<double> gen(1.000001, 5);
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine e(seed);
        vector<double> face_weight = {};
        for (auto i = 0; i < n; i++) {
            face_weight.emplace_back(gen(e));
        }
        return face_weight;
    }

    // return an arbitrary surface point and its face location within a given region.
    pair<Point, unsigned> generateArbitrarySurfacePoint(const Mesh &mesh, const AABB_tree &aabb_tree,
                                        const double &x_min, const double &x_max,
                                        const double &y_min, const double &y_max){
        unsigned seed = system_clock::now().time_since_epoch().count();
        uniform_real_distribution<double> gen_x(x_min, x_max);
        uniform_real_distribution<double> gen_y(y_min, y_max);
        default_random_engine e(seed);

        Point pp;
        while(1){
            double x_rand = gen_x(e), y_rand = gen_y(e);
            Point bot(x_rand, y_rand, 0), top(x_rand, y_rand, 1), top2(x_rand, y_rand, -1);
            Ray ray(bot, top);
            Ray_intersection intersection = aabb_tree.first_intersection(ray);
            if (intersection){
                if (boost::get<Point>(&(intersection->first))){
                    const Point* p = boost::get<Point>(&(intersection->first));
                    pp = *p;
                    break;
                }
            }
            else{
                Ray ray2(bot, top2);
                intersection = aabb_tree.first_intersection(ray2);
                if (intersection){
                    if (boost::get<Point>(&(intersection->first))){
                        const Point* p = boost::get<Point>(&(intersection->first));
                        pp = *p;
                        break;
                    }
                }
            }
        }
        auto p_location = CGAL::Polygon_mesh_processing::locate_with_AABB_tree(pp, aabb_tree, mesh);
        return make_pair(pp, p_location.first.idx());
    }

    // generate A2A queries, if inner_flag is set, the queries will be inner-box queries.
    bool generateQueriesA2A(const Mesh &mesh, const vector<double> &mesh_boundary, const AABB_tree &aabb_tree,
                            const unsigned &q_num, const unsigned &grid_num,
                            const bool &inner_flag){
        unsigned seed = system_clock::now().time_since_epoch().count();
        default_random_engine e(seed);
        uniform_int_distribution<unsigned> gen_grid_id(0, grid_num - 1);
        for (auto i = 0; i != q_num; i++){
            unsigned grid_id1 = gen_grid_id(e);
            unsigned grid_id2 = grid_id1;
            do{
                if (inner_flag) break;
                grid_id2 = gen_grid_id(e);
            } while (grid_id1 == grid_id2);

            vector<double> grid1_boundary = fastRetrieveGridBoundary(grid_id1, grid_num, mesh_boundary[0], mesh_boundary[1], mesh_boundary[2], mesh_boundary[3]);
            vector<double> grid2_boundary = fastRetrieveGridBoundary(grid_id2, grid_num, mesh_boundary[0], mesh_boundary[1], mesh_boundary[2], mesh_boundary[3]);

            auto s = generateArbitrarySurfacePoint(mesh, aabb_tree, grid1_boundary[0], grid1_boundary[1], grid1_boundary[2], grid1_boundary[3]),
                t = generateArbitrarySurfacePoint(mesh, aabb_tree, grid2_boundary[0], grid2_boundary[1], grid2_boundary[2], grid2_boundary[3]);
            A2A_query.emplace_back(s.first, t.first);
            A2A_fid.emplace_back(s.second, t.second);
        }
    }


    pair<Point, int> generateArbitrarySurfacePoint(Mesh &mesh) {
        auto fid = rand() % mesh.num_faces();
//        cout << "fid = " << fid << endl;
        auto fd = *(mesh.faces().begin() + fid);

        uniform_real_distribution<double> gen(0, 1);
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine e(seed);
        double a = gen(e), b = gen(e), c = gen(e);
        vector<Point> p;
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd))) {
            p.push_back(mesh.points()[vd]);
        }
        Point A = p[0] + (p[1] - p[0]) * a, B = p[0] + (p[2] - p[0]) * b;
        Point ret_point = A + (B - A) * c;

        return make_pair(ret_point, fid);
    }

    vector<pair<Point, Point> > generateQueriesA2A(string &file_name, int q_num, vector<pair<int, int> > &A2A_fid) {
        vector<pair<Point, Point> > A2A_query = {};
        A2A_fid.clear();
        srand((int) time(0));
        Mesh surface_mesh;
        ifstream fin(file_name);
        cout << "mesh file = " << file_name << endl;
        fin >> surface_mesh;
        cout << "V= " << surface_mesh.num_vertices() << " E= " << surface_mesh.num_edges() << " F= " << surface_mesh.num_faces() << endl;
        AABB_tree aabb_tree;
        CGAL::Polygon_mesh_processing::build_AABB_tree(surface_mesh, aabb_tree);
        double x_min = 1e60, x_max = -1e60, y_min = 1e60, y_max = -1e60;
        for (auto vd: surface_mesh.vertices()){
            double x = surface_mesh.points()[vd].x(),
                y = surface_mesh.points()[vd].y();
            if (doubleCmp(x - x_min) < 0){
                x_min = x;
            }
            if (doubleCmp(x - x_max) > 0){
                x_max = x;
            }
            if (doubleCmp(y - y_min) < 0){
                y_min = y;
            }
            if (doubleCmp(y - y_max) > 0){
                y_max = y;
            }
        }
        int cnt = 0;
        for (auto i = 0; i < q_num; i++) {
            auto p1_pair = generateArbitrarySurfacePoint(surface_mesh);
            auto p2_pair = generateArbitrarySurfacePoint(surface_mesh);
            A2A_query.emplace_back(p1_pair.first, p2_pair.first);
            A2A_fid.emplace_back(p1_pair.second, p2_pair.second);
        }
        return A2A_query;
    }

    void loadQueriesA2A(vector<pair<Point, Point> > &A2A_query, vector<pair<int, int> > &A2A_fid) {
        ifstream fin("A2A.query");
        Point p1, p2;
        int fid1, fid2;
        A2A_query.clear();
        A2A_fid.clear();
        while (fin >> p1 >> fid1 >> p2 >> fid2) {
            A2A_query.emplace_back(p1, p2);
            A2A_fid.emplace_back(fid1, fid2);
        }
    }

    void loadFaceWeight(vector<double> &face_weight) {
        face_weight.clear();
        ifstream fin("face_weight.query");
        double w;
        while (fin >> w) {
            face_weight.emplace_back(w);
        }
    }

}

#endif //CODE_BASE_H
