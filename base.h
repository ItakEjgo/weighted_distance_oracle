//
// Created by huang on 2021/7/6.
//

#ifndef CODE_BASE_H
#define CODE_BASE_H

//#define PrintDetails

#include <bits/stdc++.h>
#include "getopt.h"
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


//boost headers
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

//Base contains the basic functions
namespace Base{

    using namespace std;
    using namespace chrono;

    using Kernel = CGAL::Simple_cartesian<double>;
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

    const double eps = 1e-6;
    const double PI = acos(-1.0);
    const double unreachable = numeric_limits<double>::max();


    int doubleCmp(const double &x){
        if (fabs(x) < eps) return 0;
        return x < 0 ? -1 : 1;
    }
    

    void getOpt(int argc, char **argv, string &file_name, double &eps, int &type, int &point_num);

    // Return the point rotate around a vector with given angle.
    Point pointRotationAroundVector(Point p, Point vec_source, Point vec_target, double theta);

    double distanceSnell(Mesh &mesh, vector<double> &face_weight, Point p1, int fid1, Point p2, int fid2, int eid);

    void getOpt(int argc, char **argv, string &file_name, double &eps, int &type, int &point_num, int &level) {
        struct option long_options[] = {
                {"mesh", required_argument, 0, 'm'},
                {"eps", required_argument, 0, 'e'},
                {"type", required_argument, 0, 't'},
                {"pointnum", required_argument, 0, 'p'},
                {"level", required_argument, 0, 'l'},
                {0, 0, 0, 0}
        };
        int opt;
        const char* opt_string = "m:e:t:p:l:";
        int option_index = 0;
        while ((opt = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1){
            char ch = (char)opt;
            switch (ch){
                case 'm': file_name = optarg; break;
                case 'e': eps = atof(optarg); break;
                case 't': type = atoi(optarg); break;
                case 'p': point_num = atoi(optarg); break;
                case 'l': level = atoi(optarg); break;
                default: break;
            }
        }
    }

    double distanceSnell(Mesh &mesh, vector<double> &face_weight, Point p1, int fid1, Point p2, int fid2, int eid) {
        //  same face
        if (!(fid1 - fid2)){
            return face_weight[fid1] * sqrt(CGAL::squared_distance(p1, p2));
        }
            //  same point
        else if(!doubleCmp(p1.x() - p2.x()) && !doubleCmp(p1.y() - p2.y()) && !doubleCmp(p1.z() - p2.z())){
            return 0.0;
        }
        //  different points on different faces
        auto fd1 = *(mesh.faces().begin() + fid1), fd2 = *(mesh.faces().begin() + fid2);
        vector<Point> points_f1 = {}, points_f2 = {};
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd1))){
            points_f1.emplace_back(mesh.points()[vd]);
        }
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd2))){
            points_f2.emplace_back(mesh.points()[vd]);
        }
        auto norm1 = CGAL::normal(points_f1[0], points_f1[1], points_f1[2]);
        auto norm2 = CGAL::normal(points_f2[0], points_f2[1], points_f2[2]);

        double theta = Base::PI * CGAL::approximate_angle(norm1, norm2) / 180;
        auto common_ed = *(mesh.edges().begin() + eid);
        Point vec_source = mesh.points()[mesh.source(mesh.halfedge(common_ed))];
        Point vec_target = mesh.points()[mesh.target(mesh.halfedge(common_ed))];

        Point p2_rot1 = pointRotationAroundVector(p2, vec_source, vec_target, theta);
        Point p2_rot2 = pointRotationAroundVector(p2, vec_source, vec_target, 2 * PI - theta);

        Point p2_rot; // the further one is the rotated one in the same plane.
        if (doubleCmp(CGAL::squared_distance(p1, p2_rot1) - CGAL::squared_distance(p1, p2_rot2)) > 0){
            p2_rot = p2_rot1;
        }
        else{
            p2_rot = p2_rot2;
        }

        Segment common_edge(vec_source, vec_target);
        Line common_edge_line(common_edge.supporting_line());
        Point point_l(common_edge_line.projection(p1)), point_r(common_edge_line.projection(p2_rot));

        double dis1, dis2;
        if (!common_edge.has_on(point_l)){
            dis1 = CGAL::squared_distance(point_l, common_edge.source());
            dis2 = CGAL::squared_distance(point_l, common_edge.target());
            if (dis1 < dis2){
                point_l = common_edge.source();
            }
            else{
                point_l = common_edge.target();
            }
        }
        if (!common_edge.has_on(point_r)){
            dis1 = CGAL::squared_distance(point_r, common_edge.source());
            dis2 = CGAL::squared_distance(point_r, common_edge.target());
            if (dis1 < dis2){
                point_r = common_edge.source();
            }
            else{
                point_r = common_edge.target();
            }
        }
        const double n1 = face_weight[fid1], n2 = face_weight[fid2];
        const double ratio = n2 / n1;
        Point point_mid;

        for (auto i = 0; i != 50; i++){
            Vector t_vec(point_r - point_l);
            point_mid = point_l + t_vec * 0.5;
            double sin_w1 = sqrt(CGAL::squared_distance(point_l, point_mid)) / sqrt(CGAL::squared_distance(p1, point_mid));
            double sin_w2 = sqrt(CGAL::squared_distance(point_r, point_mid)) / sqrt(CGAL::squared_distance(p2_rot, point_mid));
            if (doubleCmp(sin_w1 - sin_w2 * ratio) >= 0){
                point_l = point_mid;
            }
            else{
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
        Base::Aff_transformation_3 trans(
                M[0][0],
                M[0][1],
                M[0][2],
                M[0][3],
                M[1][0],
                M[1][1],
                M[1][2],
                M[1][3],
                M[2][0],
                M[2][1],
                M[2][2],
                M[2][3]
        );
        return trans.transform(p);
    }

    void dijkstra_SSAD(Graph &g, int source_idx, vector<double> &d, vector<int> &fa){
        d.resize(num_vertices(g));
        vector<vertex_descriptor> p(num_vertices(g));
        auto s = *(vertices(g).first + source_idx);
        boost::dijkstra_shortest_paths(g, s,
                                       predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g)))
                                       .distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));
        fa.resize(num_vertices(g));
        for (auto i = 0; i != num_vertices(g); i++){
            fa[i] = p[i];
        }
    }

    pair<Point, int> generateArbitrarySurfacePoint(Mesh &mesh, AABB_tree &aabb_tree){
        auto fid = rand() % mesh.num_faces();
//        cout << "fid = " << fid << endl;
        auto fd = *(mesh.faces().begin() + fid);

        uniform_real_distribution<double> gen(0.0, 1.0);
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine e(seed);
        double a = gen(e), b = gen(e), c = gen(e);
        vector<Point> p;
        for (auto vd: mesh.vertices_around_face(mesh.halfedge(fd))){
            p.push_back(mesh.points()[vd]);
        }
        Point A = p[0] + (p[1] - p[0]) * a, B = p[0] + (p[2] - p[0]) * b;
        Point ret_point = A + (B - A) * c;
//        cout << x_gen << " || " << y_gen << " || " << z_cor_min << endl;
        return make_pair(ret_point, fid);
    }
}




#endif //CODE_BASE_H
