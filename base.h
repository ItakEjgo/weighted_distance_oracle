//
// Created by huang on 2021/7/6.
//

#ifndef CODE_BASE_H
#define CODE_BASE_H

#define PrintDetails

#include <bits/stdc++.h>
#include "getopt.h"
#include "chrono"

//CGAL headers
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Surface_mesh_shortest_path.h>

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
    using Plane = Kernel::Plane_3;
    using Ray = Kernel::Ray_3;
    using Line = Kernel::Line_3;
    using Segment = Kernel::Segment_3;
    using Mesh = CGAL::Surface_mesh<Point>;
    using Location_property = Mesh::Property_map<Mesh::Vertex_index, Point>;

    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, boost::property<boost::edge_weight_t, double> >;
    using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
    using Traits = CGAL::Surface_mesh_shortest_path_traits<Kernel, Mesh>;
    using Surface_mesh_shortest_path = CGAL::Surface_mesh_shortest_path<Traits>;

    const double eps = 1e-6;
    int doubleCmp(const double &x){
        if (fabs(x) < eps) return 0;
        return x < 0 ? -1 : 1;
    }
    

    void getOpt(int argc, char **argv, string &file_name, double &eps, int &type, int &point_num);

    void getOpt(int argc, char **argv, string &file_name, double &eps, int &type, int &point_num) {
        struct option long_options[] = {
                {"mesh", required_argument, 0, 'm'},
                {"eps", required_argument, 0, 'e'},
                {"type", required_argument, 0, 't'},
                {"pointnum", required_argument, 0, 'p'},
                {0, 0, 0, 0}
        };
        int opt;
        const char* opt_string = "m:e:t:p:";
        int option_index = 0;
        while ((opt = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1){
            char ch = (char)opt;
            switch (ch){
                case 'm': file_name = optarg; break;
                case 'e': eps = atof(optarg); break;
                case 't': type = atoi(optarg); break;
                case 'p': point_num = atoi(optarg); break;
                default: break;
            }
        }
    }


}




#endif //CODE_BASE_H
