//
// Created by huang on 2021/7/6.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H
#define WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H

#include "base.h"

namespace WeightedDistanceOracle{

    using namespace std;

    double getFaceMaxLength(Base::Mesh &mesh, vector<double> &face_max_length);


    double placeBisectorPointsFixed(Base::Mesh &mesh,
                                    const int &point_num,
                                    map<int, vector<int> > &edge_bisector_map,
                                    map<int, vector<int> > &bisector_point_map,
                                    map<int, int> &point_face_map,
                                    map<int, Base::Point> &point_location_map);

    double getFaceMaxLength(Base::Mesh &mesh, vector<double> &face_max_length) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        face_max_length.resize(mesh.num_faces(), 0.0);
        double max_length = -1.0, e_length;
        for (auto fd: mesh.faces()){
            max_length = -1.0;
            for (auto hed: mesh.halfedges_around_face(mesh.halfedge(fd))){
                e_length = sqrt(CGAL::squared_distance(mesh.points()[mesh.source(hed)], mesh.points()[mesh.target(hed)]));
                max_length = Base::doubleCmp(e_length - max_length) > 0 ? e_length : max_length;
            }
            face_max_length[fd.idx()] = max_length;
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

    double placeBisectorPointsFixed(Base::Mesh &mesh, const int &point_num,
                                    map<int, vector<int>> &edge_bisector_map,
                                    map<int, vector<int>> &bisector_point_map,
                                    map<int, int> &point_face_map,
                                    map<int, Base::Point> &point_location_map) {
        auto start_time = chrono::_V2::system_clock::now();  //  timer

        int num_vertices = static_cast<int>(mesh.num_vertices());
        int bisector_num = 0;
        for (auto fd: mesh.faces()){
            for (auto hed:mesh.halfedges_around_face(mesh.halfedge(fd))){
                int eid = static_cast<int>(mesh.edge(hed).idx()),
                        eid2 = static_cast<int>(mesh.edge(mesh.prev(hed)).idx());
                vector<Base::Point> p(3);
                vector<int> pid(3);
                for (auto i = 0; i != 3; hed = mesh.next(hed), i++){
                    p[i] = mesh.points()[mesh.source(hed)];
                    pid[i] = static_cast<int>(mesh.source(hed).idx());
                }
                vector<int> cur_bisector = {};
                cur_bisector.push_back(pid[0]);

                double len1 = sqrt(CGAL::squared_distance(p[1], p[0])),
                        len2 = sqrt(CGAL::squared_distance(p[2], p[0]));
                Base::Point p_end(p[1] + Base::Vector(p[2] - p[1]) * len1 / (len1 + len2));
                Base::Vector vec_bisector(p_end - p[0]);
                double limit_distance = sqrt(CGAL::squared_distance(p[0], p_end));
                double cur_distance = 0;
                int k = 0;
#ifdef PrintDetails
                cout << "current V = " << pid[0] << endl;
#endif
                while (Base::doubleCmp(cur_distance - limit_distance) < 0){
                    Base::Point bisector_p = p[0] + static_cast<double>(++k) / point_num * vec_bisector;
#ifdef PrintDetails
                    cout << "k = " << k << ", p = " << bisector_p << endl;
#endif
                    point_location_map[num_vertices] = bisector_p;
                    point_face_map[num_vertices] = static_cast<int>(fd.idx());
                    cur_bisector.push_back(num_vertices++);
                    cur_distance += limit_distance / point_num;
                }
                bisector_point_map[bisector_num] = cur_bisector;
                edge_bisector_map[eid].push_back(bisector_num);
                edge_bisector_map[eid2].push_back(bisector_num++);
            }
        }

        auto end_time = chrono::_V2::system_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        return static_cast<double>(duration.count());
    }

}

#endif //WEIGHTED_DISTANCE_ORACLE_WEIGHTED_DISTANCE_ORACLE_H
