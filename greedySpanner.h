//
// Created by huang on 2022/1/23.
//

#ifndef WEIGHTED_DISTANCE_ORACLE_GREEDYSPANNER_H
#define WEIGHTED_DISTANCE_ORACLE_GREEDYSPANNER_H

#include "base.h"
#include "k_skip.h"
#include "weighted_distance_oracle.h"

namespace GreedySpanner{
    using namespace std;

    float s_value = -1.0;
    int boundary_vertex_flag = -1;

    void setSeparationValue(const float &s){
        s_value = s;
    }

    void setBoundaryVertexFlag(const int &b){
        boundary_vertex_flag = b;
    }

    struct pair_hash{
        template<class T1, class T2>
        size_t operator() (const pair<T1, T2> &p) const{
            auto h1 = hash<T1>{}(p.first);
            auto h2 = hash<T2>{}(p.second);
            return h1 ^ h2;
        }
    };

    bool generateGreedySpanner(const vector<int> &boundary_vertices, kSkip::Graph &spanner, map<int, int> &new_id){
        vector<pair<int, int> > boundary_vertices_pairs;
        for (auto i = 0; i != boundary_vertices.size(); i++){
            for (auto j = i + 1; j != boundary_vertices.size(); j++){
                boundary_vertices_pairs.emplace_back(boundary_vertices[i], boundary_vertices[j]);
            }
        }
        auto seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(boundary_vertices_pairs.begin(), boundary_vertices_pairs.end(), default_random_engine(seed));
        unordered_set<pair<int, int>, pair_hash> pair_mask(boundary_vertices_pairs.begin(), boundary_vertices_pairs.end());
        vector<int> covered_s, covered_t;
        for (auto i = 0; i != boundary_vertices_pairs.size(); i++){
            if (i % 1000 == 0){
                cout << i << "/" << boundary_vertices_pairs.size() << endl;
            }
            auto current_pair = boundary_vertices_pairs[i];
//            cout << "Current pair is: " << current_pair.first << " " << current_pair.second << endl;
            if (pair_mask.find(current_pair) == pair_mask.end()){
//                cout << "Skip one pair" << endl;
                continue;
            }
            auto dijk_ret = kSkip::dijkstra(kSkip::my_base_graph, current_pair.first, current_pair.second);
            float dis = dijk_ret.first;
//            cout << "dis = " << dis << endl;
            float r = dis / (2 + 2 * s_value);
//            cout << "Radius = " << r << endl;
            kSkip::collect_covered_vertices(kSkip::my_base_graph, current_pair.first, r, covered_s, boundary_vertex_flag);
            kSkip::collect_covered_vertices(kSkip::my_base_graph, current_pair.second, r, covered_t, boundary_vertex_flag);
            for (auto u: covered_s){
                for (auto v: covered_t){
                    pair<int, int> inner_pair(min(u, v), max(u, v));
                    if (pair_mask.find(inner_pair) == pair_mask.end()){
                        continue;
                    }
                    pair_mask.erase(inner_pair);
                }
            }
            spanner.addEdge(new_id[current_pair.first], new_id[current_pair.second], dis);
            spanner.addEdge(new_id[current_pair.second], new_id[current_pair.first], dis);
        }
    }
}

#endif //WEIGHTED_DISTANCE_ORACLE_GREEDYSPANNER_H
