# weighted_distance_oracle
weighted distance oracle implementation

Usage:

cd build/
../cmake-build-release/main -g 0 -m 'path-to-file' -w 0 -e 0.2 -l 5 -q 1000 -a 0 -s 5 -o 'output-name'
-g: generate queries and face weight or not.
-m: path to the mesh file (.off format).
-w: is the terrain weighted or not? 0: unweighted. 1: weighted.
-e: the value of error-bound epsilon.
-l: quadtree level.
-q: number of queries. (if g = 1, we will generate this number of queries)
-a: algorithm type. Please see annotation of function run in main.cpp for details.
-s: the number of Steiner points for each edge. This parameter does not influence Unfixed Scheme.
-o: output file name.

