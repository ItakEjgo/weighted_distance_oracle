## Efficient Arbitrary Point-to-Arbitrary Point Geodesic Distance Oracle
This repository contains the implementation of *EAR-Oracle*, *SE-Oracle*, *Fixed Scheme*, *Unfixed Scheme*, *K-Algo* and *MMP* algorithms. 
### prameters

[comment]: <> (generate arbitrary point-to-arbitrary point query or not:&#41; bool generate_flag = getarg&#40;0, "--generate")
**--generate: The flag whether generate A2A queries.** 
If this parameter set to 1, the program will be run to generate given number of queries. The queries will locate in the current working directory. The file *A2A.query* contains the generated query points and their location (the face they locate). The file *face_weight.query* contains the face weight of each face.  

[comment]: <> (string input = getarg&#40;"", "--input"&#41;,)
**--input: The path of the input file.** The input should be a *2-manifold* terrain surface in .off format. Please refer to https://en.wikipedia.org/wiki/OFF_(file_format) for more information about the .off format.

[comment]: <> (output = getarg&#40;"", "--output"&#41;;)
**--output: The path of the output file.** The output is a text file containing building time, space consumption, query time and approximate geodesic distances of the given queries.  

[comment]: <> (unsigned grid_num = getarg&#40;4, "--grid-num"&#41;;)
**--grid-num: The number of boxes of EAR-Oracle.** EAR-Oracle will generate $$\zeta^2$$ boxes for exacting highway nodes and establishing highway edges.


unsigned q_num = getarg(100, "--query-num");

bool weighted_flag = getarg(0, "--weighted");

string method_type = getarg("", "--method");
float err = getarg(0.2, "--eps");

unsigned sp_num = getarg(5, "--sp-num");

unsigned parallel_num = getarg(1, "--parallel-num");

unsigned parallel_id = getarg(0, "--parallel-id");