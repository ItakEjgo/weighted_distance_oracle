## Efficient Arbitrary Point-to-Arbitrary Point Geodesic Distance Oracle
This repository contains the implementation of *EAR-Oracle*, *SE-Oracle*, *Fixed Scheme*, *Unfixed Scheme*, *K-Algo* and *MMP* algorithms.

### Dataset
We provided the datasets used in the paper in *datasets/* folder. The datasets used in EAR-Oracle should be *2-manifold* terrain surfaces in *.off* format. 

In addition, We also provide a script for collecting real world digital elevation models (DEM) to 2-manifold terrain surfaces in .off format.
Firstly, choose a DEM from *OpenTopography* (https://portal.opentopography.org/datasets). To choose a DEM, OpenTopography provides convenient interfaces that the DEM from a region could be selected by selecting a rectangle on the terrain.
Secondly, download the GeoTiff file (in .tif) format and run *dem2off.py* for the .tif file in DEM2OFF folder. The script could be used to simplify the DEM to a off file with given number of faces.


### Prameter
[comment]: <> (generate arbitrary point-to-arbitrary point query or not:&#41; bool generate_flag = getarg&#40;0, "--generate")
**--generate={0, 1}: The flag whether generate A2A queries.** 
If this parameter set to 1, the program will be run to generate given number of queries. The queries will locate in the current working directory. The file *A2A.query* contains the generated query points and their location (the face they locate). The file *face_weight.query* contains the face weight of each face.  

[comment]: <> (string input = getarg&#40;"", "--input"&#41;,)
**--input={input file path}: The path of the input file.** The input should be a *2-manifold* terrain surface in .off format. Please refer to https://en.wikipedia.org/wiki/OFF_(file_format) for more information about the .off format.

[comment]: <> (output = getarg&#40;"", "--output"&#41;;)
**--output={output file path}: The path of the output file.** The output is a text file containing building time, space consumption, query time and approximate geodesic distances of the given queries.  

[comment]: <> (unsigned grid_num = getarg&#40;4, "--grid-num"&#41;;)
**--grid-num={int}: The number of boxes of EAR-Oracle.** EAR-Oracle will generate $\zeta^2$ boxes for constructing *highway network*. Please note that this parameter is the square of $\zeta$.

[comment]: <> (unsigned q_num = getarg&#40;100, "--query-num"&#41;;)
**--query-num={int}: The number of queries.** If the --generate flag set to be 1, the program will generate 3 * query-num queries (1/3 of them are inner-box queries, 1/3 of them are inter-box queries, and the rest are inner-inter mixed queries). 

[comment]: <> (bool weighted_flag = getarg&#40;0, "--weighted"&#41;;)
**--weighted={0, 1}: The flag whether the terrain surface is weighted or not.** If this parameter set to be 1, the terrain surface (input file) is a weighted terrain surface and thus, the face weights will be generated and be used for distance computation. Please note that only *Fixed Scheme*, *Unfixed Scheme* and *EAR-Oracle* support weighted terrain surfaces queries.  

[comment]: <> (string method_type = getarg&#40;"", "--method"&#41;;)
**--method={'FixedS', 'UnfixedS', 'KAlgo', 'SE', 'EAR', 'MMP'}: Which algorithm is going to be used.** *Fixed Scheme*, *Unfixed Scheme* and *K-Algo* are state-of-the-art on-the-fly algorithms. *SE-Oracle* is state-of-the-art index-based algorithm. *EAR-Oracle* is the proposed *highway network-based indexing algorithm*. *MMP* is the exact geodesic distance algorithm and it is implemented by the CGAL(v5.3) library. 

[comment]: <> (float err = getarg&#40;0.2, "--eps"&#41;;)
**--eps={float}: The value of multiplicative error bound.** 

[comment]: <> (unsigned sp_num = getarg&#40;5, "--sp-num"&#41;;)
**--sp-num={int}: The number of Steiner points placed on each bisector.** This parameter will be used for the bisector-fixed scheme base graph construction. 

[comment]: <> (unsigned parallel_num = getarg&#40;1, "--parallel-num"&#41;;)
**--parallel-num={int}: The number of processes to used.** Please note that this parameter only works for on-the-fly algorithms (*Fixed Scheme, Unfixed Scheme, K-Algo and MMP*). Besides, the query-num should be a multiple of query-num.  

[comment]: <> (unsigned parallel_id = getarg&#40;0, "--parallel-id"&#41;;)
**--parallel-id={int}: Which partition of queries this process is going to deal.** Please note that the parallel-id should be less than the parallel-num.

###Example
##### Generate 1000 queries for *small_terrain_ori.off* dataset.
--generate=1 --input=../datasets/small_terrain_ori.off --output=../results/test.log --grid-num=16 --query-num=1000

##### Run *EAR-Oracle* for query processing
--generate=0 --input=../datasets/small_terrain_ori.off --output=../results/test.log --grid-num=16 --query-num=1000 --method=EAR

##### Run *Fixed Scheme* for query processing
--generate=0 --input=../datasets/small_terrain_ori.off --output=../results/test.log --grid-num=16 --query-num=1000 --method=FixedS