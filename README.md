## EAR-Oracle: On Efficient Indexing for Distance Queries between Arbitrary Points on Terrain Surface
This repository contains the implementation of *EAR-Oracle* published in SIGMOD'2023 (also includes implementation of *SE-Oracle*, *Fixed Scheme*, *Unfixed Scheme*, *K-Algo* and *MMP* (by CGAL) algorithms.). 
***

### Structure Overview (and important files)
        ./
        |-build/ # folder of compile and run the master script.
        | |-out/ # final figures in .eps format.
        | |-master_script.py # a master script for running all experiments and ploting.
        | |-script_config # several paths and configs of the master script.
        |-datasets # all datasets used.
        | |-query # querys for all experiments.
        |-exp # folder stores the final experimental results.
        |-scripts # scripts for experiments, ploting and additional .off datasets.
        |-technicalreport # full version of our paper.
        |-*.h/hpp/cpp # source code.

### Environment & Run ( tested in Ubuntu 18.04LTS and Debian GNU/Linux 11 )
0. Install the following packages: *build-essential, cmake, libcgal-dev, gnuplot, tmux* [ **Note:** *tmux* is for running long-time experiments, any other substitute is fine. Please make sure the CGAL(*libcgal-dev*) version is newer than v5.2 (for *cmake* support).]:

        apt-get update
        apt-get upgrade
        apt-get install -y build-essential cmake libcgal-dev gnuplot tmux

1. Clone this repo and go to build folder:

        git clone https://github.com/ItakEjgo/weighted_distance_oracle/
        cd weighted_distance_oracle/build

2. Compile by *cmake*:

        rm CMakeCache.txt
        cmake .. -DCMAKE_BUILD_TYPE=Release # CGAL reuqires to be built in release mode.
        make

3. Check the *script_config* and run *master script*:

        tmux
        less script_config # check config (normally nothing to change).
        python3 master_script.py script_config > run.log # store the log.
        Ctrl + B, D # leave this tmux session.
        #You can use |tmux a -t 0| to touch the previous session at any time later.

4. The experimental results will be generated in the *../exp* folder. The final figures will locate in *out/* folder.

### Note for Reproducing
0. All experimental results (except for scalability test) is tested on a machine with *256GB* memory. To support verification on machines with smaller memory space, the **default config** (locate in *build/script_config*) is only set for small datasets (i.e., four terrains in *datasets/small*). Please change **dataset_list(line5)** and **tested_dataset(line10)** in the **default config** if you want to test all datasets.
1. For large datasets, the tested algorithms requires a large memory space. You may need to use *ulimit* command to adjust the memory limit.
2. The *UnfixedS* and *KAlgo* are on-the-fly algorithms, therefore they are **very slow** for median and large datasets (since they need to place a lot of auxiliary points and distance quires are independent to each other). If you can access multiple machines, we suggest separating the experiments into multiple machines (For example, each machine runs a single experiment. You can achieve this by *comment/uncomment **line 650-670*** of the master script).
**Fully** run the on-the-fly algorithms on large datasets may **require a couple of weeks**. We suggest to add a time limit for these on-the-fly algorithms(UnfixedS, KAlgo). You can achieve this by modifying the *exp_\*.sh* in folder *scripts/exp/* (such as using *timelimit* command).
3. The output figures are in *.eps format. To check them, you can use **Epstool** on Windows (Please remember to adjust the portait size (2000 * 1000 is large enough)). Or you can simply create a *.tex file in overleaf, insert the *.eps file and compile (view) it:

        \begin{figure*}
            \centering
            \includegraphics[width=\linewidth]{default-2row.eps}
            \caption{Test Fig}
            \label{fig:test_fig}
        \end{figure*}

4. All experiments are affected by the query points generated. Therefore **randomness** influences the final figures to some extent, especially for the *query time* and *relative error*. We set the y-axis of relative error figure in the range [0.00, 0.06]. In case of the relative error curve does not appeared in the figure completely, you can change the y-axis scale by modifying the following code in the corresponding plot file (in *scripts/figures/\*.plot*):

        set yrange [0.000 : 0.06001]
5. Generate the query for distance gap is quite time consuming (we need to calculte many distances and repeat if they are out of a certain range). Therefore, we provide distance queries in the *datasets/query* folder. For scalability test, it becomes **extremly slow** for on-the-fly algorithms. Therefore we reduce the number of queries to 15. 
### Scripts
#### Dataset
We provided the datasets used in the paper in **./datasets** folder. The datasets used in *EAR-Oracle* must be *2-manifold* terrain surfaces in *.off* format.

In addition, We also provide a script *dem2off.py* in **./scripts/DEM2OFF** folder for processing digital elevation models (DEM) to 2-manifold terrain surfaces in .off format. [**Note:** additional *gdalwarp* package must be installed.]:

(1) Choose a DEM from *OpenTopography* (https://portal.opentopography.org/datasets). To choose a DEM, OpenTopography provides convenient interfaces that the DEM from a region could be selected by selecting a rectangle on the terrain.

(2) Download the GeoTiff file (in .tif format) and run *dem2off.py* for the .tif file. The script could be used to simplify the DEM to a .off file with a certain number of faces.
#### Experiments
The scripts of the experiments are located in **./scripts/exp** folder. 
        
        clean.py # clean the results (*.log files) to *.cln files.
        generate_query.sh # generate distance queries.
        exp_*.sh # scripts for all experiments.
These scripts are called by the *master_script* (build/master_scripts.py). Normally doing experiments just need to change the *config* (build/script_config) file.

#### plots
The scripts of plotting are located in **./scripts/figures** folder.
We use **gnuplot** to plot the figures in the paper. They can be called by the *master script*. But you can also run them individually (once you have the data):

        gnuplot -c [PLOT_SCRIPT_PATH] [DATA_PATH]
The data used for plotting is generated by the *master script* and locate in *figures/data* folder.

### Running Prameter
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
### Running Example
##### Generate 1000 queries for *example.off* dataset.
--generate=1 --input=../datasets/example.off --output=../results/test.log --grid-num=16 --query-num=1000

##### Run *EAR-Oracle* for query processing
--generate=0 --input=../datasets/example.off --output=../results/test.log --grid-num=16 --query-num=1000 --method=EAR

##### Run *Fixed Scheme* for query processing
--generate=0 --input=../datasets/example.off --output=../results/test.log --grid-num=16 --query-num=1000 --method=FixedS
***
###Reference
If you use our code, please cite our paper:

        @article{10.1145/3588694,
        author = {Huang, Bo and Wei, Victor Junqiu and Wong, Raymond Chi-Wing and Tang, Bo},
        title = {EAR-Oracle: On Efficient Indexing for Distance Queries between Arbitrary Points on Terrain Surface},
        year = {2023},
        issue_date = {May 2023},
        publisher = {Association for Computing Machinery},
        address = {New York, NY, USA},
        volume = {1},
        number = {1},
        url = {https://doi.org/10.1145/3588694},
        doi = {10.1145/3588694},
        journal = {Proc. ACM Manag. Data},
        month = {may},
        articleno = {14},
        numpages = {26},
        }
    