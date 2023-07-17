#/bin/bash

datasets_dir=$1;
query_dir=$2;
gridnum=$3;
query_num=$4;
type=$5;

for file in `ls $datasets_dir`; do
    #generate query
    ./main --generate=1 --input=$datasets_dir/$file --grid-num=$gridnum --query-num=$query_num
    mv A2A.query $query_dir/$file-A2A-$type.query
    mv face_weight.query $query_dir/$file-face_weight-$type.query

    ./main --generate=1 --input=$datasets_dir/$file --weighted=1 --grid-num=$gridnum --query-num=$query_num 
    mv A2A.query $query_dir/$file-A2A-$type.query
    mv face_weight.query $query_dir/$file-face_weight-$type.query
    
    # cp A2A.query query/$file-A2A-$flag.query
    # cp face_weight.query faceweight/$file-face_weight-$flag.query

    # cp $query_dir/$file-A2A-default.query A2A.query
    # cp $query_dir/$file-face_weight-default.query face_weight.query 

    # ./main --input=$datasets_dir/$file --output=$output_dir/$file-EAR-default.log --grid-num=$gridnum --query-num=$query_num --method=$algorithm
    # python3 $cleaner $output_dir/$file-EAR-default.log $query_num

    # rm A2A.query
    # rm face_weight.query
done
