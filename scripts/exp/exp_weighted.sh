#/bin/bash

datasets_dir=$1;
query_dir=$2;
output_dir=$3;
algorithm=$4;
gridnum=$5;
cleaner=$6;


for file in `ls $datasets_dir`; do
    cp $query_dir/$file-A2A-weighted.query A2A.query
    cp $query_dir/$file-face_weight-weighted.query face_weight.query 

    ./main --weighted=1 --input=$datasets_dir/$file --output=$output_dir/$file-$algorithm-weighted.log --method=$algorithm --grid-num=$gridnum
    python3 $cleaner $output_dir/$file-$algorithm-weighted.log 

    rm A2A.query
    rm face_weight.query
done
