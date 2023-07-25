#/bin/bash

datasets_dir=$1;
query_dir=$2;
weight_flag=$3;
gridnum=$4;
type=$5;
disgap=$6

for file in `ls $datasets_dir`; do
    #generate query
    ./main --generate=1 --weighted=$weight_flag --input=$datasets_dir/$file --grid-num=$gridnum --distance-gap=$disgap
    mv A2A.query $query_dir/$file-A2A-$type.query
    mv face_weight.query $query_dir/$file-face_weight-$type.query
done
