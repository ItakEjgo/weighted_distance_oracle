#/bin/bash

datasets_dir=$1;
query_dir=$2;
output_dir=$3;
algorithm=$4;
gridnum=$5;
cleaner=$6;

for file in `ls $datasets_dir`; do
    #generate query
    # ../cmake-build-release/main --generate=1 --input=$dir/$file --grid-num=$gridnum --query-num=$qnum
    
    # cp A2A.query query/$file-A2A-$flag.query
    # cp face_weight.query faceweight/$file-face_weight-$flag.query

    cp $query_dir/$file-A2A-disgap.query A2A.query
    cp $query_dir/$file-face_weight-disgap.query face_weight.query 

    ./main --input=$datasets_dir/$file --output=$output_dir/$file-$algorithm-disgap.log --method=$algorithm --grid-num=$gridnum
    python3 $cleaner $output_dir/$file-$algorithm-disgap.log 

    rm A2A.query
    rm face_weight.query
done
  
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum FixedS 5 $flag
    
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum KAlgo 1 $flag
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum UnfixedS 1 $flag
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum EAR 1 $flag
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum SE 1 default
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum MMP 5 $flag

    #fixedS 
    # ../cmake-build-release/main --input=$dir/$file --output=../results/unweighted/$file-FixedS-default.log --grid-num=16 --query-num=1000 --method=FixedS
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum FixedS 4 default

    #unfixedS 
    # ../cmake-build-release/main --input=$dir/$file --output=../results/unweighted/$file-UnfixedS-default.log --grid-num=16 --query-num=1000 --method=UnfixedS
    # bash parallel_run.sh $dir $file $gridnum 0.05 5 $qnum UnfixedS 10 eps-005
    # bash parallel_run.sh $dir $file $gridnum 0.10 5 $qnum UnfixedS 10 eps-010
    # bash parallel_run.sh $dir $file $gridnum 0.15 5 $qnum UnfixedS 10 eps-015
    # bash parallel_run.sh $dir $file $gridnum 0.20 5 $qnum UnfixedS 10 eps-020
    # bash parallel_run.sh $dir $file $gridnum 0.25 5 $qnum UnfixedS 10 eps-025

    #KAlgo
    # ../cmake-build-release/main --input=$dir/$file --output=../results/unweighted/$file-KAlgo-default.log --grid-num=16 --query-num=1000 --method=KAlgo


    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum MMP 5 $flag

    # bash parallel_run.sh $dir $file $gridnum 0.05 5 $qnum KAlgo 10 eps-005
    # bash parallel_run.sh $dir $file $gridnum 0.10 5 $qnum KAlgo 10 eps-010
    # bash parallel_run.sh $dir $file $gridnum 0.15 5 $qnum KAlgo 10 eps-015
    # bash parallel_run.sh $dir $file $gridnum 0.20 5 $qnum KAlgo 10 eps-020
    # bash parallel_run.sh $dir $file $gridnum 0.25 5 $qnum KAlgo 10 eps-025
    
    #SE
    # ../cmake-build-release/main --input=$dir/$file --output=../results/unweighted/$file-SE-default.log --grid-num=16 --query-num=1000 --method=SE
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum SE 1 default
    

    #EAR-default setting eps=0.2, m=5
    # ../cmake-build-release/main --input=$dir/$file --output=../results/unweighted/$file-EAR-default.log --grid-num=16 --query-num=1000 --method=EAR
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum EAR 1 default

    #MMP
    # ../cmake-build-release/main --input=$dir/$file --output=../results/unweighted/$file-MMP-default.log --grid-num=16 --query-num=1000 --method=MMP

# done
