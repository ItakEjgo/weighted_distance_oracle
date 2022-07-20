#/bin/bash

dir=$1;
gridnum=$2;
qnum=$3;
flag=$4;
processnum=$5;


for file in `ls $dir`; do
    #generate query
    # ../cmake-build-release/main --weighted=1 --generate=1 --input=$dir/$file --grid-num=$gridnum --query-num=$qnum
    
    # cp A2A.query query/$file-A2A-$flag.query
    # cp face_weight.query faceweight/$file-face_weight-$flag.query

    cp query/$file-A2A-$flag.query A2A.query
    cp faceweight/$file-face_weight-$flag.query face_weight.query 

    bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum FixedS 1 $flag
    
    bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum KAlgo 1 $flag
    bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum UnfixedS 1 $flag
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum EAR 1 $flag
    bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum SE 1 $flag
    # bash parallel_run.sh $dir $file $gridnum 0.2 5 $qnum MMP 1 $flag 

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

done
