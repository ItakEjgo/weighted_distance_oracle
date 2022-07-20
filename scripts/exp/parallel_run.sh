#/bin/bash

dir=$1;
file=$2;
gridnum=$3;
eps=$4
spnum=$5
qnum=$6;
method=$7;
parallelnum=$8;
flag=$9;


#generate query
# ../cmake-build-release/main --generate=1 --input=$dir/$file --grid_num=16 --query-num=1000

# cp A2A.query $file-A2A.query
# cp face_weight.query $file-face_weight.query

# cp $file-A2A.query A2A.query
# cp $file-face_weight.query face_weight.query

# ../cmake-build-release/main --input=$dir/$file --output=../results/unweighted/test-$file-$method.log --grid_num=16 --query-num=1000 --method=$method --parallel-num=1 --parallel-id=0 &

#fixedS 
for ((i=0; i<$parallelnum; i++)) do
    # ../cmake-build-release/main --weighted=1 --input=$dir/$file --output=../exp/$file-$method-part$i.log --grid-num=$gridnum --eps=$eps --sp-num=$spnum --query-num=$qnum --method=$method --parallel-num=$parallelnum --parallel-id=$i &
    ../cmake-build-release/main --input=$dir/$file --output=../exp/$file-$method-part$i.log --grid-num=$gridnum --eps=$eps --sp-num=$spnum --query-num=$qnum --method=$method --parallel-num=$parallelnum --parallel-id=$i 
done
wait
echo "parallel run $dir/$file $method finished."  
python3 deal.py $file-$method $parallelnum $flag $qnum


