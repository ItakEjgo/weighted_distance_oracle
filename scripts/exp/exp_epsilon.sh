datasets_dir=$1;
query_dir=$2;
output_dir=$3;
algorithm=$4;
eps_val=$5;
flag=$6;
cleaner=$7;

for file in `ls $datasets_dir`; do
    cp $query_dir/$file-A2A-default.query A2A.query
    cp $query_dir/$file-face_weight-default.query face_weight.query 

    ./main --input=$datasets_dir/$file --output=$output_dir/$file-$algorithm-$flag.log --method=$algorithm --eps=$eps_val 
    python3 $cleaner $output_dir/$file-$algorithm-$flag.log 

    rm A2A.query
    rm face_weight.query
done