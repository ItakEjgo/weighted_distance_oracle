#/bin/bash

bash exp_weighted.sh ../datasets/final_dataset/small 16 100 weighted 1
bash exp_weighted.sh ../datasets/final_dataset/median 256 100 weighted 1
bash exp_weighted.sh ../datasets/final_dataset/large 256 100 weighted 1
# bash exp_default.sh ../datasets/bk 16 5 default