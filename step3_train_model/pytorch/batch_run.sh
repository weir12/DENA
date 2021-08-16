#!/bin/bash

motifs=("AAACA" "AAACT" "AGACC" "GAACA" "GAACT" "GGACC" "AAACC" "AGACA" "AGACT" "GAACC" "GGACA" "GGACT")

for motif in ${motifs[@]}
do
python ./train.py -c config.json -i 2_2${motif} --motif ${motif}
done