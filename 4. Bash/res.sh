# !/bin/bash

for i in 2 #1
do
for j in heavisine #"blocks" "bumps" "heavisine" "doppler"
do
for k in 5 #1 2 5 10
do
   Rscript "2. R/2_Results.r" $i $j $k
done
done
done