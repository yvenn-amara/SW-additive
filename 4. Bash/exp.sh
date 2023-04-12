# !/bin/bash

for i in 2 #1 2 5 10
do
for j in heavisine #"blocks" "bumps" "heavisine" "doppler"
do
for k in {1..10}
do
for l in 5 #1 2 5 10
do
	echo "num_weeks" $i "; lambda_num" $j "; iter" $k "; num_effects" $l
	Rscript "2. R/1_Simulated_Arrivals_v6.R" $i $j $k $l
done
done
done
done