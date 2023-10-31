#!/bin/bash

for K in 0 150 300 450 600
do
    echo $K
    echo 8
    mpirun -n 4 ./09 -f whisker_small.xyz -o whisker_small_out_${K}_10_8.xyz -x whisker_small_traj_${K}_10_8.xyz -t 5 -d 3000000 -l whisker_small_log_${K}_10_8.txt -v -s whisker_small_strain_stress_${K}_10_8.txt --tau 500 --temp $K --relax 10000 --rate 8 >> whisker_small_output_${K}_10_8.txt
    echo 8.5
    mpirun -n 4 ./09 -f whisker_small.xyz -o whisker_small_out_${K}_10_8_5.xyz -x whisker_small_traj_${K}_10_8_5.xyz -t 5 -d 1500000 -l whisker_small_log_${K}_10_8_5.txt -v -s whisker_small_strain_stress_${K}_10_8_5.txt --tau 500 --temp $K --relax 10000 --rate 8 --factor 5 >> whisker_small_output_${K}_10_8_5.txt
    echo 9
    mpirun -n 4 ./09 -f whisker_small.xyz -o whisker_small_out_${K}_10_9.xyz -x whisker_small_traj_${K}_10_9.xyz -t 5 -d 500000 -l whisker_small_log_${K}_10_9.txt -v -s whisker_small_strain_stress_${K}_10_9.txt --tau 500 --temp $K --relax 10000 --rate 9 >> whisker_small_output_${K}_10_9.txt
done