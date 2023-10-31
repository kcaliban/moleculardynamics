#!/bin/bash

for K in 0 150 300 450 600
do
    echo $K
    echo 8
    mpirun ./09 -f whisker_large.xyz -o whisker_large_out_${K}_10_8.xyz -x whisker_large_traj_${K}_10_8.xyz -t 5 -d 3000000 -l whisker_large_log_${K}_10_8.txt -v -s whisker_large_strain_stress_${K}_10_8.txt --tau 500 --temp $K --relax 10000 --rate 8 --dx 84 --dy 84 --dz 289 >> whisker_large_output_${K}_10_8.txt
    echo 8.5
    mpirun ./09 -f whisker_large.xyz -o whisker_large_out_${K}_10_8_5.xyz -x whisker_large_traj_${K}_10_8_5.xyz -t 5 -d 1500000 -l whisker_large_log_${K}_10_8_5.txt -v -s whisker_large_strain_stress_${K}_10_8_5.txt --tau 500 --temp $K --relax 10000 --rate 8 --factor 5 --dx 84 --dy 84 --dz 289 >> whisker_large_output_${K}_10_8_5.txt
    echo 9
    mpirun ./09 -f whisker_large.xyz -o whisker_large_out_${K}_10_9.xyz -x whisker_large_traj_${K}_10_9.xyz -t 5 -d 500000 -l whisker_large_log_${K}_10_9.txt -v -s whisker_large_strain_stress_${K}_10_9.txt --tau 500 --temp $K --relax 10000 --rate 9 --dx 84 --dy 84 --dz 289 >> whisker_large_output_${K}_10_9.txt
done
