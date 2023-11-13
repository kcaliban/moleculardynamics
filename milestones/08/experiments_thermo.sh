cd "${0%/*}" # set working directory to directory of this script

./08 -f cluster_5083_transformed.xyz -o cluster_5083_transformed_out_1.xyz -t 10 -d 100000 -l cluster_5083_transformed_log_1.txt -v --thermostat --t0 300 --tau 1000 
mpirun -n 2 ./08 -f cluster_5083_transformed.xyz -o cluster_5083_transformed_out_2.xyz -t 10 -d 100000 -l cluster_5083_transformed_log_2.txt -v --thermostat --t0 300 --tau 1000 
mpirun -n 4 ./08 -f cluster_5083_transformed.xyz -o cluster_5083_transformed_out_4.xyz -t 10 -d 100000 -l cluster_5083_transformed_log_4.txt -v --thermostat --t0 300 --tau 1000 
mpirun -n 8 ./08 -f cluster_5083_transformed.xyz -o cluster_5083_transformed_out_8.xyz -t 10 -d 100000 -l cluster_5083_transformed_log_8.txt -v --thermostat --t0 300 --tau 1000 