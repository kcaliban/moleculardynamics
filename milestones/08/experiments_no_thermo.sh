./08 -f cluster_5083_transformed.xyz -o cluster_5083_transformed_out_1.xyz -t 10 -d 100000 -l cluster_5083_transformed_log_1.txt -v 
mpirun -n 2 ./08 -f cluster_5083_transformed.xyz -o cluster_5083_transformed_out_2.xyz -t 10 -d 100000 -l cluster_5083_transformed_log_2.txt -v
mpirun -n 4 ./08 -f cluster_5083_transformed.xyz -o cluster_5083_transformed_out_4.xyz -t 10 -d 100000 -l cluster_5083_transformed_log_4.txt -v
mpirun -n 8 ./08 -f cluster_5083_transformed.xyz -o cluster_5083_transformed_out_8.xyz -t 10 -d 100000 -l cluster_5083_transformed_log_8.txt -v