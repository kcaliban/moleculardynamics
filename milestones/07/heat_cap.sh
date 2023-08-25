#!/bin/bash

for i in {1..50}
do
    echo "Step $i"
    previous="$((i-1))"
    ./07 -f cluster_923_5000_$previous.xyz -o cluster_923_5000_$i.xyz -x cluster_923_5000_${i}_traj.xyz -l cluster_923_5000_${i}_log.txt -t 1 -q -k 75 -r 4000 -m 1000
done
