#!/bin/bash

sz=3871

for i in {1..600}
do
    echo "Step $i"
    previous="$((i-1))"
    ./07 -f cluster_${sz}_out_$previous.xyz -o cluster_${sz}_out_$i.xyz -x cluster_${sz}_out_${i}_traj.xyz -l cluster_${sz}_out_${i}_log.txt -t 10 -q -k 1 -r 10000 -m 10000 --et cluster_${sz}_out_et.txt
done
