#!/bin/bash

cd "${0%/*}" # set working directory to directory of this script

sz=3871

# first step: ignore velocities, add no heat
./07 -f cluster_${sz}_0.xyz -o cluster_${sz}_out_1.xyz -x cluster_${sz}_out_0_traj.xyz -l cluster_${sz}_out_0_log.txt -t 10 -r 10000 -m 10000 --et cluster_${sz}_out_et.txt -v

for i in {2..600}
do
    echo "Step $i"
    previous="$((i-1))"
    ./07 -f cluster_${sz}_out_$previous.xyz -o cluster_${sz}_out_$i.xyz -x cluster_${sz}_out_${i}_traj.xyz -l cluster_${sz}_out_${i}_log.txt -t 10 -q -k 1 -r 10000 -m 10000 --et cluster_${sz}_out_et.txt
done
