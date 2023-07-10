#!/bin/bash

# Measure time for different simulation runs
for i in {2..8}
do
    output=$(TIMEFORMAT="%R"; { time ./05 $i; } 2>&1)
    echo -e "$i\t$output"
done
