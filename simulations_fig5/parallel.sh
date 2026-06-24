#!/bin/bash
for replica_index in 6 7 8 9 10; do
    sbatch run.sh $replica_index
    sleep 1
done
