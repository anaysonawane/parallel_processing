#!/bin/bash

MAX_NUM_THREADS=28
NUM_THREADS=0

while [ $NUM_THREADS -lt $MAX_NUM_THREADS ]
do

NUM_THREADS=`expr $NUM_THREADS + 1`

sbatch --ntasks-per-node $NUM_THREADS job_openmpi.sh

done

