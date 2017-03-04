#!/bin/bash
##File name: gpu_start.sh
##Author: Anay Sonawane
#SBATCH -N 2
#SBATCH -p GPU
#SBATCH --ntasks-per-node 28
#SBATCH -t 5:00:00
#SBATCH --gres=gpu:4

#echo commands to stdout
set -x

pgcc -acc -Minfo=accel -ta=tesla:cuda8.0 floyd_openacc.c

echo 100 | ./a.out
