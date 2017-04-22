#!/bin/bash
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH -N 2
#SBATCH --ntasks-per-node 28

#echo commands to stdout
set -x

module unload icc

module load pgi
module load mpi/pgi_openmpi

mpicc -acc -ta=telsa:cuda8.0  linear_mpiacc.c

for var in 2 4 16 
do

	for p in 256 512 1024 2048 4096
	do
		mpirun -n $var ./a.out $p
	done

done


