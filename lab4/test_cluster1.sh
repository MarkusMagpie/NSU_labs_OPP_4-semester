#!/bin/bash

#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -m n

cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
echo "=== Running with $MPI_NP processes ==="

mpicxx -O3 parallel.cpp -o parallel
mpirun -machinefile $PBS_NODEFILE -np $MPI_NP ./parallel