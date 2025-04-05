#!/bin/bash

#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -m n

cd $PBS_O_WORKDIR
## количество строк в файле $PBS_NODEFILE, то есть количество узлов, на которых будет выполняться задача
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
echo "Number of MPI processes: $MPI_NP"
echo "File $PBS_NODEFILE:"
cat $PBS_NODEFILE
echo "==========================="

/mnt/storage/home/hpcusers/hpcuser09/mpe/bin/mpecc -mpilog parallel.cpp -o parallel
mpirun -machinefile $PBS_NODEFILE -np $MPI_NP ./parallel

echo "Программа завершена. Проверьте файл parallel.clog2 для логов MPE."
