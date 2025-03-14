#!/bin/bash

#PBS -N test_parallel
#PBS -l walltime=00:10:00      
#PBS -l select=1:ncpus=16:mpiprocs=16:mem=8000m

cd $PBS_O_WORKDIR

# настройка окружения: загрузка переменных для ITAC и Intel компилятора
source /opt/intel/itac/8.1.3.037/bin/itacvars.sh
source /opt/intel/composerxe/bin/compilervars.sh intel64
export I_MPI_CC=icc

NP_VALUES=(1 2 4 6 8 10 12 14 16)

mpic++ -O3 -o parallel1 parallel1.cpp

RESULTS_MPI="results_mpi.txt"
> "$RESULTS_MPI"

echo "Запуск параллельного варианта с различным числом процессов"
for np in "${NP_VALUES[@]}"; do
    echo "==========================="
    echo "mpirun -np $np ./parallel1"
    mpirun -np $np ./parallel1 | tee -a "$RESULTS_MPI"
    echo ""
done

echo "Success! Результаты замеров записаны в файл $RESULTS_MPI"
echo "Запуск визуализации с помощью matplotlib:"
python3 makegraphs.py
