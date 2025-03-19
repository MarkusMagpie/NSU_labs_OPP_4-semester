#!/bin/bash

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
# echo "Запуск визуализации с помощью matplotlib:"
# python3 makegraphs.py
