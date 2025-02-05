#!/bin/bash
# скрипт запустит программы с разным числом потоков, результаты будут записаны в results.txt

g++ -fopenmp -o parallel parallel1.cpp

# очистить файл с результатами
RESULTS_FILE="results.txt"
> "$RESULTS_FILE"

THREADS=(1 2 4 8)

for t in "${THREADS[@]}"; do
    echo "========================================"
    echo "Запуск с OMP_NUM_THREADS=$t"
    export OMP_NUM_THREADS=$t
    # вывод в файл результатов
    ./parallel | tee -a "$RESULTS_FILE"
    echo ""
done

echo "Success! Результаты записаны в $RESULTS_FILE"
