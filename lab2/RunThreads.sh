#!/bin/bash
# скрипт запустит программы с разным числом потоков, результаты будут записаны в results.txt

g++ -fopenmp -Ofast -o parallel1 parallel1.cpp

# очистить файл с результатами
RESULTS_FILE="results1.txt"
> "$RESULTS_FILE"

THREADS=(1 2 4 8)

echo "Несколько параллельных секций"
for t in "${THREADS[@]}"; do
    echo "========================================"
    echo "Запуск с OMP_NUM_THREADS=$t"
    export OMP_NUM_THREADS=$t
    # вывод в файл результатов
    ./parallel1 | tee -a "$RESULTS_FILE"
    echo ""
done

# g++ -fopenmp -o parallel2 parallel2.cpp

# RESULTS_FILE2="results2.txt"
# > "$RESULTS_FILE2"

# echo ""
# echo "Одна параллельная секция"
# for t in "${THREADS[@]}"; do
#     echo "========================================"
#     echo "Запуск с OMP_NUM_THREADS=$t"
#     export OMP_NUM_THREADS=$t
#     ./parallel2 | tee -a "$RESULTS_FILE2"
#     echo ""
# done

# echo "Success! Результаты записаны в $RESULTS_FILE и в $RESULTS_FILE2"

echo "Запуск визуализации с помощью matplotlib"
python3 MakeGraphs.py