#!/bin/bash
# скрипт запустит программы с разным числом потоков, результаты будут записаны в results.txt

THREADS=(1 2 4 6 8 12 14 16 18 20)

# попытка равномерно распределить потоки между двумя группами ядер:
# export OMP_PLACES="{0:6}, {6:8}" # на каких группаах ядер выполнять потоки(?) с 0 по 5 и с 6 по 14
# export OMP_PROC_BIND=spread # равномерно распределять потоки по тем ядрам, которые заданы переменной OMP_PLACES

RESULTS_FILE2="results2.txt"
> "$RESULTS_FILE2"

g++ -fopenmp -O3 -o parallel2 parallel2.cpp

echo ""
echo "Одна параллельная секция"
for t in "${THREADS[@]}"; do
    echo "==========================="
    echo "Запуск с OMP_NUM_THREADS=$t"
    export OMP_NUM_THREADS=$t
    ./parallel2 | tee -a "$RESULTS_FILE2"
    echo ""
done

# очистить файл с результатами
RESULTS_FILE="results1.txt"
> "$RESULTS_FILE"

g++ -fopenmp -O3 -o parallel1 parallel1.cpp

echo "Несколько параллельных секций"
for t in "${THREADS[@]}"; do
    echo "==========================="
    echo "Запуск с OMP_NUM_THREADS=$t"
    export OMP_NUM_THREADS=$t
    # вывод в файл результатов
    ./parallel1 | tee -a "$RESULTS_FILE"
    echo ""
done

echo "Success! Результаты записаны в $RESULTS_FILE и в $RESULTS_FILE2"

echo "Запуск визуализации с помощью matplotlib"
python3 MakeGraphs.py