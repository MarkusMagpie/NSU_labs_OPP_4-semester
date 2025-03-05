#!/bin/bash

g++ -O3 -fopenmp parallel.cpp -o parallel

# 1 - минимальный размер для разных порогов
echo "1 - определить минимальный выгодный размера массива"
echo "size,threshold,threads,time" > results1.csv
threshold=1000
threads=1
for size in 1000 10000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000 ; do
    ./parallel $size $threshold $threads >> results1.csv
done

echo "2 - показать зависимость времени сортировки от числа потоков"
echo "size,threshold,threads,time" > results2.csv
size=100000000 # надо большой (1e8 - примерно 10 секунд при 1 треде)
threshold=1000
for threads in 1 2 4 8 10 12 14 16; do
    ./parallel $size $threshold $threads >> results2.csv
done

echo "3 - определить оптимальный порог (threshold)"
echo "size,threshold,threads,time" > results3.csv
fixed_size=100000000
threads=6
for threshold in 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000 200000000 500000000; do
    ./parallel $fixed_size $threshold $threads >> results3.csv
done

echo "Success! Результаты записаны в results1.csv, results2.csv и results3.csv"
echo "Запуск визуализации с помощью matplotlib"
python3 makegraphs.py