#!/bin/bash

g++ -O3 -fopenmp parallel.cpp -o parallel

# 1 задачка: поиск минимального размера для разных порогов
echo "Задача 1: определить минимальный выгодный размера массива"
echo "size,threshold,threads,time" > results1.csv
for size in 1000 10000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000 ; do
    ./parallel $size 1000 1 >> results1.csv
done

echo "Задача 2: Зависимость времени сортировки от числа потоков"
echo "size,threshold,threads,time" > results2.csv
fixed_size=400000000 # надо большой (1e8 - примерно 10 секунд при 1 треде)
fixed_threshold=1000
for threads in 1 2 4 8 12 16; do
    ./parallel $fixed_size $fixed_threshold $threads >> results2.csv
done

echo "Success! Результаты записаны в results1.csv и в results2.csv"

echo "Запуск визуализации с помощью matplotlib"
python3 makegraphs.py