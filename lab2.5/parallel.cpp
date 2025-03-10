#include <iostream>
#include <ctime>
#include <vector>
#include <chrono>

#include <omp.h>

// const int threshold = 1000; // ниже этого размера задачи выполняются последовательно

struct Params {
    int size;
    int threshold;
    int num_threads;
};

void partition(float *v, int low, int high, int &i, int &j) {
    float pivot = v[(low + high) / 2];
    i = low;
    j = high;
    while (i <= j) {
        while (v[i] < pivot) i++;
        while (v[j] > pivot) j--;
        if (i <= j) {
            std::swap(v[i], v[j]);
            i++;
            j--;
        }
    }
}

void quicksort_parallel(float *v, int low, int high, int threshold) {
    if (low >= high) return;
    
    int i, j;
    partition(v, low, high, i, j);
    
    // массив и подмассивы разделенные пивотом должны быть >= threshhold ВСЕ чтобы использовать задачи
    #pragma omp task shared(v) firstprivate(low, j, threshold) if((j - low) >= threshold)
    {
        quicksort_parallel(v, low, j, threshold);
    }
    
    #pragma omp task shared(v) firstprivate(i, high, threshold) if((high - i) >= threshold)
    {
        quicksort_parallel(v, i, high, threshold);
    }

    // #pragma omp taskwait
}

bool is_sorted(float *v, int low, int high) {
    for (int i = 1; i < high; i++) {
        if (v[i-1] > v[i]) return false;
    }
    return true;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cout << "Неверное число параметров" << std::endl;
        return 0;
    }
    
    int runs = 5;
    float best_time = std::numeric_limits<float>::max();

    Params params;
    params.size = std::stoi(argv[1]);
    params.threshold = std::stoi(argv[2]);
    params.num_threads = std::stoi(argv[3]);

    omp_set_num_threads(params.num_threads);

    std::vector<float> arr(params.size);
    srand(static_cast<unsigned int>(time(0))); // устанавливаем значение системных часов в качестве стартового числа
    for (int i = 0; i < params.size; ++i) {
        arr[i] = static_cast<float>(rand()) / RAND_MAX * 10000.0f;
    }
    
    for (int i = 0; i < runs; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        // параллельная область, внутри которой единственный поток начинает сортировку
        #pragma omp parallel
        {
            #pragma omp single
            {
                quicksort_parallel(arr.data(), 0, params.size - 1, params.threshold);
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        if (elapsed.count() < best_time) {
            best_time = elapsed.count();
        }
    }

    if (!is_sorted(arr.data(), 0, params.size)) {
        std::cout << "неправильно отсортировал" << std::endl;
    }

    std::cout << params.size << "," << params.threshold << "," << params.num_threads << "," << best_time << "\n";
    
    return 0;
}

// g++ -o parallel parallel.cpp -fopenmp
// ./parallel