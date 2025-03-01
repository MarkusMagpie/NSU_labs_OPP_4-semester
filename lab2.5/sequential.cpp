#include <iostream>
#include <chrono>
#include <ctime>

// разделение массива: выбираем пивот и делим массив на две части
void partition(float *v, int &i, int &j, int low, int high) {
    float pivot = v[(low + high) / 2];
    i = low;
    j = high;
    do {
        while (v[i] < pivot) i++;
        while (v[j] > pivot) j--;
        if (i <= j) {
            std::swap(v[i], v[j]);
            i++;
            j--;
        }
    } while (i <= j);
}

void quicksort_sequential(float *v, int low, int high) {
    if (low >= high) return;

    int i, j;
    partition(v, i, j, low, high);
    
    if (low < j)
        quicksort_sequential(v, low, j);
    if (i < high)
        quicksort_sequential(v, i, high);
}

bool is_sorted(float *v, int low, int high) {
    for (int i = 1; i < high; i++) {
        if (v[i-1] > v[i]) return false;
    }
    return true;
}

int main() {
    int size = 10000000;
    float* arr = new float[size];
    srand(static_cast<unsigned int>(time(0))); // устанавливаем значение системных часов в качестве стартового числа
    for (int i = 0; i < size; ++i) {
        arr[i] = static_cast<float>(rand()) / RAND_MAX * 10000.0f;
    }

    auto start = std::chrono::high_resolution_clock::now();
    quicksort_sequential(arr, 0, size - 1);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // можно вывести отсортированный массив или проверить корректность сортировки
    // if (!is_sorted(arr, 0, size)) {
    //     std::cout << "неправильно отсортировал" << std::endl;
    // }

    std::cout << elapsed.count() << std::endl;

    delete[] arr;
    return 0;
}

// g++ -o sequential sequential.cpp 
// ./sequential