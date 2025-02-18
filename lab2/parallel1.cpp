#include <iostream>
#include <vector>
#include <cstring>    // std::memcpy
#include <cmath>
#include <omp.h>      // OpenMP для параллельных вычислений
#include <fstream>
#include <chrono>

const int MAX_ITERATIONS = 100000;
const double EPSILON = 1e-2;
const double TAU = -0.01;
const int N = 2500; // размерность мтатрицы A из фаайлика Матвееыва

bool loadBinary(const std::string& filename, std::vector<float>& data, size_t expected_size) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cout << "ошибка при открытии файла: " << filename << std::endl;
        return false;
    }

    file.read(reinterpret_cast<char*>(data.data()), expected_size * sizeof(float));
    if (!file) {
        std::cerr << "ошибка чтения файла: " << filename << std::endl;
        return false;
    }

    return true;
}

void iterate(std::vector<float>& matrix_a, std::vector<float>& vector_b, std::vector<float>& vector_x, int& iterations_count) {
    float b_norm = 0;
    for (int i = 0; i < N; ++i) {
        b_norm += vector_b[i] * vector_b[i];
    }
    b_norm = std::sqrt(b_norm);
    
    float current_norm;

    for (; iterations_count < MAX_ITERATIONS; ++iterations_count) {
        current_norm = 0;

        // параллелизируем цикл по строкам матрицы
        // каждый поток обрабатывает разные строки (по индексу i) независимо
        #pragma omp parallel for reduction(+ : current_norm)
        for (int i = 0; i < N; ++i) {
            float sum = -vector_b[i];
            for (int j = 0; j < N; ++j) {
                sum += matrix_a[i * N + j] * vector_x[j];
            }
            vector_x[i] -= TAU * sum;
            current_norm += sum * sum;
        }

        float rel_error = std::sqrt(current_norm) / b_norm;
        // std::cout << "iteration: " << iterations_count << std::endl;
        if (rel_error < EPSILON) {
            // std::cout << "rel_error = " << rel_error << " ; ";
            break;
        }
    }
}

int main() {
    std::vector<float> matrix_a(N * N);
    std::vector<float> vector_b(N);
    std::vector<float> vector_x(N, 0.f); // инициализируем вектор_x нулями
    // std::vector<float> correct_x(N);

    if (!loadBinary("matA.bin", matrix_a, N * N) ||
        !loadBinary("vecB.bin", vector_b, N)) {
        return 0;
    } // ошибка при загрузке

    int iterations_count = 0;

    auto start_time = std::chrono::high_resolution_clock::now();
    iterate(matrix_a, vector_b, vector_x, iterations_count);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    // double diff = 0;
    // for (int i = 0; i < N; ++i) {
    //     diff += (vector_x[i] - correct_x[i]) * (vector_x[i] - correct_x[i]);
    // }
    // diff = std::sqrt(diff);
    // std::cout << "diff = " << diff << std::endl;

    std::cout << elapsed.count() << std::endl;

    return 0;
}

// компилируй так: g++ -fopenmp -O3 -o parallel1 parallel1.cpp