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
    double b_norm = 0.0;
    for (int i = 0; i < N; ++i) {
        b_norm += vector_b[i] * vector_b[i];
    }
    b_norm = std::sqrt(b_norm);

    bool stop = false;
    float current_norm, local_norm;

    #pragma omp parallel 
    {
        while (!stop) {
            local_norm = 0.f;
            #pragma omp for schedule(static) reduction(+ : local_norm)
            for (int i = 0; i < N; ++i) {
                float sum = -vector_b[i];
                for (int j = 0; j < N; ++j) {
                    sum += matrix_a[i * N + j] * vector_x[j];
                }
                vector_x[i] -= TAU * sum;
                local_norm += sum * sum;
            }

            // только master тред обновляет переменные и проверяет условие выхода; если оно выполняется, то
            // все потоки выйдут из цикла в следующей итерации
            #pragma omp master
            {
                current_norm = local_norm;
                ++iterations_count;
                // std::cout << "Iteration " << iterations_count << ": " << rel_error << std::endl;

                double rel_error = std::sqrt(current_norm) / b_norm;
                if (rel_error < EPSILON) {
                    stop = true;
                }
            }
            #pragma omp barrier // синхронизация потоков между итерациями - все ждут выполнения одной итерации и переходят к следующей
        }
    }
}

int main() {
    std::vector<float> matrix_a(N * N);
    std::vector<float> vector_b(N);
    std::vector<float> vector_x(N, 0.f); // инициализируем вектор_x нулями

    if (!loadBinary("matA.bin", matrix_a, N * N) ||
        !loadBinary("vecB.bin", vector_b, N)) {
        return 0;
    }

    int iterations_count = 0;

    auto start_time = std::chrono::high_resolution_clock::now();
    iterate(matrix_a, vector_b, vector_x, iterations_count);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << elapsed.count() << std::endl;
    // std::cout << "Iterations: " << iterations_count << std::endl;

    return 0;
}

// компилируй так: g++ -fopenmp -O1 -o parallel2 parallel2.cpp