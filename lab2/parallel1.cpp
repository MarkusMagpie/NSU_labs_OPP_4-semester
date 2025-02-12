#include <iostream>
#include <vector>
#include <cstring>    // std::memcpy
#include <cmath>
#include <omp.h>      // OpenMP для параллельных вычислений
#include <fstream>

const int MAX_ITERATIONS = 100000;
const double EPSILON = 1e-3;
const double TAU = -0.03;
const int N = 2500; // размерность мтатрицы A из фаайлика Матвееыва

// void initialize(std::vector<float>& matrix, std::vector<double>& vector_b, std::vector<double>& vector_x) {
//     for (int i = 0; i < N; ++i) {
//         vector_x[i] = 0;
//         vector_b[i] = N + 1;
//         for (int j = 0; j < N; ++j) {
//             if (i == j) {
//                 matrix[i * N + j] = 2.0;
//             } else {
//                 matrix[i * N + j] = 1.0;
//             }
//         }
//     }
// }

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

float multMatrixVector(const std::vector<float>& matrix, const std::vector<float>& vector) {
    float result = 0;
    for (int i = 0; i < N; ++i) {
        result += matrix[i] * vector[i];
    }
    return result;
}

int main() {
    // функция из OpenMP для подсчета времени
    double start_time = omp_get_wtime();

    std::vector<float> matrix_a(N * N);
    std::vector<float> vector_b(N);
    std::vector<float> vector_x(N);

    // std::vector<float> correct_x(N);

    if (!loadBinary("matA.bin", matrix_a, N * N) ||
        !loadBinary("vecB.bin", vector_b, N)) {
        return 0;
    } // типа ошибка при загрузке
    
    std::fill(vector_x.begin(), vector_x.end(), 0.f);

    int iterations_count = 0;
    float current_norm = 0;

    float b_norm = 0;
    for (int i = 0; i < N; ++i) {
        b_norm += vector_b[i] * vector_b[i];
    }
    b_norm = std::sqrt(b_norm);

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
        // std::cout << "Iteration " << iterations_count << ": " << rel_error << std::endl;
        if (rel_error < EPSILON) {
            break;
        }
    }

    double end_time = omp_get_wtime();
    double elapsed = end_time - start_time;
    
    // std::cout << "iterations = " << iterations_count << " ; ";
    // std::cout << "rel_error = " << rel_error << " ; ";
    // std::cout << "Time taken: " << elapsed << " sec.\n";

    // double diff = 0;
    // for (int i = 0; i < N; ++i) {
    //     for (int j = 0; j < N; ++j) {
    //         diff += (vector_x[i] - correct_x[i]) * (vector_x[i] - correct_x[i]);
    //     }
    // }
    // diff = std::sqrt(diff);
    // std::cout << "diff = " << diff << std::endl;

    std::cout << elapsed << std::endl;

    return 0;
}

// компилируй так: g++ -fopenmp -O1 -o parallel1 parallel1.cpp