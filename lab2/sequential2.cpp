#include <iostream>
#include <vector>
#include <cstring> // memcpy
#include <chrono>
#include <cmath>
#include <fstream> // std::ifstream

const int MAX_ITERATIONS = 10000;
const double EPSILON = 1e-8;
const double TAU = 0.004;
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

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<float> matrix_a(N * N);
    std::vector<float> vector_b(N);
    std::vector<float> vector_x(N);
    std::vector<float> next_x(N);

    // initialize(matrix_a, vector_b, next_x);

    if (!loadBinary("matA.bin", matrix_a, N * N) ||
        !loadBinary("vecB.bin", vector_b, N) ||
        !loadBinary("vecX.bin", vector_x, N)) {
        return 0;
    } // типа ошибка при загрузке

    std::memcpy(&next_x[0], &vector_x[0], N * sizeof(float));

    int iterations_count = 0;
    float current_norm = 0;
    
    float b_norm = 0;
    for (int i = 0; i < N; ++i) {
        b_norm += vector_b[i] * vector_b[i];
    }
    b_norm = std::sqrt(b_norm);
    
    for (; iterations_count < MAX_ITERATIONS; ++iterations_count) {
        current_norm = 0;
        for (int i = 0; i < N; ++i) {
            // sum = -b + A*x^n
            float sum = -vector_b[i];
            for (int j = 0; j < N; ++j) {
                sum += matrix_a[i * N + j] * vector_x[j];
            }
            // x^(n+1) = x^n - tau * sum
            next_x[i] = vector_x[i] - TAU * sum;
            // обновил нормы
            current_norm += sum * sum; 
        }
        // обновляем вектор приближения x
        std::memcpy(&vector_x[0], &next_x[0], N * sizeof(float));
        
        float rel_error = std::sqrt(current_norm) / b_norm;
        if (rel_error < EPSILON) {
            std::cout << "Relative error: " << rel_error << ", EPS: " << EPSILON << std::endl;
            break;
        } else {
            std::cout << "Iteration " << iterations_count << ": " << rel_error << std::endl;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // std::cout << next_x[0] << " ; ";
    std::cout << "iterations = " << iterations_count << " ; ";
    std::cout << "Time taken: " << elapsed.count() << " sec.\n";
    
    return 0;
}

// компилируй так: g++ -O1 -o sequential2 sequential2.cpp 