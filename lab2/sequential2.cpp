#include <iostream>
#include <vector>
#include <cstring> // memcpy
#include <chrono>
#include <cmath>
#include <fstream> // std::ifstream

const int MAX_ITERATIONS = 10000;
const double EPSILON = 1e-2;
const double TAU = -0.01;
const int N = 2500; // размерность мтатрицы A из фаайлика Матвееыва

// void initialize(std::vector<float>& matrix, std::vector<float>& vector_b, std::vector<float>& vector_x) {
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

void iterate(std::vector<float>& matrix_a, std::vector<float>& vector_b, std::vector<float>& vector_x, int iterations_count) {
    float b_norm = 0;
    for (int i = 0; i < N; ++i) {
        b_norm += vector_b[i] * vector_b[i];
    }
    b_norm = std::sqrt(b_norm);

    float current_norm;
    for (; iterations_count < MAX_ITERATIONS; ++iterations_count) {
        current_norm = 0;

        for (int i = 0; i < N; ++i) {
            float sum = -vector_b[i];
            for (int j = 0; j < N; ++j) {
                sum += matrix_a[i * N + j] * vector_x[j];
            }
            vector_x[i] -= TAU * sum;
            current_norm += sum * sum;
        }

        float rel_error = std::sqrt(current_norm) / b_norm;
        // std::cout << "iteration: " << iterations_count << "; rel_error: " << rel_error << std::endl;
        if (rel_error < EPSILON) {
            break;
        }
    }
}

int main() {
    std::vector<float> matrix_a(N * N);
    std::vector<float> vector_b(N);
    std::vector<float> vector_x(N, 0.f); // инициализируем вектор_x нулями

    // std::vector<float> correct_x(N);
    // loadBinary("vecX.bin", vector_x, N);

    if (!loadBinary("matA.bin", matrix_a, N * N) ||
        !loadBinary("vecB.bin", vector_b, N)) {
        return 0;
    } // ошибка при загрузке

    auto start_time = std::chrono::high_resolution_clock::now();
    iterate(matrix_a, vector_b, vector_x, 0);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    // double diff = 0;
    // for (int i = 0; i < N; ++i) {
    //     for (int j = 0; j < N; ++j) {
    //         diff += (vector_x[i] - correct_x[i]) * (vector_x[i] - correct_x[i]);
    //     }
    // }
    // diff = std::sqrt(diff);
    // std::cout << "diff = " << diff << std::endl;

    std::cout << elapsed.count() << std::endl;

    return 0;
}

// компилируй так: g++ -O1 -o sequential2 sequential2.cpp 