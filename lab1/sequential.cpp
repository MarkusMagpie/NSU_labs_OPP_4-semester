#include <iostream>
#include <vector>
#include <cstring> // memcpy
#include <chrono>
#include <cmath>
#include <fstream> // std::ifstream
#include <mpi.h> // по приколу

const int MAX_ITERATIONS = 10000;
const double EPSILON = 1e-2;
const double TAU = -0.01;
const int N = 2500;

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

    return file.good();
}

int iterate(std::vector<float>& matrix_a, std::vector<float>& vector_b, std::vector<float>& vector_x, int iterations_count) {
    float b_norm = 0.0f;
    for (float val : vector_b) {
        b_norm += val * val;
    }
    b_norm = std::sqrt(b_norm);

    float current_norm;
    std::vector<float> diffs(N, 0.0f); // остатки каждой строки матрицы сохраняем в этот вектор

    for (; iterations_count < MAX_ITERATIONS; ++iterations_count) {
        current_norm = 0;
        

        for (int i = 0; i < N; ++i) {
            float sum = -vector_b[i];
            for (int j = 0; j < N; ++j) {
                sum += matrix_a[i * N + j] * vector_x[j];
            }
            diffs[i] = sum;
            current_norm += sum * sum;
        }

        float rel_error = std::sqrt(current_norm) / b_norm;
        // std::cout << "iteration: " << iterations_count << "; rel_error: " << rel_error << std::endl;
        if (rel_error < EPSILON) {
            break;
        }

        for (int i = 0; i < N; ++i) {
            vector_x[i] -= TAU * diffs[i];
        }
    }

    return iterations_count;
}

int main() {
    std::vector<float> matrix_a(N * N);
    std::vector<float> vector_b(N);
    std::vector<float> vector_x(N, 0.f); // инициализируем вектор_x нулями

    if (!loadBinary("matA.bin", matrix_a, N * N) ||
        !loadBinary("vecB.bin", vector_b, N)) {
        return 0;
    } // ошибка при загрузке

    auto start_time = std::chrono::high_resolution_clock::now();
    int iterations_count = iterate(matrix_a, vector_b, vector_x, 0);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    // double diff = 0;
    // for (int i = 0; i < N; ++i) {
    //     diff += (vector_x[i] - correct_x[i]) * (vector_x[i] - correct_x[i]);
    // }
    // diff = std::sqrt(diff);
    // std::cout << "diff = " << diff << std::endl;

    std::cout << elapsed.count() << "; " << iterations_count << std::endl;

    return 0;
}

// компилируй так: g++ -O3 -o sequential sequential.cpp 