#include <iostream>
#include <vector>
#include <cstring> // memcpy
#include <chrono>
#include <cmath>

const int MAX_ITERATIONS = 10000;
const double EPSILON = 0.000001;
const double TAU = 0.001;
const int N = 1990; // размерность мтатрицы A

void initialize(std::vector<double>& matrix, std::vector<double>& vector_b, std::vector<double>& vector_x) {
    for (int i = 0; i < N; ++i) {
        vector_x[i] = 0;
        vector_b[i] = N + 1;

        for (int j = 0; j < N; ++j) {
            if (i == j) {
                matrix[i * N + j] = 2.0;
            } else {
                matrix[i * N + j] = 1.0;
            }
        }
    }
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> matrix_a(N * N);
    std::vector<double> vector_b(N);
    std::vector<double> vector_x(N);
    std::vector<double> next_x(N);

    initialize(matrix_a, vector_b, next_x);
    std::memcpy(&vector_x[0], &next_x[0], N * sizeof(double));

    int iterations_count = 0;
    double current_norm = 0;
    
    double b_norm = 0;
    for (int i = 0; i < N; ++i) {
        b_norm += vector_b[i] * vector_b[i];
    }
    b_norm = std::sqrt(b_norm);
    
    for (; iterations_count < MAX_ITERATIONS; ++iterations_count) {
        current_norm = 0;
        for (int i = 0; i < N; ++i) {
            // sum = -b + A*x^n
            long double sum = -vector_b[i];
            for (int j = 0; j < N; ++j) {
                sum += matrix_a[i * N + j] * vector_x[j];
            }
            // x^(n+1) = x^n - tau * sum
            next_x[i] = vector_x[i] - TAU * sum;
            // обновил нормы
            current_norm += sum * sum; 
        }
        // обновляем вектор приближения x
        std::memcpy(&vector_x[0], &next_x[0], N * sizeof(double));
        
        if ((std::sqrt(current_norm) / std::sqrt(b_norm)) < EPSILON) {
            break;
        } else {
            std::cout << "Iteration " << iterations_count << ": " << std::sqrt(current_norm) / std::sqrt(b_norm) << std::endl;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // std::cout << next_x[0] << " ; ";
    std::cout << "iterations = " << iterations_count << " ; ";
    std::cout << "Time taken: " << elapsed.count() << " sec.\n";
    
    return 0;
}

// компилируй так: g++ -O1 -o sequential sequential.cpp 