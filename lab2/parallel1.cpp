#include <iostream>
#include <vector>
#include <cstring>    // std::memcpy
#include <cmath>
#include <omp.h>      // OpenMP для параллельных вычислений

const int MAX_ITERATIONS = 10000;
const double EPSILON = 0.000001;
const double TAU = 0.001;
const int N = 1990;

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
    // функция из OpenMP для подсчета времени
    double start_time = omp_get_wtime();

    std::vector<double> matrix_a(N * N);
    std::vector<double> vector_b(N);
    std::vector<double> vector_x(N);
    std::vector<double> next_x(N);

    initialize(matrix_a, vector_b, vector_x);
    // std::memcpy(&next_x[0], &vector_x[0], N * sizeof(double));
    std::fill(next_x.begin(), next_x.end(), 0);

    int iterations_count = 0;
    double current_norm = 0;

    double b_norm = 0;
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
            long double sum = -vector_b[i];
            for (int j = 0; j < N; ++j) {
                sum += matrix_a[i * N + j] * vector_x[j];
            }
            next_x[i] = vector_x[i] - TAU * sum;
            current_norm += sum * sum;
        }

        std::memcpy(&vector_x[0], &next_x[0], N * sizeof(double));

        double rel_error = std::sqrt(current_norm) / b_norm;
        // std::cout << "Iteration " << iterations_count << ": " << rel_error << std::endl;
        if (rel_error < EPSILON) {
            break;
        }
    }

    double end_time = omp_get_wtime();
    double elapsed = end_time - start_time;

    // std::cout << next_x[0] << " ; ";
    
    // std::cout << "iterations = " << iterations_count << " ; ";
    // std::cout << "Time taken: " << elapsed << " sec.\n";

    std::cout << elapsed << std::endl;

    return 0;
}

// компилируй так: g++ -fopenmp -O1 -o parallel parallel1.cpp