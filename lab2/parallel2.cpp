#include <iostream>
#include <vector>
#include <cstring>    // std::memcpy
#include <cmath>
#include <omp.h>      // OpenMP для параллельных вычислений

const int MAX_ITERATIONS = 10000;
const double EPSILON = 0.00001;
const double TAU = 0.001;
const int N = 1990;

void initialize(std::vector<double>& matrix,
                std::vector<double>& vector_b,
                std::vector<double>& vector_x) {
    for (int i = 0; i < N; ++i) {
        vector_x[i] = 0;
        vector_b[i] = N + 1;
        for (int j = 0; j < N; ++j) {
            matrix[i * N + j] = (i == j) ? 2.0 : 1.0;
        }
    }
}

int main() {
    double start_time = omp_get_wtime();

    std::vector<double> matrix_a(N * N);
    std::vector<double> vector_b(N);
    std::vector<double> vector_x(N);

    initialize(matrix_a, vector_b, vector_x);
    std::fill(next_x.begin(), next_x.end(), 0.0);

    double b_norm = 0;
    for (int i = 0; i < N; ++i) {
        b_norm += vector_b[i] * vector_b[i];
    }
    b_norm = std::sqrt(b_norm);

    int iterations_count = 0;
    double current_norm = 0;

    bool stop = false;  // новый флаг для завершения итераций
    double local_norm;  // не приватна для каждого потока

    // одна параллельная секция, охватывающая весь итерационный цикл
    #pragma omp parallel
    {
        while (!stop) {
            local_norm = 0.0;
            #pragma omp for schedule(static) reduction(+ : local_norm)
            for (int i = 0; i < N; ++i) {
                double sum = -vector_b[i];
                for (int j = 0; j < N; ++j) {
                    sum += matrix_a[i * N + j] * vector_x[j];
                }
                next_x[i] = vector_x[i] - TAU * sum;
                local_norm += sum * sum;
            }

            #pragma omp master
            // только master тред обновляет переменные и проверяет условие выхода; если оно выполняется, то
            // все потоки выйдут из цикла в следующей итерации
            {
                current_norm = local_norm;
                std::memcpy(vector_x.data(), next_x.data(), N * sizeof(double));
                iterations_count++;

                double rel_error = std::sqrt(current_norm) / b_norm;
                // std::cout << "Iteration " << iterations_count << ": " << rel_error << std::endl;
                if (rel_error < EPSILON) {
                    stop = true;
                }
            }

            #pragma omp barrier // синхронизация потоков между итерациями - все ждут выполнения одной итерации и переходят к следующей
        }
    }

    double end_time = omp_get_wtime();
    double elapsed = end_time - start_time;

    std::cout << elapsed << std::endl;
    // std::cout << "iterations = " << iterations_count << std::endl;

    return 0;
}
