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

int main() {
    // функция из OpenMP для подсчета времени
    double start_time = omp_get_wtime();

    std::vector<float> matrix_a(N * N);
    std::vector<float> vector_b(N);
    std::vector<float> vector_x(N);

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
                vector_x[i] -= TAU * sum;
                local_norm += sum * sum;
            }

            #pragma omp master
            // только master тред обновляет переменные и проверяет условие выхода; если оно выполняется, то
            // все потоки выйдут из цикла в следующей итерации
            {
                current_norm = local_norm;
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

    return 0;
}
