#include <iostream>
#include <vector>
#include <cstring>    // std::memcpy
#include <cmath>
#include <mpi.h>      // MPI для параллельных вычислений
#include <fstream>
// #include <chrono>

const int MAX_ITERATIONS = 100000;
const double EPSILON = 1e-2;
const double TAU = -0.01;
const int N = 2500;

bool loadBinary(const std::string& filename, std::vector<float>& data, size_t size) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cout << "ошибка при открытии файла: " << filename << std::endl;
        return false;
    }

    file.read(reinterpret_cast<char*>(data.data()), size * sizeof(float));
    if (!file) {
        std::cerr << "ошибка чтения файла: " << filename << std::endl;
        return false;
    }

    return true;
}

void iterate(std::vector<float>& matrix_a, std::vector<float>& vector_b, std::vector<float>& vector_x, int& iterations_count,
            int local_N, int local_offset, std::vector<int>& recvcounts, std::vector<int>& offsets) {
    float b_norm = 0;
    for (int i = 0; i < N; ++i) {
        b_norm += vector_b[i] * vector_b[i];
    }
    b_norm = std::sqrt(b_norm);

    std::vector<float> diffs(local_N, 0.0f); // остатки каждой строки матрицы сохраняем в этот вектор

    bool run = true;
    for (; iterations_count < MAX_ITERATIONS; ++iterations_count) {
        float local_norm = 0.0f;
        
        for (int i = 0; i < local_N; ++i) {
            int global_i = local_offset + i;
            float sum = -vector_b[global_i];
            for (int j = 0; j < N; ++j) {
                sum += matrix_a[global_i * N + j] * vector_x[j];
            }
            diffs[i] = sum;
            local_norm += sum * sum;
        }

        float global_norm;
        MPI_Allreduce(&local_norm, &global_norm, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        float rel_error = std::sqrt(global_norm) / b_norm;
        // std::cout << "iteration: " << iterations_count << "; rel_error: " << rel_error << std::endl;
        if (rel_error < EPSILON) {
            return;
        }

        for (int i = 0; i < local_N; ++i) {
            int global_i = local_offset + i;
            vector_x[global_i] -= TAU * diffs[i];
        }
        // сборка блоков x от всех процессов в один; получатели - все процессы вообще, но: 
        // recvcounts.data() - массив, где каждому процессу говорим сколько элементов он отправил
        // offsets.data() - массив, где каждому процессу указано его смещение
        MPI_Allgatherv(vector_x.data() + local_offset, local_N, MPI_FLOAT, 
                        vector_x.data(), recvcounts.data(), offsets.data(), MPI_FLOAT, MPI_COMM_WORLD);

    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // номер (ранг) процесса, вызвавшего функцию
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); // количчество процессов в области связи коммуникатора MPI_COMM_WORLD 

    //-----------------------------------------------------------------------------------------------------
    // rows_per_proc - столько строк получил бы каждый процесс если бы распределение было ровным
    int rows_per_proc = N / num_procs;
    // остаток от деления
    int remainder = N % num_procs;
    // число строк которое обрабатывает текущий процесс
    int local_N = rows_per_proc + (rank < remainder ? 1 : 0);
    // индекс первой строки, которую обрабатывает текущий процесс
    // либо rank*(rows_per_proc + 1), либо сдвиг с учетом доп строк
    int local_offset = (rank < remainder) ? rank * (rows_per_proc + 1)
                                          : remainder * (rows_per_proc + 1) + (rank - remainder) * rows_per_proc;
    // int local_offset = rank * rows_per_proc + (rank < remainder ? rank : remainder);

    // отдельно создание двух массивов для Allgatherv - модифицированный Allgather
    // recvcounts[p] – число строк вектора x, обновляемых процессом p
    // offsets[p] – смещение в итоговом массиве vector_x для данных от процесса p
    std::vector<int> recvcounts(num_procs), offsets(num_procs);
    int offset = 0;
    for (int p = 0; p < num_procs; p++) {
        int rows = rows_per_proc + (p < remainder ? 1 : 0);
        recvcounts[p] = rows;
        offsets[p] = offset;
        offset += rows;
    }
    //-----------------------------------------------------------------------------------------------------

    std::vector<float> matrix_a(N * N);
    std::vector<float> vector_b(N); // 
    std::vector<float> vector_x(N, 0.f); // глобальный вектор x

    if (!loadBinary("matA.bin", matrix_a, N * N) ||
        !loadBinary("vecB.bin", vector_b, N)) {
        MPI_Finalize();
        return 0;
    }

    int iterations_count = 0;

    auto start_time = MPI_Wtime();
    iterate(matrix_a, vector_b, vector_x, iterations_count, local_N, local_offset, recvcounts, offsets);
    auto end_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD); // здесь синхронизация чтобы вывод в консоль был ПОСЛЕДНЕЙ строкой
    if (rank == 0) {
        std::cout << "time: " << end_time - start_time << "; iterations: " << iterations_count  << std::endl;
    }

    MPI_Finalize(); // завершаем MPI (закрыли все MPI процесы, ликвидация всех областей связи)
    return 0;
}

/* 
компилируй так: mpic++ -O3 -o parallel1 parallel1.cpp
                mpirun ./parallel1

дебаг:          mpic++ -O3 -o parallel1 parallel1.cpp -g
                mpirun -np 4 xterm -fa 'Monospace' -fs 14 -e gdb -ex run parallel1
*/ 