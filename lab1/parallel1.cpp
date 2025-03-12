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

bool loadBinaryPart(const std::string& filename, std::vector<float>& data, size_t offset, size_t count) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cout << "ошибка при открытии файла: " << filename << std::endl;
        return false;
    }
    
    // НОВОЕ - перемещаю указатель на offset * sizeof(float)
    file.seekg(offset * sizeof(float), std::ios::beg);
    data.resize(count);

    file.read(reinterpret_cast<char*>(data.data()), count * sizeof(float));
    if (!file) {
        std::cerr << "ошибка чтения файла: " << filename << std::endl;
        return false;
    }

    return file.good(); // https://cplusplus.com/reference/ios/ios/good/
}

int iterate(const std::vector<float>& local_a, std::vector<float>& b, std::vector<float>& x, 
             int& iterations_count, int local_n, int offset, const std::vector<int>& recvcounts, 
             const std::vector<int>& displs) {
    float b_norm = 0.0f;
    for (float val : b) {
        b_norm += val * val;
    }
    b_norm = std::sqrt(b_norm);

    std::vector<float> local_diffs(local_n); // остатки каждой строки матрицы сохраняем в этот вектор

    for (; iterations_count < MAX_ITERATIONS; ++iterations_count) {
        float local_norm = 0.0f;
        for (int i = 0; i < local_n; ++i) {
            float sum = -b[offset + i];
            for (int j = 0; j < N; ++j) {
                sum += local_a[i * N + j] * x[j];
            }
            local_diffs[i] = sum;
            local_norm += sum * sum;
        }

        // собираю суммой нормы на всех процессах и проверяю условие выхода из итерационного цикла
        float total_norm;
        MPI_Allreduce(&local_norm, &total_norm, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        float rel_error = std::sqrt(total_norm) / b_norm;
        // std::cout << "iteration: " << iterations_count << "; rel_error: " << rel_error << std::endl;
        if (rel_error < EPSILON) { 
            break;
        }

        for (int i = 0; i < local_n; ++i) {
            x[offset + i] -= TAU * local_diffs[i];
        }

        // сбор x на всех процессах
        MPI_Allgatherv(&x[offset], local_n, MPI_FLOAT, 
            x.data(), recvcounts.data(), displs.data(), MPI_FLOAT, MPI_COMM_WORLD);
    }

    return iterations_count;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv); // создаем предопределеннную область связи, содержащую все процессы MPI программы, с ней связывается коммуникатор MPI_COMM_WORLD

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // номер (ранг) процесса, вызвавшего функцию
    MPI_Comm_size(MPI_COMM_WORLD, &size); // количчество процессов в области связи коммуникатора MPI_COMM_WORLD 

    // распределение строк матрицы A между процессами
    int base = N / size;
    int rem = N % size;
    std::vector<int> recvcounts(size), displs(size);
    int offset = 0;
    for (int p = 0; p < size; ++p) {
        recvcounts[p] = base + (p < rem);
        displs[p] = offset;
        offset += recvcounts[p];
    }
    int local_n = recvcounts[rank]; // число строк которое обрабатывает текущий процесс
    int local_offset = displs[rank]; // индекс первой строки, которую обрабатывает текущий процесс

    // загрузка части матрицы A, local_offset и local_n зависят от процечча
    std::vector<float> local_a(local_n * N);
    if (!loadBinaryPart("matA.bin", local_a, local_offset * N, local_n * N)) {
        MPI_Finalize();
        return 0;
    }

    // загрузка вектора b и его рассылка 
    std::vector<float> b(N), x(N, 0.0f);
    if (rank == 0 && !loadBinaryPart("vecB.bin", b, 0, N)) {
        MPI_Finalize();
        return 0;
    }
    MPI_Bcast(b.data(), N, MPI_FLOAT, 0, MPI_COMM_WORLD);

    int iterations_count = 0;
    double start = MPI_Wtime();
    iterations_count = iterate(local_a, b, x, iterations_count, local_n, local_offset, recvcounts, displs);
    double elapsed = MPI_Wtime() - start;

    // MPI_Barrier(MPI_COMM_WORLD); // здесь синхронизация чтобы вывод в консоль был ПОСЛЕДНЕЙ строкой
    
    if (rank == 0) {
        std::cout << elapsed << "; " << iterations_count << std::endl;
    }

    MPI_Finalize();
    return 0;
}

/* 
компилируй так: mpic++ -O3 -o parallel1 parallel1.cpp
                mpirun ./parallel1

дебаг:          mpic++ -O3 -o parallel1 parallel1.cpp -g
                mpirun -np 4 xterm -fa 'Monospace' -fs 14 -e gdb -ex run parallel1
*/ 
