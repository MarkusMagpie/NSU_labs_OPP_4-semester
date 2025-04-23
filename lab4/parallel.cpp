#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <cmath>

const float a = 1e2; // параметр уравнения
const float epsilon = 1e-8; // порог сходимости
const float D_x = 2.0, D_y = 2.0, D_z = 2.0; // область моделирования: [-1, 1] x [-1, 1] x [-1, 1]
const float x_0 = -1.0, y_0 = -1.0, z_0 = -1.0;
const int Nx = 100, Ny = 100, Nz = 100; // параметры сетки
const float hx = D_x / (Nx - 1); // шаги сетки по различным осям (расстояния между узлами сетки)
const float hy = D_y / (Ny - 1);
const float hz = D_z / (Nz - 1);
const int max_iterarions = 10000;

float phi(float x, float y, float z) {
    return x*x + y*y + z*z;
}

float rho(float x, float y, float z) { 
    return 6.0 - a * phi(x, y, z); 
}

int idx(int i, int j, int k) {
    return i * Ny * Nz + j * Nz + k;
}

void for_slice(int i, 
               float* phi_old, float* phi_new, float* rho_vals,
               float& max_diff, int start_x, int shadow) {
    for (int j = 1; j < Ny - 1; ++j) {
        for (int k = 1; k < Nz - 1; ++k) {
            float term_x = (phi_old[idx(i+1,j,k)] + phi_old[idx(i-1,j,k)]) / (hx*hx);
            float term_y = (phi_old[idx(i,j+1,k)] + phi_old[idx(i,j-1,k)]) / (hy*hy);
            float term_z = (phi_old[idx(i,j,k+1)] + phi_old[idx(i,j,k-1)]) / (hz*hz);

            float numerator = term_x + term_y + term_z - rho_vals[idx(i,j,k)];
            float denom = 2.0/(hx*hx) + 2.0/(hy*hy) + 2.0/(hz*hz) + a;
            phi_new[idx(i,j,k)] = numerator / denom;

            float diff = std::abs(phi_new[idx(i,j,k)] - phi_old[idx(i,j,k)]);
            max_diff = std::max(max_diff, diff);
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv); // создаем предопределеннную область связи, содержащую все процессы MPI программы, с ней связывается коммуникатор MPI_COMM_WORLD
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // номер (ранг) процесса, вызвавшего функцию
    MPI_Comm_size(MPI_COMM_WORLD, &size); // количество процессов в области связи коммуникатора MPI_COMM_WORLD 

    // 1 - декомпозиция пространства - метод распараллеливания 
    int chunk_size = Nx / size;
    int rem = Nx % size; // если Nx % num_procs != 0, первые remainder процессов получают +1 узел
    int start_x = rank * chunk_size + std::min(rank, rem); // начальный индекс узлов процесса в глобальной сетке
    int end_x = start_x + chunk_size - 1;
    if (rank < rem) end_x += 1;
    if (rank == size - 1) end_x = Nx - 1; // коректирую последний процесс

    std::cout << "start_x = " << start_x << ", end_x = " << end_x << std::endl;

    // 2 - инициализация данных
    // выделим память для локальных массивов phi_old и phi_new с учетом 2 теневых слоев 
    // теневые слои хранят данные соседних процессов (1 сдева и 1 справа)
    int shadow = 1;
    int local_nx = end_x - start_x + 1 + 2*shadow; // строки которые обрабатывает текущий процесс + 2 теневых слоя

    // 3d массивы 
    float* phi_old = new float[local_nx * Ny * Nz]();
    float* phi_new = new float[local_nx * Ny * Nz]();
    float* rho_vals = new float[local_nx * Ny * Nz]();

    // инициализация краевых узлов в 3d массивах - ЕДИНОЖДЫ
    // УСЛОВИЕ: если узел на границе, то его значение будет равно phi(x,y,z)
    for (int i = 0; i < local_nx; ++i) {
        // глобальная координата X 
        int global_i = start_x + (i - shadow);
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                float x = x_0 + global_i * hx;
                float y = y_0 + j * hy;
                float z = z_0 + k * hz;
                rho_vals[idx(i,j,k)] = rho(x, y, z);
                bool boundary = (global_i == 0 || global_i == Nx - 1 || j == 0 || j == Ny - 1 || k == 0 || k == Nz - 1);
                if (boundary) {
                    phi_old[idx(i, j, k)] = phi(x, y, z);
                }
            }
        }
    }

    // 3 - итерационный процесс
    int iterations = 0;
    float max_diff, global_max_diff = epsilon + 1.0;
    MPI_Request requests[4];
    int request_count;

    float start_time = MPI_Wtime();

    do {
        request_count = 0;
        // обмен граничными слоями - пары неблокирующих опрераций

        // (<-) я левому соседу отправляю мой первый реальный слой (shadow). левый сосед примет этот слой в свой правый ghost слой (local_nx-1)
        // (->) левый сосед отправляет мне свой последний реальный слой (local_nx-shadow-1) и я его пишу в свой левый ghost слой (0)
        if (rank > 0) {
            MPI_Isend(&phi_old[shadow * Ny * Nz], Ny * Nz, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Irecv(&phi_old[0], Ny * Nz, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
        }

        // (->) я правому соседу отправляю мой последний реальный слой (local_nx-shadow-1). правый сосед примет его в левый ghost слой (0)
        // (<-) правый сосед отправляет мне свой первый реальный слой (shadow) и я его пишу в свой правый ghost слой (local_nx-1)
        if (rank < size - 1) {
            MPI_Isend(&phi_old[(local_nx - shadow - 1) * Ny * Nz], Ny * Nz, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Irecv(&phi_old[(local_nx - 1) * Ny * Nz], Ny * Nz, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
        }

        // вычисление внутренних узлов (по оси X: с shadow+1 по local_nx-shadow-1; )
        max_diff = 0.0;
        int i;
        for (i = shadow + 2; i < local_nx - shadow - 2; ++i) {
            for_slice(i, phi_old, phi_new, rho_vals, max_diff, start_x, shadow);
        }

        // ждёт завершения всех запросов (отправок/приёмов) чтобы мог безопасно читать/писать в переданные буфера
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);

        i = shadow; // i=1
        for_slice(i, phi_old, phi_new, rho_vals, max_diff, start_x, shadow);

        i = local_nx - shadow - 1; // i=local_nx-2
        for_slice(i, phi_old, phi_new, rho_vals, max_diff, start_x, shadow);

        // по операции MPI_MAX собирает значения max_diff со всех процессов и выдаёт global_max_diff
        MPI_Allreduce(&max_diff, &global_max_diff, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
        // if (rank == 0) std::cout << "Iteration: " << iterations << " Max diff: " << global_max_diff << std::endl;

        // std::swap(phi_old, phi_new);
        for (int i = shadow; i < local_nx - shadow; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int k = 1; k < Nz - 1; ++k) {
                    phi_old[idx(i,j,k)] = phi_new[idx(i,j,k)];
                }
            }
        }
        ++iterations;
    } while (global_max_diff > epsilon && iterations < max_iterarions);

    float end_time = MPI_Wtime();

    // 4 - вывод результатов и завершение MPI программы
    if (rank == 0) {
        std::cout << "итерации: " << iterations << std::endl;
        std::cout << "max diff: " << global_max_diff << "; epsilon: " << epsilon << std::endl;
        std::cout << "time: " << end_time - start_time << std::endl;
    }

    delete[] phi_old;
    delete[] phi_new;
    delete[] rho_vals;

    MPI_Finalize();
    return 0;
}

/* 
компилируй так: mpic++ -O3 -o parallel parallel.cpp
                mpirun  -np 4 ./parallel

дебаг:          mpic++ -O3 -o parallel parallel.cpp -g
                mpirun -np 4 xterm -fa 'Monospace' -fs 14 -e gdb -ex run parallel
*/ 
