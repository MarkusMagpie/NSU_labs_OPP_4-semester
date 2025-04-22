#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <cmath>

const double a = 1e2; // параметр уравнения
const double epsilon = 1e-16; // порог сходимости
const double D = 2.0; // область моделирования: [-1, 1] x [-1, 1] x [-1, 1]
const int Nx = 100, Ny = 100, Nz = 100; // параметры сетки
const double hx = D / (Nx - 1); // шаги сетки по различным осям (расстояния между узлами сетки)
const double hy = D / (Ny - 1);
const double hz = D / (Nz - 1);
const int max_iterarions = 10000;

double phi(double x, double y, double z) {
    return x*x + y*y + z*z;
}

double rho(double x, double y, double z) { 
    return 6.0 - a * phi(x, y, z); 
}

int idx(int i, int j, int k) {
    return i * Ny * Nz + j * Nz + k;
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

    // if (rank == 0) {
    //     for (int i = 0; i < Nx; i++) {
    //         for (int j = 0; j < Ny; j++) {
    //             for (int k = 0; k < Nz; k++) {
    //             std::cout << "(" << i << "," << j << "," << k << ") -> " << idx(i,j,k) << std::endl;
    //             }
    //         }
    //     }
    // }

    // 2 - инициализация данных
    // выделим память для локальных массивов phi_old и phi_new с учетом 2 теневых слоев 
    // теневые слои хранят данные соседних процессов (1 сдева и 1 справа)
    int shadow = 1;
    int local_nx = end_x - start_x + 1 + 2*shadow; // строки которые обрабатывает текущий процесс + 2 теневых слоя

    // 3d массивы 
    double* phi_old = new double[local_nx * Ny * Nz]();
    double* phi_new = new double[local_nx * Ny * Nz]();

    // инициализация краевых узлов в 3d массивах - ЕДИНОЖДЫ
    // УСЛОВИЕ: если узел на границе, то его значение будет равно phi(x,y,z)
    for (int i = 0; i < local_nx; ++i) {
        // глобальная координата X 
        int global_i = start_x + (i - shadow);
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                bool is_boundary = (global_i == 0 || global_i == Nx - 1 || j == 0 || j == Ny - 1 || k == 0 || k == Nz - 1);
                if (is_boundary) {
                    double x = -1.0 + global_i * hx;
                    double y = -1.0 + j * hy;
                    double z = -1.0 + k * hz;
                    phi_old[idx(i, j, k)] = phi(x, y, z);
                }
            }
        }
    }

    // 3 - итерационный процесс
    int iterations = 0;
    double max_diff, global_max_diff = epsilon + 1.0;
    MPI_Request requests[4];
    int request_count;

    double start_time = MPI_Wtime();

    for (; iterations < max_iterarions && global_max_diff > epsilon; ++iterations) {
        request_count = 0;
        // обмен граничными слоями - пары неблокирующих опрераций

        // (<-) я левому соседу отправляю мой первый реальный слой (shadow). левый сосед примет этот слой в свой правый ghost слой (local_nx-1)
        // (->) левый сосед отправляет мне свой последний реальный слой (local_nx-shadow-1) и я его пишу в свой левый ghost слой (0)
        if (rank > 0) {
            MPI_Isend(&phi_old[shadow * Ny * Nz], Ny * Nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Irecv(&phi_old[0], Ny * Nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
        }

        // (->) я правому соседу отправляю мой последний реальный слой (local_nx-shadow-1). правый сосед примет его в левый ghost слой (0)
        // (<-) правый сосед отправляет мне свой первый реальный слой (shadow) и я его пишу в свой правый ghost слой (local_nx-1)
        if (rank < size - 1) {
            MPI_Isend(&phi_old[(local_nx - shadow - 1) * Ny * Nz], Ny * Nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Irecv(&phi_old[(local_nx - 1) * Ny * Nz], Ny * Nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
        }

        // ждёт завершения всех запросов (отправок/приёмов) чтобы мог безопасно читать/писать в переданные буфера
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);

        // вычисление внутренних узлов (по оси X: с shadow+1 по local_nx-shadow-1; )
        max_diff = 0.0;
        for (int i = shadow + 1; i < local_nx - shadow; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int k = 1; k < Nz - 1; ++k) {
                    // глобальная координата X
                    int global_i = start_x + (i - shadow);

                    // формула для итерационного процесса метода Якоби
                    double term_x = (phi_old[idx(i+1,j,k)] + phi_old[idx(i-1,j,k)]) / (hx*hx);
                    double term_y = (phi_old[idx(i,j+1,k)] + phi_old[idx(i,j-1,k)]) / (hy*hy);
                    double term_z = (phi_old[idx(i,j,k+1)] + phi_old[idx(i,j,k-1)]) / (hz*hz);

                    // координаты узла с индексами i, j, k в физическом пространстве
                    double x = -1.0 + global_i * hx;
                    double y = -1.0 + j * hy;
                    double z = -1.0 + k * hz;

                    // double rho_new_x = (phi_old[(i+1)*Ny*Nz + j*Nz + k] - 2*phi_old[i*Ny*Nz + j*Nz + k] + phi_old[(i-1)*Ny*Nz + j*Nz + k]) / (hx*hx);
                    // double rho_new_y = (phi_old[i*Ny*Nz + (j+1)*Nz + k] - 2*phi_old[i*Ny*Nz + j*Nz + k] + phi_old[i*Ny*Nz + (j-1)*Nz + k]) / (hy*hy);
                    // double rho_new_z = (phi_old[i*Ny*Nz + j*Nz + (k+1)] - 2*phi_old[i*Ny*Nz + j*Nz + k] + phi_old[i*Ny*Nz + j*Nz + (k-1)]) / (hz*hz);
                    // double rho_new = (rho_new_x + rho_new_y + rho_new_z) - a * phi_old[i*Ny*Nz + j*Nz + k]; 

                    // формула якоби
                    double numerator = term_x + term_y + term_z - rho(x,y,z);
                    double denom = 2.0/(hx*hx) + 2.0/(hy*hy) + 2.0/(hz*hz) + a;
                    phi_new[idx(i,j,k)] = numerator / denom;

                    double diff = std::abs(phi_new[idx(i,j,k)] - phi_old[idx(i,j,k)]);
                    max_diff = std::max(max_diff, diff);
                }
            }
        }

        // по операции MPI_MAX собирает значения max_diff со всех процессов и выдаёт global_max_diff
        MPI_Allreduce(&max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        // if (rank == 0) std::cout << "Iteration: " << iterations << " Max diff: " << global_max_diff << std::endl;

        // std::swap(phi_old, phi_new); - ошибка
        // обновляем ТОЛЬКО внутренние узлы. границы гнать нахуй
        for (int i = 1; i < local_nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int k = 1; k < Nz - 1; ++k) {
                    phi_old[idx(i,j,k)] = phi_new[idx(i,j,k)];
                }
            }
        }
    }

    double end_time = MPI_Wtime();

    // 4 - вывод результатов и завершение MPI программы
    if (rank == 0) {
        std::cout << "итерации: " << iterations << std::endl;
        std::cout << "max diff: " << global_max_diff << "; epsilon: " << epsilon << std::endl;
        std::cout << "time: " << end_time - start_time << std::endl;
    }

    delete[] phi_old;
    delete[] phi_new;

    MPI_Finalize();
    return 0;
}

/* 
компилируй так: mpic++ -O3 -o parallel parallel.cpp
                mpirun  -np 4 ./parallel

дебаг:          mpic++ -O3 -o parallel parallel.cpp -g
                mpirun -np 4 xterm -fa 'Monospace' -fs 14 -e gdb -ex run parallel
*/ 
