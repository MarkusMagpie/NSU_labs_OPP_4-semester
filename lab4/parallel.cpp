#include <mpi.h>
#include <algorithm>
#include <iostream>

const double a = 1e5; // параметр уравнения
const double epsilon = 1e-8; // порог сходимости
const double D = 2.0; // область моделирования: [-1, 1] x [-1, 1] x [-1, 1]
const int Nx = 100, Ny = 100, Nz = 100; // параметры сетки
const double hx = D / (Nx - 1); // шаги сетки по различным осям
const double hy = D / (Ny - 1);
const double hz = D / (Nz - 1);
const int max_iterarions = 200;

double phi(double x, double y, double z) {
    return x*x + y*y + z*z;
}

double rho(double x, double y, double z) { 
    return 6.0 - a * phi(x, y, z); 
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

    double* phi_old = new double[local_nx * Ny * Nz]();
    double* phi_new = new double[local_nx * Ny * Nz]();

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            // если узел на границе глобальной области - 2 варианта - либо он слева, либо справа
            // левый крайний ли этот процесс?
            if (start_x == 0) {
                double x = -1.0 + start_x * hx;
                double y = -1.0 + j * hy;
                double z = -1.0 + k * hz;
                // 1 - первый реальный слой, а 0 - теневой. поэтому и (shadow * Ny * Nz)
                phi_old[shadow * Ny * Nz + j * Nz + k] = phi(x, y, z);
            } 

            // это последний/правый крайний процесс? да - по индексу последнего реального слоя пишу фи
            if (end_x == Nx - 1) {
                double x = -1.0 + end_x * hx;
                double y = -1.0 + j * hy;
                double z = -1.0 + k * hz;
                // индекс: (local_nx - shadow - 1) - индекс последнего реального слоя (после него - теневой)
                phi_old[(local_nx - shadow - 1) * Ny * Nz + j * Nz + k] = phi(x, y, z);
            }   
        }
    }

    // 3 - итерационный процесс
    int iterations = 0;
    double max_diff, global_max_diff;
    MPI_Request requests[4];
    int request_count;

    double start_time = MPI_Wtime();

    do {
        request_count = 0;
        // обмен граничными слоями
        // отсылка левого теневого слоя соседу слева (shadow-1) и получение от него
        if (rank > 0) {
            // отправляем левую границу (индекс shadow) соседу слева
            MPI_Isend(&phi_old[shadow * Ny * Nz], Ny * Nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
            // получаем левый теневой (индекс 0) от соседа слева
            MPI_Irecv(&phi_old[0], Ny * Nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
        }

        // отправка правой границы соседу справа (rank+1) и получение от него
        if (rank < size - 1) {
            // отправляем правую границу (индекс local_nx - shadow - 1) соседу справа
            MPI_Isend(&phi_old[(local_nx - shadow - 1) * Ny * Nz], Ny * Nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
            // получаем правый теневой (индекс local_nx - 1)
            MPI_Irecv(&phi_old[(local_nx - 1) * Ny * Nz], Ny * Nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
        }

        // вычисление внутренних узлов (по оси X: с shadow+1 по local_nx-shadow-1; )
        max_diff = 0.0;
        for (int i = shadow + 1; i < local_nx - shadow; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int k = 1; k < Nz - 1; ++k) {
                    // глобальная координата X
                    int global_i = start_x + (i - shadow);

                    // формула для итерационного процесса метода Якоби
                    double term_x = (phi_old[(i+1)*Ny*Nz + j*Nz + k] + phi_old[(i-1)*Ny*Nz + j*Nz + k]) / (hx*hx);
                    double term_y = (phi_old[i*Ny*Nz + (j+1)*Nz + k] + phi_old[i*Ny*Nz + (j-1)*Nz + k]) / (hy*hy);
                    double term_z = (phi_old[i*Ny*Nz + j*Nz + (k+1)] + phi_old[i*Ny*Nz + j*Nz + (k-1)]) / (hz*hz);

                    // вычисление координат точки в физическом пространстве
                    double x = -1.0 + global_i * hx;
                    double y = -1.0 + j * hy;
                    double z = -1.0 + k * hz;

                    // double rho_new_x = (phi_old[(i+1)*Ny*Nz + j*Nz + k] - 2*phi_old[i*Ny*Nz + j*Nz + k] + phi_old[(i-1)*Ny*Nz + j*Nz + k]) / (hx*hx);
                    // double rho_new_y = (phi_old[i*Ny*Nz + (j+1)*Nz + k] - 2*phi_old[i*Ny*Nz + j*Nz + k] + phi_old[i*Ny*Nz + (j-1)*Nz + k]) / (hy*hy);
                    // double rho_new_z = (phi_old[i*Ny*Nz + j*Nz + (k+1)] - 2*phi_old[i*Ny*Nz + j*Nz + k] + phi_old[i*Ny*Nz + j*Nz + (k-1)]) / (hz*hz);
                    // double rho_new = (rho_new_x + rho_new_y + rho_new_z) - a * phi_old[i*Ny*Nz + j*Nz + k]; 
                    
                    phi_new[i*Ny*Nz + j*Nz + k] = (term_x + term_y + term_z - rho(x, y, z)) / (2/(hx*hx) + 2/(hy*hy) + 2/(hz*hz) + a);

                    double diff = std::abs(phi_new[i*Ny*Nz + j*Nz + k] - phi_old[i*Ny*Nz + j*Nz + k]);
                    max_diff = std::max(max_diff, diff);
                }
            }
        }

        // ждёт завершения отправок/приёмов теневых слоёв прежде чем перейти к следующей итерации
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
        // по операции MPI_MAX собирает значения max_diff со всех процессов и выдаёт global_max_diff
        MPI_Allreduce(&max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) std::cout << "Iteration: " << iterations << " Max diff: " << global_max_diff << std::endl;

        std::swap(phi_old, phi_new);
        ++iterations;
    } while (global_max_diff > epsilon && iterations < max_iterarions);

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
