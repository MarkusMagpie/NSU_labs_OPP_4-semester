#include <mpi.h>
#include <mpe.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib> // atoi

const double EPS = 1e-5; // критерий завершения
const double TAU = 0.01; // коэффициент

// умножение подматрицы A на вектор x
std::vector<double> MatVectMult(const std::vector<std::vector<double>>& A, 
                                const std::vector<double>& x, int rows) {
    int cols = x.size();
    std::vector<double> result(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

int main(int argc, char** argv) {
    // инициализация MPI и профайлера MPE
    MPI_Init(&argc, &argv);
    MPE_Init_log();

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // ранг процесса
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); // общее число процессов

    int mode = 1;
    if (argc > 1) {
        mode = std::atoi(argv[1]);
        if (mode != 1 && mode != 2) {
            if (rank == 0) {
                std::cerr << "Invalid mode. Use 1 or 2." << std::endl;
            }
            MPI_Finalize();
            return 0;
        }
    }
    
    if (rank == 0) {
        std::cout << "Mode: " << mode << std::endl;
        if (mode == 1) {
            std::cout << "Vectors x and b are duplicated on all MPI processes." << std::endl;
        } else {
            std::cout << "Vectors x and b are partitioned among MPI processes." << std::endl;
        }
    }

    // идентификаторы событий для профилирования
    int evtid_begin_init = MPE_Log_get_event_number();
    int evtid_end_init   = MPE_Log_get_event_number();
    int evtid_begin_iter = MPE_Log_get_event_number();
    int evtid_end_iter   = MPE_Log_get_event_number();

    // описание фаз профилирования
    MPE_Describe_state(evtid_begin_init, evtid_end_init, "Initialization", "green");
    MPE_Describe_state(evtid_begin_iter, evtid_end_iter, "Iteration", "blue");

    // ФАЗА ИНИЦИАЛИЗАЦИИ
    MPE_Log_event(evtid_begin_init, rank, (char*)"" );

    int N = 100;              // размерность матрицы A (N x N)
    int local_N = N / num_procs; // количество строк на каждый процесс

    // Каждый процесс выделяет память только для своей части матрицы A
    std::vector<std::vector<double>> A(local_N, std::vector<double>(N, 1.0));
    int global_row_start = rank * local_N;
    for (int i = 0; i < local_N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (global_row_start + i == j) {
                A[i][j] = 2.0;
            } else {
                A[i][j] = 1.0;
            }
        }
    }

    // u-вектор, элементы которого вычисляются по формуле: u_i = sin(2 * Pi * i / N)
    std::vector<double> u(N);
    for (int i = 0; i < N; ++i) {
        u[i] = std::sin(2 * M_PI * i / N);
    }

    // ыычисляем локальный вектор b = A * u
    std::vector<double> b_local = MatVectMult(A, u, local_N);
    // собираем глобальный вектор b на всех процессах
    std::vector<double> global_b(N, 0.0);
    MPI_Allgather(b_local.data(), local_N, MPI_DOUBLE, global_b.data(), local_N, MPI_DOUBLE, MPI_COMM_WORLD);

    // Инициализируем вектор x (начальное приближение x = 0)
    std::vector<double> global_x(N, 0.0);
    MPI_Bcast(global_x.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // конец фазы инициализации
    MPE_Log_event(evtid_end_init, rank, (char*)"" );

    // ФАЗА ИТЕРАЦИОННОГО ПРОЦЕССА
    // начало фазы фикируем
    MPE_Log_event(evtid_begin_iter, rank, (char*)"" );

    // начало замера времени итерационного процесса
    double start_time = MPI_Wtime();

    bool converged = false;
    int iter = 0;
    while (!converged) {
        // локальное произведение матрицы и вектора: Ax_local = A * global_x
        std::vector<double> Ax_local = MatVectMult(A, global_x, local_N);

        // все процесс получают полный вектор Ax
        std::vector<double> Ax(N, 0.0);
        MPI_Allgather(Ax_local.data(), local_N, MPI_DOUBLE, Ax.data(), local_N, MPI_DOUBLE, MPI_COMM_WORLD);

        // обновляем x по формуле: x^(n+1) = x^n - TAU * (Ax^n - b)
        for (int i = 0; i < N; ++i) {
            global_x[i] = global_x[i] - TAU * (Ax[i] - global_b[i]);
        }

        // НИЖЕ - УСЛОВИЯ ЗАВЕРШЕНИЯ ИТЕРАЦИИ
        // Вычисляем норму ||Ax - b||
        double norm_Ax_minus_b = 0.0;
        for (int i = 0; i < N; ++i) {
            norm_Ax_minus_b += (Ax[i] - global_b[i]) * (Ax[i] - global_b[i]);
        }
        // Вычисляем норму ||b||
        double norm_b = 0.0;
        for (int i = 0; i < N; ++i) {
            norm_b += global_b[i] * global_b[i];
        }

        double global_norm_Ax_minus_b, global_norm_b;
        MPI_Allreduce(&norm_Ax_minus_b, &global_norm_Ax_minus_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&norm_b, &global_norm_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double relative_error = std::sqrt(global_norm_Ax_minus_b) / std::sqrt(global_norm_b);
        if (relative_error < EPS) {
            converged = true;
        }

        // считать итерации на одном процессе
        if (rank == 0) {
            iter++;
        }
    }

    double end_time = MPI_Wtime();
    double total_time = end_time - start_time;

    // Фиксируем конец итерационного этапа
    MPE_Log_event(evtid_end_iter, rank, (char*)"" );

    // собираем время работы (минимальное среди процессов)
    double total_time_all;
    MPI_Reduce(&total_time, &total_time_all, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Converged in " << iter << " iterations." << std::endl;
        std::cout << "Total time: " << total_time_all << " seconds." << std::endl;
    }

    // завершаем профилирование и записываем лог в файл "profile.clog2"
    MPE_Finish_log("profile");

    MPI_Finalize();
    return 0;
}
