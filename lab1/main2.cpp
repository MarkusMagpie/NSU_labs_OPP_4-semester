#include <mpi.h>
#include <mpe.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib> // atoi

const int DEFAULT_N = 100;
const double DEFAULT_TAU = 0.01;
const double EPS = 1e-5; // критерий завершения
const int DEFAULT_MAX_ITER = 1e+5; // максимальное число итерации

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

// функция для вычисления начального распределения строк (с учетом остатка)
// Возвращает пару: (local_N, global_start)
// local_N = число строк, принадлежащих текущему процессу,
// global_start = глобальный индекс первой строки, принадлежащей процессу.
std::pair<int, int> computeLocalRows(int N, int num_procs, int rank) {
    int base = N / num_procs;
    int rem = N % num_procs;
    int local_N = (rank < rem) ? base + 1 : base;
    int global_start = (rank < rem) ? rank * (base + 1)
                                    : rem * (base + 1) + (rank - rem) * base;
    return {local_N, global_start};
}

// 1 вариант - векторы x и b дублируются (полны на всех процессах)
int solveIterationVariant1(int N, double tau, int max_iter,
                           const std::vector<std::vector<double>>& A,
                           std::vector<double>& x,        // полный вектор x (размер N)
                           const std::vector<double>& b,  // полный вектор b (размер N)
                           int local_N) {
    int iter = 0;
    bool converged = false;
    while (!converged && iter < max_iter) {
        // каждый процесс вычисляет свою часть произведения: Ax_local = A * x
        std::vector<double> Ax_local = MatVectMult(A, x, local_N);
        // собираем полный вектор Ax, его получат все процессы
        std::vector<double> Ax(N, 0.0);
        MPI_Allgather(Ax_local.data(), local_N, MPI_DOUBLE,
                      Ax.data(), local_N, MPI_DOUBLE, MPI_COMM_WORLD);
        // обновляем полный вектор x: x = x - tau*(Ax - b)
        for (int i = 0; i < N; ++i) {
            x[i] = x[i] - tau * (Ax[i] - b[i]);
        }
        // ||Ax - b||
        double local_norm = 0.0;
        for (int i = 0; i < N; ++i) {
            double diff = Ax[i] - b[i];
            local_norm += diff * diff;
        }
        double global_norm = 0.0;
        MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        double b_norm = 0.0;
        for (double val : b) {
            b_norm += val * val;
        }

        double rel_err = std::sqrt(global_norm) / std::sqrt(b_norm);
        if (rel_err < EPS) {
            converged = true;
        }

        iter++;
    }

    return iter;
}

int solveIterationVariant2(int N, double tau, int max_iter,
                           const std::vector<std::vector<double>>& A,
                           std::vector<double>& x_local,   // локальный вектор x (размер local_N)
                           const std::vector<double>& b_local, // -//-
                           int local_N, int num_procs) {
    int iter = 0;
    bool converged = false;
    while (!converged && iter < max_iter) {
        // СНАЧАЛА собрать полный вектор x, уже затем можем A*full_x
        std::vector<double> full_x(N, 0.0);
        MPI_Allgather(x_local.data(), local_N, MPI_DOUBLE,
                      full_x.data(), local_N, MPI_DOUBLE, MPI_COMM_WORLD);


        std::vector<double> Ax_local = MatVectMult(A, full_x, local_N);

        // x_local = x_local - tau*(Ax_local - b_local)
        for (int i = 0; i < local_N; ++i) {
            x_local[i] = x_local[i] - tau * (Ax_local[i] - b_local[i]);
        }

        // вычисляем локальную норму ошибки
        double local_norm = 0.0;
        for (int i = 0; i < local_N; ++i) {
            double diff = Ax_local[i] - b_local[i];
            local_norm += diff * diff;
        }

        double global_norm = 0.0;
        MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        if (std::sqrt(global_norm) < EPS) {
            converged = true;
        }

        iter++;
    }

    return iter;
}

int main(int argc, char** argv) {
    // инициализация MPI и профайлера MPE
    MPI_Init(&argc, &argv);
    MPE_Init_log();

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // ранг процесса
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); // общее число процессов

    // считываем режим работы и параметры из командной строки
    int mode = (argc > 1) ? std::atoi(argv[1]) : 1;
    int N = (argc > 2) ? std::atoi(argv[2]) : DEFAULT_N;
    double tau = (argc > 3) ? std::atof(argv[3]) : DEFAULT_TAU;
    int max_iter = (argc > 4) ? std::atoi(argv[4]) : DEFAULT_MAX_ITER;

    if (mode != 1 && mode != 2) {
        if (rank == 0) {
            std::cerr << "Invalid mode. Use 1 or 2!" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        std::cout << "Mode: " << mode << std::endl;
        std::cout << "Matrix size N = " << N << ", tau = " << tau << ", max_iter = " << max_iter << std::endl;
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

    // Вычисляем число строк для текущего процесса с учётом остатка
    auto [local_N, global_start] = computeLocalRows(N, num_procs, rank);

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

    // вектора x и вектора b зависят от режима:
    // 1: векторы дублируются – собирается полный вектор b и x (размер N)
    // 2: векторы разрезаются – каждый процесс хранит свою часть
    std::vector<double> full_b; // variant1
    std::vector<double> full_x; // variant1
    std::vector<double> x_local; // variant2

    if (mode == 1) {
        full_b.resize(N, 0.0);
        // собираем полный вектор b из b_local
        MPI_Allgather(b_local.data(), local_N, MPI_DOUBLE,
                      full_b.data(), local_N, MPI_DOUBLE, MPI_COMM_WORLD);
        full_x.resize(N, 0.0);
        // рассылаем начальное приближение full_x (все нули)
        MPI_Bcast(full_x.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        x_local.resize(local_N, 0.0);
        // b_local и x_local остаются локальными
    }

    MPE_Log_event(evtid_end_init, rank, (char*)"" );

    // ФАЗА ИТЕРАЦИОННОГО ПРОЦЕССА
    MPE_Log_event(evtid_begin_iter, rank, (char*)"" );
    double start_time = MPI_Wtime();
    int iter = 0;

    if (mode == 1) {
        iter = solveIterationVariant1(N, tau, max_iter, A, full_x, full_b, local_N);
    } else {
        iter = solveIterationVariant2(N, tau, max_iter, A, x_local, b_local, local_N, num_procs);
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
