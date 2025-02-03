#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

const double EPS = 1e-5; // критерий завершения
const double TAU = 0.01; // коэффициент

// умножение подматрицы A на вектор x
std::vector<double> MatVectMult(const std::vector<std::vector<double>>& A, const std::vector<double>& x, int start_row, int rows) {
    int cols = x.size();
    std::vector<double> result(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i] += A[start_row + i][j] * x[j];
        }
    }
    return result;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ранг текущего процеесса в коммуникаторе
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); // общее количество процессов в коммуникаторе

    int N = 100; // размерность матрицы A (N на N)
    int local_N = N / num_procs; // количество строк на каждый процесс

    // инициализация данных в матрице А (элементы главной диагонали = 2.0, остальные = 1.0)
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

    // // рассылаем матрицу A на всех процессах помимо главного процесса 0
    // MPI_Bcast(A.data(), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // вектор u, элементы которого заполняются по формуле u_i = sin(2 * Pi * i) / N
    std::vector<double> u(N);
    for (int i = 0; i < N; ++i) {
        u[i] = std::sin(2 * M_PI * i / N);
    }

    // вектор b = A * u (матрица A умножается на вектор u)
    // std::vector<double> global_b(N, 0.0);
    // if (rank == 0) {
    //     global_b = MatVectMult(A, u, 0, N);
    // }

    // // разделяем вектор b между процессами
    // MPI_Bcast(global_b.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Вычисляем локальный b = A * u
    std::vector<double> b_local = MatVectMult(A, u, 0, local_N);
    // Собираем глобальный вектор b на всех процессах:
    std::vector<double> global_b(N, 0.0);
    MPI_Allgather(b_local.data(), local_N, MPI_DOUBLE, global_b.data(), local_N, MPI_DOUBLE, MPI_COMM_WORLD);
    
    std::vector<double> global_x(N, 0.0); // начальные значения элементов x = 0

    // разделяем вектор x на всех процессах (для варианта 1 или 2???)
    MPI_Bcast(global_x.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // начало замера времени работы итерационной части программы
    double start_time = MPI_Wtime();

    // итерации метода
    bool converged = false;
    int iter = 0;
    while (!converged) {
        // локальное произведение матрицы и вектора
        std::vector<double> Ax_local = MatVectMult(A, global_x, 0, local_N); // A содержит только local_N строк

        // собираем Ax на всех процессах
        std::vector<double> Ax(N, 0.0);
        MPI_Allgather(Ax_local.data(), local_N, MPI_DOUBLE, Ax.data(), local_N, MPI_DOUBLE, MPI_COMM_WORLD);

        // обновляем x (формула: x^(n+1) = x^n - TAU * (Ax^n - b))
        for (int i = 0; i < N; ++i) {
            double new_x = global_x[i] - TAU * (Ax[i] - global_b[i]);
            global_x[i] = new_x;
        }

        // вычисление нормы ||Ax^n - b||
        double norm_Ax_minus_b = 0.0;
        for (int i = 0; i < N; ++i) {
            norm_Ax_minus_b += (Ax[i] - global_b[i]) * (Ax[i] - global_b[i]);
        }

        // вычисление нормы ||b||
        double norm_b = 0.0;
        for (int i = 0; i < N; ++i) {
            norm_b += global_b[i] * global_b[i];
        }

        // проверка сходимости
        double global_norm_Ax_minus_b;
        double global_norm_b;
        MPI_Allreduce(&norm_Ax_minus_b, &global_norm_Ax_minus_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&norm_b, &global_norm_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double relative_error = std::sqrt(global_norm_Ax_minus_b) / std::sqrt(global_norm_b);
        if (relative_error < EPS) {
            converged = true;
        }

        if (rank == 0) {
            std::cout << "Iteration " << iter++ << ": relative error = " << relative_error << std::endl;
        }
    }

    double end_time = MPI_Wtime(); // Конец времени
    double total_time = end_time - start_time;

    // сбор времени работы всех процессов
    double total_time_all;
    MPI_Reduce(&total_time, &total_time_all, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Converged in " << iter << " iterations." << std::endl;
        // std::cout << "Time taken for " << num_procs << " processes: " << total_time_all << " seconds." << std::endl;
        std::cout << "Total time: " << total_time_all << " seconds." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
