#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

const double EPS = 1e-5; // епсилон - критерий завершения счета
const double TAU = 0.01; // коэффициент в вычислении шага метода простой итерации

// Функция для умножения подматрицы A на вектор x
std::vector<double> mat_vec_mult(const std::vector<std::vector<double>>& A, const std::vector<double>& x, int start_row) {
    int rows = A.size();
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
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int N = 1000; // размерность системы
    int local_N = N / num_procs; // количество строк для каждого процесса

    // инициализация матрицы A и b (каждый процесс хранит только свою часть)
    std::vector<std::vector<double>> A(local_N, std::vector<double>(N, 1.0));
    std::vector<double> b(local_N, N + 1);
    std::vector<double> x(N, 0.0); // начальное приближение x

    // Запуск итераций
    bool converged = false;
    int iter = 0;
    while (!converged) {
        // Вычисляем A * x (локально)
        std::vector<double> Ax_local = mat_vec_mult(A, x, rank * local_N);

        // Собираем Ax из всех процессов
        std::vector<double> Ax(N, 0.0);
        MPI_Allgather(Ax_local.data(), local_N, MPI_DOUBLE, Ax.data(), local_N, MPI_DOUBLE, MPI_COMM_WORLD);

        // Обновляем x
        double error = 0.0;
        for (int i = 0; i < N; ++i) {
            double new_x = x[i] - TAU * (Ax[i] - b[i]);
            error += (new_x - x[i]) * (new_x - x[i]);
            x[i] = new_x;
        }

        // Проверяем сходимость
        double global_error;
        MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (std::sqrt(global_error) < EPS) {
            converged = true;
        }

        if (rank == 0) {
            std::cout << "Iteration " << iter++ << ", error: " << std::sqrt(global_error) << std::endl;
        }
    }

    if (rank == 0) {
        std::cout << "Converged in " << iter << " iterations." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
