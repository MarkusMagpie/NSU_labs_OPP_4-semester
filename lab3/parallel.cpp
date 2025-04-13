#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <random>

void fill_matrix(float* matrix, int n1, int n2) {
    std::mt19937 gen(100); // если seed одинаковый, то и последовательность рандомных чисел одинаковая при любом запуске
    std::uniform_real_distribution<float> dis(0.0f, 100.0f); // настройка равномерного распределения 

    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            matrix[i * n2 + j] = dis(gen); // генератор gen создает случайное число, dis преобразует его в float в диапазоне [0, 100)
        }
    }
}

void print_matrix(float* matrix, int n1, int n2) {
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            std::cout << matrix[i * n2 + j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void mult_matrix(float* local_A, float* local_B, float* local_C, int local_n1, int n2, int local_n3) {
    for (int i = 0; i < local_n1; ++i) {
        for (int k = 0; k < n2; ++k) {
            for (int j = 0; j < local_n3; ++j) {
                local_C[i * local_n3 + j] += local_A[i * n2 + k] * local_B[k * local_n3 + j];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv); // создаем предопределеннную область связи, содержащую все процессы MPI программы, с ней связывается коммуникатор MPI_COMM_WORLD

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // номер (ранг) процесса, вызвавшего функцию
    MPI_Comm_size(MPI_COMM_WORLD, &size); // количество процессов в области связи коммуникатора MPI_COMM_WORLD 

    // 1 - автоматическое определение размеров решетки
    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims); // перестраивает все процессы из MPI_COMM_WORLD в двумерную решётку
    int p1 = dims[0], p2 = dims[1]; // параметры решетки 
    
    // размеры матриц
    int n1 = 4, n2 = 4, n3 = 4;

    // 2 - создание декартовой решетки - коммуникатора для 2D 
    int periods[2] = {0, 0};
    MPI_Comm comm_grid;
    // создание нового коммуникатора comm_grid с заданной декартовой топологией 
    // (все процессы из исходного коммуникатора MPI_COMM_WORLD перестраиваются в двумерную сетку)
    // https://www.opennet.ru/docs/RUS/mpi-1/node129.html
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_grid);

    // 3 - перевод ранга в координаты декартовой решетки
    int coords[2];
    // определение декартовых координат процесса по его рангу rank в группе
    MPI_Cart_coords(comm_grid, rank, 2, coords);
    int x = coords[0], y = coords[1];

    // 4 - создание подкоммуникаторов 
    MPI_Comm comm_row, comm_col;
    // в MPI координаты процесса (x,y) <=> (номер строки, номер столбца)
    int remain_dims_row[2] = {false, true}; 
    int remain_dims_col[2] = {true, false};
    MPI_Cart_sub(comm_grid, remain_dims_row, &comm_row); // все процессы в одной строке (с одним x) в одном подкоммуникаторе -> фиксирую x, разрешаю менять y
    MPI_Cart_sub(comm_grid, remain_dims_col, &comm_col);

    // 5 - размер локальных блоков
    int local_n1 = n1 / p1; // число строк в блоке
    int local_n3 = n3 / p2; // число столбцов

    // 6 - инициализация матриц
    std::vector<float> A, B, C;

    if (rank == 0) {
        A.resize(n1 * n2);
        B.resize(n2 * n3);
        C.resize(n1 * n3);
        fill_matrix(A.data(), n1, n2);
        fill_matrix(B.data(), n2, n3);
        std::cout << "Матрица А:" << std::endl;
        print_matrix(A.data(), n1, n2);
        std::cout << "Матрица B:" << std::endl;
        print_matrix(B.data(), n2, n3);
    }

    double start = MPI_Wtime();

    // 7 - распределить матрицу А по горизонтальным полоскам в локальные матрицы partA процессов
    // ИЗ УСЛОВИЯ: Матрица А распределяется по горизонтальным полосам вдоль координаты (x,0)
    std::vector<float> local_A(local_n1 * n2);
    if (y == 0) {
        MPI_Scatter(A.data(), local_n1 * n2, MPI_FLOAT, 
                    local_A.data(), local_n1 * n2, MPI_FLOAT, 
                    0, comm_col);
    }
    // ИЗ УСЛОВИЯ: Полосы А распространяются в измерении y.
    MPI_Bcast(local_A.data(), local_n1 * n2, MPI_FLOAT, 0, comm_row);

    // 8 - распределить матрицу B по вертикальным полоскам в локальные матрицы partB процессов
    // ИЗ УСЛОВИЯ: Матрица B распределяется по вертикальным полосам вдоль координаты (0,y).
    MPI_Datatype column_type, resized_column_type;
    /*
    MPI_Type_vector - определяет новый тип данных, состоящий из указанного числа блоков указанного размера
    параметры: 
        n2 - количество блоков в созданном векторе
        local_n3 - размер каждого блока
        n3 - шаг/количество элементов между началами соседних блоков
        MPI_FLOAT - тип данных каждого блока
        &column_type - указатель на созданный тип
    */
    MPI_Type_vector(n2, local_n3, n3, MPI_FLOAT, &column_type);
    // задаю новую длину типа данных column_type
    MPI_Type_create_resized(column_type, 0, local_n3 * sizeof(float), &resized_column_type);
    MPI_Type_commit(&resized_column_type);

    std::vector<float> local_B(n2 * local_n3);
    if (x == 0) {
        MPI_Scatter(B.data(), 1, resized_column_type,
                    local_B.data(), n2 * local_n3, MPI_FLOAT,
                    0, comm_row);
    }
    // Из УСЛОВИЯ: Полосы B распространяются в измерении х.
    MPI_Bcast(local_B.data(), n2 * local_n3, MPI_FLOAT, 0, comm_col);

    // 9 - локальное умножение подматриц
    // ИЗ УСЛОВИЯ: Каждый процесс вычисляет одну подматрицу произведения.
    std::vector<float> local_C(local_n1 * local_n3, 0.0);
    mult_matrix(local_A.data(), local_B.data(), local_C.data(), local_n1, n2, local_n3);

    // 10 - сбор - объединение локальных подматриц local_C в итоговую матрицу C 
    // block_type - описывает подматрицу размера local_n1 * local_n3, с шагом в памяти равным n3
    MPI_Datatype block_type, resized_block_type;
    MPI_Type_vector(local_n1, local_n3, n3, MPI_FLOAT, &block_type);
    MPI_Type_create_resized(block_type, 0, local_n3 * sizeof(float), &resized_block_type);
    MPI_Type_commit(&resized_block_type);

    std::vector<int> recvcounts(size, 1);       // массив количества элементов, полученных из каждого процесса <=> каждой подматрицы в С
    std::vector<int> displs(size);              // массив смещений для каждой подматрицы в С (ранг - смещение)
    for (int i = 0; i < p1; ++i) {
        for (int j = 0; j < p2; ++j) {
            int global_rank;
            int crd[2] = {i, j};
            MPI_Cart_rank(comm_grid, crd, &global_rank);
            displs[global_rank] = i * p2 * local_n1 + j;
        }
    }

    MPI_Gatherv(local_C.data(), local_n1 * local_n3, MPI_FLOAT,
               C.data(), recvcounts.data(), displs.data(), resized_block_type, 
               0, comm_grid);

    double end = MPI_Wtime();
    double local_elapsed = end - start;

    // сбор времени выполнения со всех процессов
    std::vector<double> all_elapsed(size);
    MPI_Gather(&local_elapsed, 1, MPI_DOUBLE, 
              all_elapsed.data(), 1, MPI_DOUBLE, 
              0, comm_grid);

    // 11 - вывод 
    if (rank == 0) {
        std::cout << "матрица C после сборки подматриц:" << std::endl;
        print_matrix(C.data(), n1, n3);
        double avg_elapsed = std::accumulate(all_elapsed.begin(), all_elapsed.end(), 0.0) / size;
        std::cout << avg_elapsed << std::endl;
    }

    MPI_Type_free(&column_type);
    MPI_Type_free(&resized_column_type);

    MPI_Type_free(&block_type);
    MPI_Type_free(&resized_block_type);

    MPI_Comm_free(&comm_row);
    MPI_Comm_free(&comm_col);
    MPI_Comm_free(&comm_grid);

    MPI_Finalize();
    return 0;
}

/* 
компилируй так: mpic++ -O3 -o parallel parallel.cpp
                mpirun  -np 4 ./parallel

дебаг:          mpic++ -O3 -o parallel parallel.cpp -g
                mpirun -np 4 xterm -fa 'Monospace' -fs 14 -e gdb -ex run parallel
*/ 
