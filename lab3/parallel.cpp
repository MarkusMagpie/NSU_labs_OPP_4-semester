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
    std::uniform_real_distribution<float> dis(0.0f, 1000.0f); // настройка равномерного распределения 

    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            matrix[i * n2 + j] = dis(gen); // генератор gen создает случайное число, dis преобразует его в float в диапазоне [0, 10000)
        }
    }
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

void print_matrix(float* matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << matrix[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv); // создаем предопределеннную область связи, содержащую все процессы MPI программы, с ней связывается коммуникатор MPI_COMM_WORLD

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // номер (ранг) процесса, вызвавшего функцию
    MPI_Comm_size(MPI_COMM_WORLD, &size); // количчество процессов в области связи коммуникатора MPI_COMM_WORLD 

    // 1 - автоматическое определение размеров решетки
    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);
    int p1 = dims[0], p2 = dims[1]; // параметры решетки 
    
    // размеры матриц (должны быть кратны p1 и p2)
    int n1 = 4, n2 = 4, n3 = 4;

    if (n1 % p1 != 0 || n3 % p2 != 0) {
        if (rank == 0) {
            std::cout << "размеры матриц должны быть кратны числу процессов:" << std::endl;
            std::cout << "n1 = " << n1 << ", p1 = " << p1 << std::endl;
            std::cout << "n3 = " << n3 << ", p2 = " << p2 << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    // 2 - создание декартовой решетки - коммуникатора для 2D 
    int periods[2] = {0, 0};
    MPI_Comm comm_grid;
    // создание нового коммуникатора comm_grid с заданной декартовой топологией 
    // (все процессы из исходного коммуникатора MPI_COMM_WORLD перестраиваются в двумерную сетку)
    // https://www.opennet.ru/docs/RUS/mpi-1/node129.html
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_grid);

    // 3 - получение координат текущего процесса
    int coords[2];
    // определение декартовых координат процесса по его рангу rank в группе
    MPI_Cart_coords(comm_grid, rank, 2, coords);
    int x = coords[0], y = coords[1];

    // 4 - создание подкоммуникаторов 
    MPI_Comm comm_row, comm_col;
    int remain_dims_row[2] = {false, true}; // изменения по x не сохраняются, по y сохраняются
    int remain_dims_col[2] = {true, false};
    MPI_Cart_sub(comm_grid, remain_dims_row, &comm_row);
    MPI_Cart_sub(comm_grid, remain_dims_col, &comm_col);

    // 5 - подматрицы
    int local_n1 = n1 / p1;
    int local_n3 = n3 / p2;

    // std::vector<float> A, B;
    // if (rank == 0) {
    //     A.resize(n1 * n2, 1.0);
    //     B.resize(n2 * n3, 1.0);
    // } 

    // 6 - инициализация данных на нулеввых процессах подкоммуникаторов
    std::vector<float> A, B;
    int row_rank, col_rank;
    MPI_Comm_rank(comm_row, &row_rank); // получаю ранг процесса, вызвавшего функцию
    MPI_Comm_rank(comm_col, &col_rank);

    if (row_rank == 0) {                // если процесс в подкоммуникаторе строчных процессов равен нулю, то делаю resize
        A.resize(n1 * n2, 1.0);
        fill_matrix(A.data(), n1, n2);
    }
    if (col_rank == 0) {
        B.resize(n2 * n3, 1.0);
        fill_matrix(B.data(), n2, n3);
    }    

    if (row_rank == 0) {
        std::cout << "матрица А:" << std::endl;
        print_matrix(A.data(), n1, n2);
    }
    std::cout << std::endl;
    if (col_rank == 0) {
        std::cout << "матрица B:" << std::endl;
        print_matrix(B.data(), n2, n3);
    }

    // 7 - распределить матрицу А по горизонтальным полоскам
    // ЛОГИКА: хочу чтобы процессы ВДОЛЬ СТРОК РЕШЕТКИ разобрали блоки матрицы А 
    // ИЗ УСЛОВИЯ: Матрица А распределяется по горизонтальным полосам вдоль координаты (x,0).
    std::vector<float> local_A(local_n1 * n2);
    MPI_Scatter(A.data(), local_n1 * n2, MPI_FLOAT, 
                local_A.data(), local_n1 * n2, MPI_FLOAT, 
                0, comm_row);

    // 8 - распределить матрицу B по вертикальным полоскам
    // ИЗ УСЛОВИЯ: Матрица B распределяется по вертикальным полосам вдоль координаты (0,y). 
    MPI_Datatype column_type;
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
    MPI_Type_commit(&column_type);

    std::vector<float> local_B(n2 * local_n3);
    if (col_rank == 0) {
        MPI_Scatter(B.data(), 1, column_type,
                    local_B.data(), n2 * local_n3, MPI_FLOAT,
                    0, comm_col);
    }

    // 9 - рассылка данных (стадии вычисления 3-4 из условия)
    MPI_Bcast(local_A.data(), local_n1 * n2, MPI_FLOAT, col_rank, comm_col);
    MPI_Bcast(local_B.data(), n2 * local_n3, MPI_FLOAT, row_rank, comm_row);

    // 10 - локальное умножение подматриц
    std::vector<float> local_C(local_n1 * local_n3, 0.0);
    auto start = std::chrono::high_resolution_clock::now();
    mult_matrix(local_A.data(), local_B.data(), local_C.data(), local_n1, n2, local_n3);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // 11 - сборака результатов с учетом смещений
    std::vector<float> C;
    std::vector<int> recvcounts(size, local_n1 * local_n3); // массив количества элементов для каждого процесса в С
    std::vector<int> displs(size);                          // массив смещений для каждого процесса в С (ранг - смещение)
    for (int i = 0; i < p1; ++i) {
        for (int j = 0; j < p2; ++j) {
            int global_rank;
            int current_coords[2] = {i, j};
            MPI_Cart_rank(comm_grid, current_coords, &global_rank); // получаю ранг процесса в коммуникаторе сетки
            displs[global_rank] = (i * local_n1) * n3 + (j * local_n3);
        }
    }

    if (rank == 0) C.resize(n1 * n3);

    MPI_Gatherv(local_C.data(), local_n1 * local_n3, MPI_FLOAT,
                C.data(), recvcounts.data(), displs.data(), MPI_FLOAT,
                0, comm_grid);

    // 12 - вывод 
    if (rank == 0) {
        std::cout << "матрица C:" << std::endl;
        print_matrix(C.data(), n1, n3);
        std::cout << elapsed.count() << std::endl;
    }

    MPI_Type_free(&column_type);
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
