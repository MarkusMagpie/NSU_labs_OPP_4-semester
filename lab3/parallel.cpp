#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <chrono>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv); // создаем предопределеннную область связи, содержащую все процессы MPI программы, с ней связывается коммуникатор MPI_COMM_WORLD

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // номер (ранг) процесса, вызвавшего функцию
    MPI_Comm_size(MPI_COMM_WORLD, &size); // количчество процессов в области связи коммуникатора MPI_COMM_WORLD 

    // параметры решетки 
    int p1 = 2, p2 = 2;
    // размеры матриц (должны быть кратны p1 и p2)
    int n1 = 1000, n2 = 1000, n3 = 1000;

    // создание декартовой решетки - коммуникатора для 2D 
    int dims[2] = {p1, p2};
    int periods[2] = {0, 0};
    MPI_Comm comm_grid;
    // создание нового коммуникатора comm_grid с заданной декартовой топологией 
    // (все процессы из исходного коммуникатора MPI_COMM_WORLD перестраиваются в двумерную сетку)
    // https://www.opennet.ru/docs/RUS/mpi-1/node129.html
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_grid);

    // получение координат текущего процесса
    int coords[2];
    // определение декартовых координат процесса по его рангу rank в группе
    MPI_Cart_coords(comm_grid, rank, 2, coords);
    int x = coords[0], y = coords[1];

    // подкоммуникаторы 
    /*
    в MPI_Comm_split в разбиении процессов на подкоммуникаторы отвечают 2 параметра:
        color - группировка процессов по логике. То есть процессы, у которых одинаковое значение color, попадают в один подкоммуникатор 
        key - задание рангов процессам уже внутри подкоммуникатора (например - по возрастанию координат)
    */
    MPI_Comm comm_row, comm_col;
    MPI_Comm_split(comm_grid, x, y, &comm_row); // процессы в одной строке
    MPI_Comm_split(comm_grid, y, x, &comm_col); // процессы в одном столбце

    // подматрицы
    int local_n1 = n1 / p1;
    int local_n3 = n3 / p2;

    // std::vector<double> A, B;
    // if (rank == 0) {
    //     A.resize(n1 * n2, 1.0);
    //     B.resize(n2 * n3, 1.0);
    // } 

    std::vector<double> A, B;
    int row_rank, col_rank;
    MPI_Comm_rank(comm_row, &row_rank); // получаю ранг процесса, вызвавшего функцию
    MPI_Comm_rank(comm_col, &col_rank);

    if (row_rank == 0) {                // если процесс в подкоммуникаторе строчных процессов равен нулю, то делаю resize
        A.resize(n1 * n2, 1.0);
    }
    if (col_rank == 0) {
        B.resize(n2 * n3, 1.0);
    }

    // распределить матрицу А по горизонтальным полоскам
    // ЛОГИКА: хочу чтобы процессы ВДОЛЬ СТРОК РЕШЕТКИ разобрали блоки матрицы А 
    // ИЗ УСЛОВИЯ: Матрица А распределяется по горизонтальным полосам вдоль координаты (x,0).
    std::vector<double> local_A(local_n1 * n2);
    MPI_Scatter(A.data(), local_n1 * n2, MPI_DOUBLE, 
                local_A.data(), local_n1 * n2, MPI_DOUBLE, 0, comm_row);

    // по аналогии с матрицей B (только теперь)
    // ИЗ УСЛОВИЯ: Матрица B распределяется по вертикальным полосам вдоль координаты (0,y)
    std::vector<double> local_B(n2 * local_n3);
    MPI_Scatter(B.data(), n2 * local_n3, MPI_DOUBLE, 
                local_B.data(), n2 * local_n3, MPI_DOUBLE, 0, comm_col);


                

    MPI_Finalize();
    return 0;
}

/* 
компилируй так: mpic++ -O3 -o parallel parallel.cpp
                mpirun  -np 4 ./parallel

дебаг:          mpic++ -O3 -o parallel parallel.cpp -g
                mpirun -np 4 xterm -fa 'Monospace' -fs 14 -e gdb -ex run parallel1
*/ 
