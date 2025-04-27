#include <iostream>
#include <mpi.h>
#include <thread>
#include <mutex>
#include <vector>
#include <condition_variable>
#include <unistd.h>
#include <cmath>
#include <chrono>

enum Measures {
    TASK_NUM = 50, // число задач на каждом процессе
    ITERATIONS = 10 // число итераций выполнения 
};

// коды статусов для обмена сообщениями
enum Status {
    FINISHED = -1, // все задачи завершены
    NO_TASKS = -2, // нет задач для передачи
    BALANCE = 1 // флаг балансировки включён
};

// потоки: для обработки сообщений и выполнения задач
std::thread messageThread;
std::thread executingThread;

// счётчик завершённых процессов
int finishCount = 0;

// потокобезопасная очередь задач
class SafeQueue {
private:
    std::vector<int> queue; // контейнер для задач (хранит время выполнения задачи)
    std::mutex mtx; // мьютекс для синхронизации доступа к очереди
    std::condition_variable cv; // условная переменная - уведомляет потоки о появлегии новых задач
    bool running = false; // флаг - работает ли очередь?
public:
    // push - добавит задачу в очередь
    void push(int n) {
        {
            std::lock_guard<std::mutex> lock(mtx); // захват мьютекса и освобождение при выходе из зоны видимости
            queue.push_back(n);
        }
        cv.notify_one(); // уведомить один ожидающий поток - что очередь изменилась и может продолжить работу
    }

    // извлекает из начала очереди задачку - сколько милисек нужно спать
    bool pop(int &n) {
        std::lock_guard<std::mutex> lock(mtx);
        if (!queue.empty()) {
            // front возвращает ссылку на первый элемент, а begin возвращает итератор который нужен параметру в erase
            n = queue.front(); // читаю первый элемент из queue
            queue.erase(queue.begin()); // удаляю его
            return true;
        }
        return false; // очередь пуста
    }

    int getSize() {
        std::lock_guard<std::mutex> lock(mtx);
        return queue.size();
    }

    bool isRunning() {
        std::lock_guard<std::mutex> lock(mtx);
        return running;
    }

    void setRunning(bool r) {
        std::lock_guard<std::mutex> lock(mtx);
        running = r;
    }

    bool isEmpty() {
        std::lock_guard<std::mutex> lock(mtx);
        return queue.empty();
    }
};

// заполниние очереди задачами: время задачи зависит от rank процесса и номера итерации
void refillTaskList(SafeQueue &queue, int size, int rank, int iteration) {
    // задержка задачи: чем дальше rank от (iteration % size), тем дольше
    // +1 - чтобы не было нуля
    // *10 - множитель чтобы длинее работало (на мое усмотрение)
    int duration = (abs(rank - (iteration % size)) + 1) * 10;
    for (int i = 0; i < TASK_NUM; i++) {
        queue.push(duration);
    }
}

// выполняем задачи из очереди
void executeTasks(SafeQueue &queue, int rank) {
    while (!queue.isEmpty()) {
        int n;
        if (queue.pop(n)) {
            // симулируем выполнение - спим n миллисекунд
            std::this_thread::sleep_for(std::chrono::milliseconds(n));
        }
    }
}

// поток для обмена сообщениями по MPI: запрос/ответ для балансировки и сигнал завершения
void runMessageThread(SafeQueue &queue, int size, int rank) {
    messageThread = std::thread([&queue, size, rank]() {
        int msg = 0;               // принимаемая метка или номер rank
        int recvRank = 0;          // кому отправляем задачи
        std::vector<int> taskList;      // временный список задач для отправки

        while (queue.isRunning()) {
            // Если получили сообщение о завершении от всех процессов, выходим
            if (finishCount == size) {
                queue.setRunning(false);
                return;
            }

            // Неблокирующий приём сообщения от любого источника с tag = size+1
            MPI_Request request;
            MPI_Irecv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, size + 1,
                      MPI_COMM_WORLD, &request);
            MPI_Status status;
            int flag;
            // Ожидаем, пока придёт сообщение или очередь остановится
            do {
                MPI_Test(&request, &flag, &status);
            } while (queue.isRunning() && !flag);

            if (!queue.isRunning()) return;

            // Если сообщение - сигнал FINISHED, увеличиваем счётчик и останавливаем работу
            if (msg == FINISHED) {
                finishCount++;
                queue.setRunning(false);
                return;
            }

            // Иначе msg = rank процесса, запрашивающего задачи
            recvRank = msg;
            int taskCount = 0;

            // Если в очереди есть избыточные задачи, соберём половину для отправки
            if (queue.getSize() > TASK_NUM / size) {
                taskList.resize(TASK_NUM / size);
                while (queue.pop(taskList[taskCount])) {
                    taskCount++;
                    if (taskCount == TASK_NUM / size) break;
                }
            }

            std::cout << "mT: proc " << rank << " can send "
                 << taskCount << " tasks to proc " << recvRank << std::endl;

            // Если не собралось ни одной задачи, отправляем NO_TASKS
            taskCount = (taskCount == 0 ? NO_TASKS : taskCount);
            MPI_Send(&taskCount, 1, MPI_INT, recvRank, rank, MPI_COMM_WORLD);

            // Если есть задачи, отправляем их массивом
            if (taskCount != NO_TASKS) {
                MPI_Send(taskList.data(), taskCount, MPI_INT,
                         recvRank, rank, MPI_COMM_WORLD);
            }
            taskList.clear();
        }
    });
}

// поток для выполнения задач и инициирования балансировки
void runExecutingThread(SafeQueue &queue, int size, int rank) {
    executingThread = std::thread([&queue, size, rank]() {
        int taskGiven = 0;

        // Повторяем ITERATIONS раз: заполнение, выполнение, балансировка, барьер
        for (int i = 0; i < ITERATIONS; i++) {
            refillTaskList(queue, size, rank, i);
            executeTasks(queue, rank);

            if (BALANCE) {
                // Обход всех процессов для запроса задач у них
                for (int j = 0; j < size; j++) {
                    if (j == rank) continue;
                    MPI_Send(&rank, 1, MPI_INT, j, size + 1, MPI_COMM_WORLD);
                    MPI_Status status1;
                    MPI_Recv(&taskGiven, 1, MPI_INT, j, j,
                             MPI_COMM_WORLD, &status1);

                    std::cout << "eT: proc " << j << " gave "
                         << taskGiven << " tasks to proc " << rank << std::endl;

                    // Если процесс дал задачи, принимаем массив и выполняем их
                    if (taskGiven != NO_TASKS) {
                        std::vector<int> receivedTasks(taskGiven);
                        MPI_Status status2;
                        MPI_Recv(receivedTasks.data(), taskGiven,
                                 MPI_INT, j, j, MPI_COMM_WORLD, &status2);
                        std::cout << "eT: proc " << rank << " got "
                             << taskGiven << " tasks from proc " << j << std::endl;
                        for (int &val : receivedTasks) queue.push(val);
                        executeTasks(queue, rank);
                    }
                }
            }
            // синхронизация всех процессов перед следующей итерацией
            MPI_Barrier(MPI_COMM_WORLD);
        }

        // После всех итераций посылаем FINISHED всем остальным
        int execFinished = FINISHED;
        for (int i = 0; i < size; i++) {
            if (i == rank) continue;
            MPI_Send(&execFinished, 1, MPI_INT, i, size + 1, MPI_COMM_WORLD);
        }
        // Если единственный процесс, останавливаем очередь
        if (size == 1) queue.setRunning(false);
    });
}

int main(int argc, char *argv[]) {
    // инициализируем MPI с поддержкой многопоточности (MPI_Init_thread с уровнем MPI_THREAD_SINGLE = MPI_Init)
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided); // одновременная работа нескольких потоков использующих MPI-вызовы

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::cout << "size: " << size
             << ", tasks for each process: " << TASK_NUM
             << ", iterations: " << ITERATIONS << std::endl;
    }

    // создание потокобезопасной очереди и запускаем потоки учтанавливая флаг
    SafeQueue queue;
    queue.setRunning(true);

    double start = MPI_Wtime();
    runMessageThread(queue, size, rank);
    runExecutingThread(queue, size, rank);

    // ожидаем завершения двух потоков
    if (messageThread.joinable()) messageThread.join();
    if (executingThread.joinable()) executingThread.join();

    // ищем максимальное время выполнения среди всех процессов
    double finish = MPI_Wtime();
    double time  = finish - start;
    double maxTime = 0;
    MPI_Allreduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "total time: " << maxTime << " sec" << std::endl;
    }

    MPI_Finalize();
    return 0;
}

/* 
компилируй так: mpic++ -O3 -o parallel parallel.cpp
                mpirun  -np 4 ./parallel

дебаг:          mpic++ -O3 -o parallel parallel.cpp -g
                mpirun -np 4 xterm -fa 'Monospace' -fs 14 -e gdb -ex run parallel
*/ 
