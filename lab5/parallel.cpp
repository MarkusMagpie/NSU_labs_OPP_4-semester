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
    TASK_NUM = 500, // число задач на каждом процессе
    ITERATIONS = 10 // число итераций выполнения программы
};

// коды статусов для обмена сообщениями
enum Status {
    FINISHED = -1, // все задачи завершены
    NO_TASKS = -2, // нет задач для передачи
    BALANCE = 1, // ФЛАГ БАЛАНСА - делаю с без него и с ним графики 
};

// потоки: для обработки сообщений и выполнения задач
std::thread messageThread;
std::thread executingThread;

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
        std::unique_lock<std::mutex> lock(mtx);

        // ждем, пока в очереди что-то появится или пока очередь закрывается
        cv.wait(lock, [&]{ return !queue.empty() || !running; });
        
        if (queue.empty()) return false;
        n = queue.front(); // читаю первый элемент из queue
        queue.erase(queue.begin()); // удаляю его
        return true;
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

// 1 - МИН ЭФФЕКТИВНОСТЬ - заполниние очереди задачами: время задачи зависит от rank процесса и номера итерации
void refillTaskList(SafeQueue &queue, int size, int rank, int iteration) {
    // задержка задачи: чем дальше rank от (iteration % size), тем дольше
    // *10 - множитель чтобы длинее работало (на мое усмотрение)
    int duration = (abs(rank - (iteration % size)) + 1) * 10;
    for (int i = 0; i < TASK_NUM / size; i++) {
        queue.push(duration);
    }
}

// 2 - МАКС ЭФФЕКТИВНОСТЬ - нулевой процесс имеет задачи, а остальные нет 
void refillTaskList2(SafeQueue &queue, int size, int rank, int iteration) {
    int duration = (abs(rank - (iteration % size)) + 1) * 10;
    if (rank == 0) {
        for (int i = 0; i < TASK_NUM; i++) {
            queue.push(duration);
        }
    }
}

// 3 - СРЕДНЯЯ ЭФФЕКТИВНОСТЬ - первая половина процессов получает по TASK_NUM/2 задач. другая по TASK_NUM 
// void refillTaskList3(SafeQueue &queue, int size, int rank, int iteration) {
//     int duration = (abs(rank - (iteration % size)) + 1) * 10;
//     int tasks = (rank < size / 2) ? TASK_NUM/size/2 : (TASK_NUM/size + TASK_NUM/size/2);
//     for (int i = 0; i < tasks; i++) {
//         queue.push(duration);
//     }
// }

void refillTaskList3(SafeQueue &queue, int size, int rank, int iteration) {
    int duration = (abs(rank - (iteration % size)) + 1) * 10;
    int baseTasks = TASK_NUM / size;
    int extraTasks = TASK_NUM % size;
    int myTasks = baseTasks + (rank < extraTasks ? 1 : 0);

    if (rank < size / 2) {
        myTasks = myTasks / 2;
    } else {
        myTasks = myTasks + myTasks / 2;
    }

    for (int i = 0; i < myTasks; i++) {
        queue.push(duration);
    }
}

// void refillTotal(SafeQueue &queue, int size, int rank, int iteration) {
//     int duration = (abs(rank - (iteration % size)) + 1) * 10;
//     int TOTAL_TASKS = TASK_NUM;
//     int base = TOTAL_TASKS / size;
//     int extra = TOTAL_TASKS % size;
//     int myTasks = base + (rank < extra ? 1 : 0);
//     for (int i = 0; i < myTasks; ++i) {
//         queue.push(duration);
//     }
// }

// выполняем задачи из очереди
void executeTasks(SafeQueue &queue) {
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
        int msg = 0; // принимаемая метка или номер rank
        int recvRank = 0; // кому отправляем задачи
        std::vector<int> taskList; // временный список задач для отправки

        while (queue.isRunning()) {
            /* неблокирующий БЕСКОНЕЧНЫЙ приём сообщения от любого источника с tag = size+1 - УНИКАЛЬНЫЙ тег именно для 
            ообщений-запросов задач: ранг или FINISHED 
            */ 
            MPI_Request request;
            MPI_Irecv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, size + 1, MPI_COMM_WORLD, &request);
            MPI_Status status;
            int flag;
            // ждем пока придёт сообщение ИЛИ очередь остановится(running=false)
            do {
                // опрос завершения Irecv. flag=1 - завершен, flag=0 - не завершен
                MPI_Test(&request, &flag, &status);
            } while (queue.isRunning() && !flag);

            if (!queue.isRunning()) return;

            // если сообщение от executionThread - сигнал FINISHED, то останавливаем работу потока
            if (msg == FINISHED) {
                std::cout << "mT[" << rank << "]: received FINISHED" << std::endl;
                queue.setRunning(false);
                continue;
            }

            // иначе msg = rank процесса, запрашивающего задачи и нам нужно отправить в msg свободные задачи
            recvRank = msg;

            // std::cout << "mT[" << rank << "]: has " << queue.getSize() << " tasks before sending" << std::endl;
            
            int taskCount = 0;
            // если в queue есть избыточные задачи, соберём половину для отправки
            int available = queue.getSize();
            int give = available > 1 ? available / 2 : 0;
            if (give > 0) {
                // вырезаем give задач и отправляем
                taskList.resize(give);
                while (queue.pop(taskList[taskCount])) {
                    taskCount++;
                    if (taskCount == give) {
                        break;
                    }
                }
            }

            std::cout << "mT[" << rank << "]: can send "
                << taskCount << " tasks to proc " << recvRank << std::endl;

            // если так и не собралось ни одной задачи, отправляем статус NO_TASKS; иначе - число свободных задач
            taskCount = (taskCount == 0 ? NO_TASKS : taskCount);
            MPI_Send(&taskCount, 1, MPI_INT, recvRank, rank, MPI_COMM_WORLD);

            // а вот если есть задачи, то после отправки их количества (taskCount), отправляем их самих массивом
            if (taskCount != NO_TASKS) {
                MPI_Send(taskList.data(), taskCount, MPI_INT, recvRank, rank, MPI_COMM_WORLD);
            }
            taskList.clear();

            // std::cout << "mT[" << rank << "]: has " << queue.getSize() << " tasks left" << std::endl;
        }
    });
}

// поток для выполнения задач и инициирования балансировки
void runExecutingThread(SafeQueue &queue, int size, int rank) {
    executingThread = std::thread([&queue, size, rank]() {
        int taskGiven = 0; // сколько задач дал другой процесс моему (rank) процессу

        // повторяем ITERATIONS раз: заполнение, выполнение, балансировка, барьер
        for (int i = 0; i < ITERATIONS; i++) {
            refillTaskList(queue, size, rank, i);
            std::cout << "eT[" << rank << "]: initial queue size = " << queue.getSize() << " tasks" << std::endl;
            executeTasks(queue);

            if (BALANCE) {
                // обход всех процессов по size для запроса задач у них (у себя не спрашивать!)
                for (int j = 0; j < size; j++) {
                    if (j == rank) continue;
                    // отсылка ранга процессу j по особому тегу size+1
                    MPI_Send(&rank, 1, MPI_INT, j, size + 1, MPI_COMM_WORLD);
                    MPI_Status status1;
                    // получение количества избыточных задач от процесса j 
                    MPI_Recv(&taskGiven, 1, MPI_INT, j, j, MPI_COMM_WORLD, &status1);

                    // если процесс дал задачи, принимаем массив и выполняем их
                    if (taskGiven != NO_TASKS) {
                        std::vector<int> receivedTasks(taskGiven);
                        MPI_Status status2;
                        // от процесса j теперь получаю непосредственно вектор интов с избыточными задачами
                        MPI_Recv(receivedTasks.data(), taskGiven, MPI_INT, j, j, MPI_COMM_WORLD, &status2);
                        std::cout << "eT[" << rank << "]: process " << rank << " got " 
                            << taskGiven << " tasks from proc " << j << std::endl;
                        for (int &val : receivedTasks) {
                            queue.push(val);
                        }
                        executeTasks(queue);
                    }
                }
            }
            // синхронизация всех процессов перед следующей итерацией
            MPI_Barrier(MPI_COMM_WORLD);
        }

        // после всех итераций посылаем FINISHED всем остальным процессам
        std::cout << "eT[" << rank << "]: all iterations done, sending FINISHED..." << std::endl;
        int execFinished = FINISHED;
        for (int i = 0; i < size; i++) {
            if (i == rank) continue;
            MPI_Send(&execFinished, 1, MPI_INT, i, size + 1, MPI_COMM_WORLD);
        }
        
        // если единственный процесс, останавливаем очередь
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
    
    runMessageThread(queue, size, rank); // обертки для создания потоков
    runExecutingThread(queue, size, rank);
    messageThread.join();
    executingThread.join();

    // ищем максимальное время выполнения среди всех процессов
    double finish = MPI_Wtime();
    double time  = finish - start;
    double maxTime = 0;
    MPI_Allreduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

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
