import matplotlib.pyplot as plt

threads = [1, 2, 4, 8]

def parse_results(filename):
    times = []
    # with open(filename, 'r') as f:
    #     for line in f:
    #         if "Time taken:" in line:
    #             time = float(line.split("Time taken:")[1].split()[0])  # первый сепаратор - тайм тейкен, второй - пробел
    #             times.append(time)
    with open(filename, 'r') as f:
        for line in f:
            time = float(line.strip())
            times.append(time)
    return times

# Основная функция
def main():
    results_file1 = "results1.txt"
    results_file2 = "results2.txt"

    times1 = parse_results(results_file1)
    times2 = parse_results(results_file2)

    if len(times1) != len(threads):
        print("хуйня")
        return

    # вычисляю speedup и эффективность
    speedups1 = [times1[0] / t for t in times1]
    efficiencies1 = [s / t for s, t in zip(speedups1, threads)]

    speedups2 = [times2[0] / t for t in times2]
    efficiencies2 = [s / t for s, t in zip(speedups2, threads)]




    # простроение двух графиков
    plt.figure(figsize=(16, 6))

    plt.subplot(1, 3, 1)
    plt.plot(threads, times1, marker='o', linestyle='-', color='y', label="Несколько секций")
    # plt.plot(threads, times2, marker='s', linestyle='--', color='c', label="Одна секция")
    plt.xlabel("Число потоков")
    plt.ylabel("Время выполнения")
    plt.title("Графики времени выполнения распараллеленных программ")
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(threads, speedups1, marker='o', linestyle='-', color='b', label="Несколько секций")
    # plt.plot(threads, speedups2, marker='s', linestyle='--', color='g', label="Одна секция")
    plt.xlabel("Число потоков")
    plt.ylabel("Ускорение")
    plt.title("Графики ускорения распараллеленных программ")
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(threads, efficiencies1, marker='o', linestyle='-', color='r', label="Несколько секций")
    # plt.plot(threads, efficiencies2, marker='s', linestyle='--', color='m', label="Одна секция")
    plt.xlabel("Число потоков")
    plt.ylabel("Эффективность")
    plt.title("Графики эффективности распараллеленных программ")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
