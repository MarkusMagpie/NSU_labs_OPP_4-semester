import matplotlib.pyplot as plt

# Число потоков, с которыми запускались тесты
threads = [1, 2, 4, 8]

# Считываем времена работы из файла
def parse_results(filename):
    times = []
    with open(filename, 'r') as f:
        for line in f:
            if "Time taken:" in line:
                time = float(line.split("Time taken:")[1].split()[0])  # Извлекаем время
                times.append(time)
    return times

# Основная функция
def main():
    results_file = "results.txt"
    times = parse_results(results_file)

    if len(times) != len(threads):
        print("Хуйня")
        return

    # вычисляю speedup и эффективность
    t1 = times[0]  # время на первом потоке
    speedups = [t1 / t for t in times]
    efficiencies = [s / t for s, t in zip(speedups, threads)]

    # Построение графиков
    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    plt.plot(threads, speedups, marker='o', linestyle='-', color='b', label="Ускорение")
    plt.xlabel("Число потоков")
    plt.ylabel("Ускорение")
    plt.title("График ускорения")
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(threads, efficiencies, marker='o', linestyle='-', color='r', label="Эффективность")
    plt.xlabel("Число потоков")
    plt.ylabel("Эффективность")
    plt.title("График эффективности")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
