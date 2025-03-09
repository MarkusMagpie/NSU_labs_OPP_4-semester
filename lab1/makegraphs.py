import matplotlib.pyplot as plt

threads = [1, 2, 4, 6, 8, 10, 12]

def parse_results(filename):
    times = []
    with open(filename, 'r') as f:
        for line in f:
            time = float(line.strip())
            times.append(time)
    return times

def main():
    results_file1 = "results_mpi.txt"

    times1 = parse_results(results_file1)

    if len(times1) != len(threads):
        print("неправильно")
        return

    # вычисляю speedup и эффективность
    speedups1 = [times1[0] / t for t in times1]
    efficiencies1 = [s / t for s, t in zip(speedups1, threads)]

    # простроение графиков
    plt.figure(figsize=(16, 6))

    plt.subplot(1, 3, 1)
    plt.plot(threads, times1, marker='o', linestyle='-', color='b', label="Распараллеленная программа")
    plt.xlabel("Число потоков")
    plt.ylabel("Время выполнения")
    plt.title("Время выполнения")
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(threads, speedups1, marker='o', linestyle='-', color='b', label="Распараллеленная программа")
    plt.xlabel("Число потоков")
    plt.ylabel("Ускорение")
    plt.title("Ускорение")
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(threads, efficiencies1, marker='o', linestyle='-', color='b', label="Распараллеленная программа")
    plt.xlabel("Число потоков")
    plt.ylabel("Эффективность")
    plt.title("Эффективность")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig('combined_plot.png')
    
    plt.show()

if __name__ == "__main__":
    main()
