import matplotlib.pyplot as plt

threads = [1, 2, 4, 8, 16]

times_bal = {
    "min": [50.267, 37.7064, 28.7021, 24.1206, 20.4622], 
    "max": [50.2649, 44.3969, 30.7238, 19.3081, 12.1505],
    "custom": [50.967, 37.7001, 29.8236, 25.8473, 23.3321],
}

times_nobal = {
    "min": [50.2694, 50.1387, 43.8937, 41.6695, 36.5532],
    # "max": [50.2684, 75.2653, 115.695, 195.37, 276.212],
    "custom": [50.5604, 46.4544, 44.449, 42.2336, 39.6562 ],
}

def plot():
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))  # 2 строки, 3 столбца

    ax = axes[0, 0]
    for label, times in times_bal.items():
        ax.plot(threads, times, marker='o', label=label)
    ax.set_title("время (с балансировкой)")
    ax.set_xlabel("число процессов")
    ax.set_ylabel("время, с")
    ax.grid(True)
    ax.legend()

    ax = axes[0, 1]
    for label, times in times_bal.items():
        acceleration = [times[0] / t for t in times]
        ax.plot(threads, acceleration, marker='o', label=label)
    ax.set_title("ускорение (с балансировкой)")
    ax.set_xlabel("число процессов")
    ax.set_ylabel("acceleration")
    ax.grid(True)
    ax.legend()

    ax = axes[0, 2]
    for label, times in times_bal.items():
        acceleration = [times[0] / t for t in times]
        efficiency = [s / p for s, p in zip(acceleration, threads)]
        ax.plot(threads, efficiency, marker='o', label=label)
    ax.set_title("эффективность (с балансировкой)")
    ax.set_xlabel("число процессов")
    ax.set_ylabel("efficiency")
    ax.grid(True)
    ax.legend()

    ax = axes[1, 0]
    for label, times in times_nobal.items():
        ax.plot(threads, times, marker='o', label=label)
    ax.set_title("время (без балансировки)")
    ax.set_xlabel("число процессов")
    ax.set_ylabel("время, с")
    ax.grid(True)
    ax.legend()

    ax = axes[1, 1]
    for label, times in times_nobal.items():
        acceleration = [times[0] / t for t in times]
        ax.plot(threads, acceleration, marker='o', label=label)
    ax.set_title("ускорение (без балансировки)")
    ax.set_xlabel("число процессов")
    ax.set_ylabel("acceleration")
    ax.grid(True)
    ax.legend()

    ax = axes[1, 2]
    for label, times in times_nobal.items():
        acceleration = [times[0] / t for t in times]
        efficiency = [s / p for s, p in zip(acceleration, threads)]
        ax.plot(threads, efficiency, marker='o', label=label)
    ax.set_title("эффективность (без балансировки)")
    ax.set_xlabel("число процессов")
    ax.set_ylabel("efficiency")
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    plt.savefig("combined_performance.png")
    plt.show()

if __name__ == "__main__":
    plot()
