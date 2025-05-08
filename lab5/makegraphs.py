import matplotlib.pyplot as plt

threads = [1, 2, 4, 6, 8, 16]

times_bal = {
    "мин_эфф": [5.02736, 10.031, 13.3501, 16.1545, 20.8414, 28.6161],
    "макс_эфф": [5.02752, 11.291, 38.0425, 81.1965, 139.849, 307.173],
    "сред_эфф": [7.54079, 8.78556, 14.3457, 17.9495, 23.846, 30.1585]
}

times_nobal = {
    "мин_эфф": [5.02714, 10.0291, 17.5646, 24.5479, 33.545, 45.6149],
    "макс_эфф": [5.02754, 15.0535, 46.2164, 93.2236, 156.242, 330.977],
    "сред_эфф": [7.54432, 11.2948, 21.8464, 30.8323, 41.3711, 57.1084]
}

def plot():
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))  # 2 строки, 3 столбца

    ax = axes[0, 0]
    for label, times in times_bal.items():
        ax.plot(threads, times, marker='o', label=label)
    ax.set_title("Время (с балансировкой)")
    ax.set_xlabel("Число процессов")
    ax.set_ylabel("Время, с")
    ax.grid(True)
    ax.legend()

    ax = axes[0, 1]
    for label, times in times_bal.items():
        speedup = [times[0] / t for t in times]
        ax.plot(threads, speedup, marker='o', label=label)
    ax.set_title("Ускорение (с балансировкой)")
    ax.set_xlabel("Число процессов")
    ax.set_ylabel("Speedup")
    ax.grid(True)
    ax.legend()

    ax = axes[0, 2]
    for label, times in times_bal.items():
        speedup = [times[0] / t for t in times]
        efficiency = [s / p for s, p in zip(speedup, threads)]
        ax.plot(threads, efficiency, marker='o', label=label)
    ax.set_title("Эффективность (с балансировкой)")
    ax.set_xlabel("Число процессов")
    ax.set_ylabel("Efficiency")
    ax.grid(True)
    ax.legend()

    ax = axes[1, 0]
    for label, times in times_nobal.items():
        ax.plot(threads, times, marker='o', label=label)
    ax.set_title("Время (без балансировки)")
    ax.set_xlabel("Число процессов")
    ax.set_ylabel("Время, с")
    ax.grid(True)
    ax.legend()

    ax = axes[1, 1]
    for label, times in times_nobal.items():
        speedup = [times[0] / t for t in times]
        ax.plot(threads, speedup, marker='o', label=label)
    ax.set_title("Ускорение (без балансировки)")
    ax.set_xlabel("Число процессов")
    ax.set_ylabel("Speedup")
    ax.grid(True)
    ax.legend()

    ax = axes[1, 2]
    for label, times in times_nobal.items():
        speedup = [times[0] / t for t in times]
        efficiency = [s / p for s, p in zip(speedup, threads)]
        ax.plot(threads, efficiency, marker='o', label=label)
    ax.set_title("Эффективность (без балансировки)")
    ax.set_xlabel("Число процессов")
    ax.set_ylabel("Efficiency")
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    plt.savefig("combined_performance.png")
    plt.show()

if __name__ == "__main__":
    plot()
