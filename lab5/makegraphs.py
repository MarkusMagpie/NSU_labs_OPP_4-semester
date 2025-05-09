import matplotlib.pyplot as plt

threads = [1, 2, 4, 8, 16]

# times_bal = {
#     "ref1": [28.6161, 20.8414, 13.3501, 10.031, 5.02736],
#     "ref2": [365.433, 139.849, 38.0425, 11.291, 5.02752],
#     "ref3": [38.8176, 23.846, 14.3457, 8.78556, 7.54079]
# }

# times_nobal = {
#     "ref1": [52.1396, 33.545, 17.5646, 10.0291, 5.02714],
#     "ref2": [386.392, 156.242, 46.2164, 15.0535, 5.02754],
#     "ref3": [71.4146, 41.3711,  21.8464, 11.2948, 7.54432]
# }

times_bal = {
    "ref1": [50.2722, 37.7064, 28.6344, 24.3262, 20.669], 
    "ref3": [50.967, 37.7001, 29.8236, 25.8473, 22.9981],
}

times_nobal = {
    "ref1": [50.2694, 50.1387, 43.8937, 41.6695, 36.5532],
    "ref3": [50.5604, 46.4544, 44.449, 41.2336, 39.6562 ],
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
