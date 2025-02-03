import matplotlib.pyplot as plt
import numpy as np

# Значения (замените на свои данные)
n_processes = np.array([1, 2, 4, 8, ])
T = np.array([])
speedup = T[0] / T   # S(n) = T(1)/T(n)
efficiency = speedup / n_processes  # E(n) = S(n)/n

# график времени работы
plt.figure(figsize=(10, 6))
plt.subplot(3, 1, 1)
plt.plot(n_processes, T, 'o-', label='Time (s)')
plt.xlabel('Number of Processes')
plt.ylabel('Time (s)')
plt.title('Time vs Number of Processes')
plt.grid(True)
plt.legend()

# график ускорения
plt.subplot(3, 1, 2)
plt.plot(n_processes, speedup, 's-', color='green', label='Speedup')
plt.xlabel('Number of Processes')
plt.ylabel('Speedup')
plt.title('Speedup vs Number of Processes')
plt.grid(True)
plt.legend()

# график эффективности
plt.subplot(3, 1, 3)
plt.plot(n_processes, efficiency, 'd-', color='red', label='Efficiency')
plt.xlabel('Number of Processes')
plt.ylabel('Efficiency')
plt.title('Efficiency vs Number of Processes')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
