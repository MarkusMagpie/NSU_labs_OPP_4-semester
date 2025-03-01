import pandas as pd
import matplotlib.pyplot as plt

# Настройки отображения
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.figure(figsize=(20, 6))

# 1 график - показать минимальный размер массива для которого целесообразно выполнять рекурсивные вызовы
plt.subplot(1, 2, 1)
task1 = pd.read_csv('results1.csv')
plt.plot(task1['size'], task1['time'], 'b-o', linewidth=2, markersize=8)
plt.xscale('log')
plt.title('Зависимость времени от размера массива', pad=20, fontsize=14)
plt.xlabel('Размер массива (log scale)', fontsize=12)
plt.ylabel('Время (сек)', fontsize=12)
plt.grid(True)

# второй 
plt.subplot(1, 2, 2)
task2 = pd.read_csv('results2.csv')
# считаю эффективность
base_time = task2[task2['threads'] == 1]['time'].values[0]
task2['efficiency'] = (base_time / task2['time']) / task2['threads']

ax = plt.gca()
ax.plot(task2['threads'], task2['time'], 'b-o', label='Время')
ax.set_xlabel('Число потоков', fontsize=12)
ax.set_ylabel('Время (сек)', color='b', fontsize=12)
ax.tick_params(axis='y', colors='b')

ax2 = ax.twinx()
ax2.plot(task2['threads'], task2['efficiency'], 'r--o', label='Эффективность')
ax2.set_ylabel('Эффективность', color='r', fontsize=12)
ax2.tick_params(axis='y', colors='r')

plt.title('Зависимость времени и эффективности от потоков', pad=20, fontsize=14)
plt.grid(True)

# Общие настройки
plt.tight_layout()
plt.savefig('combined_plot.png')
plt.show()