import numpy as np
import matplotlib.pyplot as plt

sim_time200 = np.array([133.2, 28, 14, 11, 9])
sim_time500 = np.array([2163, 451.2, 226.8, 162, 127.8])
sim_time5000 = np.array([11988, 4068, 2093, 994, 520, 341, 194])
proc_num = np.array([1, 5, 10, 15, 20])
proc_num1 = np.array([1, 2, 4, 8, 16, 32, 64])

speedup200 = sim_time200[0] / sim_time200
speedup500 = sim_time500[0] / sim_time500
speedup5000 = sim_time5000[0] / sim_time5000

serial200 = (1 - speedup200 / proc_num) / (speedup200 - speedup200 / proc_num)

serial500 = (1 - speedup500 / proc_num) / (speedup500 - speedup500 / proc_num)


serial5000 = (1 - speedup5000 / proc_num1) / (
    speedup5000 - speedup5000 / proc_num1
)


print(serial200, serial500, serial5000, speedup200, speedup500, speedup5000)

plt.plot(proc_num1, proc_num1)
plt.plot(proc_num1, speedup5000, marker=".")
plt.plot(proc_num, speedup200, marker=".")
plt.plot(proc_num, speedup500, marker=".")

plt.show()
