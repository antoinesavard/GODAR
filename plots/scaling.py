import numpy as np
import matplotlib.pyplot as plt

sim_time200 = np.array([133.2, 28, 14, 11, 9])
sim_time500 = np.array([2163, 451.2, 226.8, 162, 127.8])
proc_num = np.array([1, 5, 10, 15, 20])

speedup200 = sim_time200[0] / sim_time200
speedup500 = sim_time500[0] / sim_time500

serial200 = np.nanmean(
    (1 - speedup200 / proc_num) / (speedup200 - speedup200 / proc_num)
)
serial500 = np.nanmean(
    (1 - speedup500 / proc_num) / (speedup500 - speedup500 / proc_num)
)

print(serial200, serial500, speedup200, speedup500)

plt.plot(proc_num, proc_num)
plt.plot(proc_num, speedup200, marker=".")
plt.plot(proc_num, speedup500, marker=".")

plt.show()
