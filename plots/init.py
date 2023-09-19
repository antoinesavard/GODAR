# file to create experience 06 initial conditions with 2500 particles

import numpy as np

# size factor
sf = 1e3
disks_num_y = 50
disks_num_x = int(disks_num_y)
disks_num = disks_num_y * disks_num_x
radius = np.random.lognormal(6, 0.25, 2500)
height = np.random.lognormal(0, 0.25, 2500) - 0.15
height = np.where(height <= 0.1, 0.1, height)

lines = sf * np.arange(1, 2 * disks_num_y + 1, 2)
lines = lines.astype(str)
lines = lines.tolist()
with open("input_files/Y.dat", "w") as f:
    for line in lines:
        i = 0
        while i < disks_num_x:
            f.write(line)
            f.write("\n")
            i += 1

lines = sf * np.arange(1, 2 * disks_num_x + 1, 2)
lines = lines.astype(str)
lines = lines.tolist()
with open("input_files/X.dat", "w") as f:
    i = 0
    while i < disks_num_y:
        for line in lines:
            f.write(line)
            f.write("\n")
        i += 1

lines = radius.astype(str).tolist()
with open("input_files/R.dat", "w") as f:
    for line in lines:
        f.write(line)
        f.write("\n")

lines = height.astype(str).tolist()
with open("input_files/H.dat", "w") as f:
    for line in lines:
        f.write(line)
        f.write("\n")
