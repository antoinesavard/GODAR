import numpy as np

# size factor
sf = 1e3
disks_num_y = 50
disks_num_x = int(disks_num_y)
disks_num = disks_num_y * disks_num_x
radius = np.random.lognormal(6, 0.25, 2500)
height = np.random.lognormal(0, 0.25, 2500) - 0.15
height = np.where(height <= 0, 0.1, height)

lines = sf * np.arange(1, 2 * disks_num_y + 1, 2)
lines = lines.astype(str)
lines = lines.tolist()
with open("src/Y.dat", "w") as f:
    for line in lines:
        i = 0
        while i < disks_num_x:
            f.write(line)
            f.write("\n")
            i += 1

lines = sf * np.arange(1, 2 * disks_num_x + 1, 2)
lines = lines.astype(str)
lines = lines.tolist()
with open("src/X.dat", "w") as f:
    i = 0
    while i < disks_num_y:
        for line in lines:
            f.write(line)
            f.write("\n")
        i += 1

lines = radius.astype(str).tolist()
with open("src/R.dat", "w") as f:
    for line in lines:
        f.write(line)
        f.write("\n")

lines = height.astype(str).tolist()
with open("src/H.dat", "w") as f:
    for line in lines:
        f.write(line)
        f.write("\n")
