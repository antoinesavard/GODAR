import numpy as np

# size factor
sf = 1e3
disks_num_y = 10
disks_num_x = int(disks_num_y / 2)
disks_num = disks_num_y * disks_num_x
radius = 1e3 / 2 * np.ones(int(disks_num))
height = np.random.normal(loc=3.3, scale=1.8, size=int(disks_num))
height = np.where(height <= 0, 0.1, height)

lines = sf * np.arange(0.5, disks_num_y + 0.5)
lines = lines.astype(str)
lines = lines.tolist()
with open("src/Y.dat", "w") as f:
    for line in lines:
        i = 0
        while i < disks_num_x:
            f.write(line)
            f.write("\n")
            i += 1

lines = sf * np.arange(1 / 2 * disks_num_x + 0.5, 3 / 2 * disks_num_x + 0.5)
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
