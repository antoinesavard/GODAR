# file to create experience 07 initial conditions with 5000 particles

import numpy as np

num = 5000
height = np.random.lognormal(0, 0.25, 5000) - 0.15
height = np.where(height <= 0.1, 0.1, height)

# x
with open("output/x.06", "r") as f:
    xlines1 = f.read().split()
    xlines1 = np.asarray(xlines1).astype(float).reshape(1000, 2500)
    xlines2 = xlines1[-1] + 20000
    xlines = np.concatenate((xlines1[-1], xlines2))

with open("input_files/X1.dat", "w") as f:
    for i, xline in enumerate(xlines):
        f.write(str(xline))
        f.write("\n")

# y
with open("output/y.06", "r") as f:
    ylines1 = f.read().split()
    ylines1 = np.asarray(ylines1).astype(float).reshape(1000, 2500)
    ylines2 = ylines1[-1] + 20000
    ylines = np.concatenate((ylines1[-1], ylines2))

with open("input_files/Y1.dat", "w") as f:
    for i, yline in enumerate(ylines):
        f.write(str(yline))
        f.write("\n")

# r
with open("output/r.06", "r") as f:
    rlines1 = f.read().split()
    rlines1 = np.asarray(rlines1).astype(float).reshape(1000, 2500)
    rlines2 = rlines1[-1] + 20000
    rlines = np.concatenate((rlines1[-1], rlines2))

with open("input_files/R1.dat", "w") as f:
    for i, rline in enumerate(rlines):
        f.write(str(rline))
        f.write("\n")

# h
hlines = height.astype(str).tolist()
with open("input_files/H1.dat", "w") as f:
    for hline in hlines:
        f.write(hline)
        f.write("\n")
