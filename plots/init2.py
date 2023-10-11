# file to create experience 07 initial conditions with 5000 particles

import numpy as np

domain = 100
size = 0.5

# x
with open("output/x.06", "r") as f:
    xlines1 = f.read().split()
    xlines1 = np.asarray(xlines1).astype(float).reshape(1000, 5000)
    xlines2 = np.ones_like(ylines2)
    xlines = np.concatenate((xlines1[-1], xlines2))
    xlines3 = np.ones_like(xlines2)
    xlines = np.concatenate((xlines, xlines3))

with open("input_files/X1.dat", "w") as f:
    for i, xline in enumerate(xlines):
        f.write(str(xline))
        f.write("\n")

# y
with open("output/y.06", "r") as f:
    ylines1 = f.read().split()
    ylines1 = np.asarray(ylines1).astype(float).reshape(1000, 2500)
    ylines2 = np.arange(start=0.5, stop=domain, step=2*size)
    ylines = np.concatenate((ylines1[-1], ylines2))

with open("input_files/Y1.dat", "w") as f:
    for i, yline in enumerate(ylines):
        f.write(str(yline))
        f.write("\n")

# r
with open("output/r.06", "r") as f:
    rlines1 = f.read().split()
    rlines1 = np.asarray(rlines1).astype(float).reshape(1000, 2500)
    rlines2 = 
    rlines = np.concatenate((rlines1[-1], rlines2))

with open("input_files/R1.dat", "w") as f:
    for i, rline in enumerate(rlines):
        f.write(str(rline))
        f.write("\n")

# h
with open("output/h.06", "r") as f:
    hlines1 = f.read().split()
    hlines1 = np.asarray(hlines1).astype(float).reshape(1000, 2500)
    hlines2 = 
    hlines = np.concatenate((hlines1[-1], hlines2))

with open("input_files/H1.dat", "w") as f:
    for i, hline in enumerate(hlines):
        f.write(str(hline))
        f.write("\n")
