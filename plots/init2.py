# file that creates initial conditions from experience 09 to experience 01 the reason why we need this is because the 5000 particles in 09 were costing too much computational time and we need to cuty the number of particles. In this I will cuty everything above 50 km (essentially cutyting the number of particles roughly in half). this should speed up the code a lot so that we can day daily simulations

import numpy as np
import subprocess

# parameters
# ---------------------------
expno = "08"
n = 5000
cuty = 30000
cutx = 10000
position = 12000
# ---------------------------

# compute the number of lines in the files to get the last one
out = subprocess.run(["wc", "-l", "output/y." + expno], capture_output=True, text=True)
num = int(out.stdout.split()[0])

# cut in y
# ---------------------------
with open("output/y." + expno, "r") as f:
    ylines = f.read().split()
    ylines = np.asarray(ylines).astype(float).reshape(num, n)
    ylast = ylines[-1]
    idx_keepy = np.argwhere(ylast < cuty)
    ylast = ylast[idx_keepy]

with open("files/y.dat", "w") as f:
    for i, yline in enumerate(ylast):
        f.write(str(yline[0]))
        f.write("\n")

vars = ["x", "h", "r", "u", "v", "theta", "omega"]
# all others
for var in vars:
    with open("output/" + var + "." + expno, "r") as f:
        lines = f.read().split()
        lines = np.asarray(lines).astype(float).reshape(num, n)
        last = lines[-1]
        last = last[idx_keepy]

    with open("files/" + var + ".dat", "w") as f:
        for i, line in enumerate(last):
            f.write(str(line[0]))
            f.write("\n")


# same but now we cut the above in x
out = subprocess.run(["wc", "-l", "files/x.dat"], capture_output=True, text=True)
num = int(out.stdout.split()[0])

# cut in x
# ---------------------------
with open("files/x.dat", "r") as f:
    xlines = f.read().split()
    xlines = np.asarray(xlines).astype(float)
    idx_keepx = np.argwhere(xlines < cutx)
    xlast = xlines[idx_keepx]

with open("files/x1.dat", "w") as f:
    for i, xline in enumerate(xlast):
        f.write(str(xline[0]))
        f.write("\n")
    # this is for 2nd plate
    for i, xline in enumerate(xlast):
        f.write(str(xline[0] + position))
        f.write("\n")

vars = ["y", "h", "r", "u", "v", "theta", "omega"]
# all others
for var in vars:
    with open("files/" + var + ".dat", "r") as f:
        lines = f.read().split()
        lines = np.asarray(lines).astype(float)
        last = lines[idx_keepx]

    with open("files/" + var + "1.dat", "w") as f:
        for i, line in enumerate(last):
            f.write(str(line[0]))
            f.write("\n")
        # this is for 2nd plate
        for i, line in enumerate(last):
            f.write(str(line[0]))
            f.write("\n")
