# file that creates initial conditions from experience 09 to experience 01 the reason why we need this is because the 5000 particles in 09 were costing too much computational time and we need to cut the number of particles. In this I will cut everything above 50 km (essentially cutting the number of particles roughly in half). this should speed up the code a lot so that we can day daily simulations

import numpy as np
import subprocess

expno = "08"
n = 5000
cut = 50000
saveas = "10"

out = subprocess.run(["wc", "-l", "output/x." + expno], capture_output=True, text=True)
num = int(out.stdout[6:8])

# y
with open("output/y." + expno, "r") as f:
    ylines = f.read().split()
    ylines = np.asarray(ylines).astype(float).reshape(num, n)
    ylast = ylines[-1]
    idx_keep = np.argwhere(ylast < cut)
    ylast = ylast[idx_keep]

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
        last = last[idx_keep]

    with open("files/" + var + ".dat", "w") as f:
        for i, line in enumerate(last):
            f.write(str(line[0]))
            f.write("\n")
