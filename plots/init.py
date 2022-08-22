import numpy as np

lines = np.arange(0.5, 100.5)
lines = lines.astype(str)
lines = lines.tolist()
with open("src/Y.dat", "w") as f:
    for line in lines:
        i = 0
        while i < 50:
            f.write(line)
            f.write("\n")
            i += 1

