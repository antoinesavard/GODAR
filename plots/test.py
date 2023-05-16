import numpy as np
import matplotlib.pyplot as plt

x = np.random.lognormal(0, 0.25, 2500) - 0.15
y = np.random.lognormal(6, 0.25, 2500)

count, bins, ignored = plt.hist(y, 20, density=1)


plt.show()
