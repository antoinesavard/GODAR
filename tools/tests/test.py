import numpy as np
import matplotlib.pyplot as plt

x = np.random.lognormal(6, 0.4, 525)

count, bins, ignored = plt.hist(x, 20)

plt.show()
