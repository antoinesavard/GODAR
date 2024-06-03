import numpy as np
import matplotlib.pyplot as plt

x = np.random.normal(450, 150, 525)
y = np.random.lognormal(6, 0.5, 525)

count, bins, ignored = plt.hist(x, 20)


plt.show()
