import numpy as np
import matplotlib.pyplot as plt

r = np.loadtxt("out.dat")

r0 = r[0]
r -= r0
plt.plot(r)
plt.show()


