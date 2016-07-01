import numpy as np
import matplotlib.pyplot as plt

k = np.genfromtxt("k_h.txt")
P = np.genfromtxt("p_k.txt")

print k.shape, P.shape

plt.loglog(k,P[0])
plt.loglog(k,P[1],ls="--")
plt.show()
