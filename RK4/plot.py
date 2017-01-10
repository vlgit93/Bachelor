import numpy as np
import matplotlib.pyplot as plt

x, y = np.genfromtxt("./data/cs_10000.dat", unpack=True)

plt.figure(1)
plt.plot(x, y, "b-", linewidth=2, label="1e4 Durchläufe")
plt.xlim(0, 100)
#plt.ylim(-0.1, 0.1)
plt.xlabel("t$\left[J_q^{-1}\\right]$")
plt.ylabel("$\langle S_0^z(t)S_0^z(0)\\rangle$")
plt.grid()
plt.legend(loc='best')
plt.savefig("plot.pdf")
