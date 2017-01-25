import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#AUTOCORRELATIONFUNCTION z-COMPONENT
t, sz = np.genfromtxt("./data/cs_iter=1e2_N=1e3_h0z=0.dat", unpack=True)

plt.figure(1)
plt.plot(t, sz, "b-", linewidth=2, label=r"1e3 Durchl\"aufe")
plt.ylabel(r"$S_\mathrm{env}(t)$")
plt.xlabel(r"t$\left[J_q^{-1}\right]$")
plt.grid()
plt.legend(loc='best')
plt.savefig("plot.pdf")


"""
t, sx, sy, sz = np.genfromtxt("./data/cs_iter=1e2_N=1e2_h0z=5.dat", unpack=True)

f, axarr = plt.subplots(3, sharex=True)
axarr[0].plot(t, sx, "b-", linewidth=2, label=r"a")
axarr[1].plot(t, sy, "b-", linewidth=2, label=r"b")
axarr[2].plot(t, sz, "b-", linewidth=2, label=r"c")

axarr[0].set_ylabel(r"$\langle S_0^x(t)S_0^x(0)\rangle$")
axarr[1].set_ylabel(r"$\langle S_0^y(t)S_0^y(0)\rangle$")
axarr[2].set_ylabel(r"$\langle S_0^z(t)S_0^z(0)\rangle$")
plt.xlabel(r"t$\left[J_q^{-1}\right]$")
axarr[0].grid()
axarr[1].grid()
axarr[2].grid()
axarr[0].legend(loc='best')
axarr[1].legend(loc='best')
axarr[2].legend(loc='best')
plt.savefig("plot.pdf")
"""
