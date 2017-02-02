import numpy as np
import matplotlib.pyplot as plt

string = np.genfromtxt("string.txt", dtype=str, unpack=True)

mean_sxsx0 = np.zeros(200000)
mean_sysx0 = np.zeros(200000)

for i in range(0, len(string)):
    t, sxsx, sysx0 = np.genfromtxt(string[i], unpack=True)
    mean_sxsx0 += sxsx
    mean_sysx0 += sysx0

mean_sxsx0 /= len(string)
mean_sysx0 /= len(string)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

t = np.linspace(0, 10000, 200000)
senv = np.zeros(200000)
senv = np.sqrt(mean_sxsx0**2 + mean_sysx0**2)/3

plt.figure(1)
f, axarr = plt.subplots(3, sharey=True)
axarr[0].plot(t, senv, 'b-', linewidth=2)
axarr[0].set_xlim(0, 100)
axarr[1].plot(t, senv, 'b-', linewidth=2)
axarr[1].set_xlim(200, 300)
axarr[2].plot(t, senv, 'b-', linewidth=2)
axarr[2].set_xlim(900, 1000)

#plt.plot(t, senv, 'b-', markersize=2, label="r$S_\text{env}(t)$")
f.text(0.04, 0.5, r"$S_\mathrm{env}(t)$", va='center', rotation='vertical')
plt.xlabel(r"$t[J_q^{-1}]$")
axarr[0].grid()
axarr[1].grid()
axarr[2].grid()
axarr[0].legend(loc='best')
axarr[1].legend(loc='best')
axarr[2].legend(loc='best')
plt.savefig("senv.pdf")
