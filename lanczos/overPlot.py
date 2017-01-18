import numpy as np
import matplotlib.pyplot as plt

p0, p50, p200, p1000 = np.genfromtxt("./data/overhauserDist.dat", unpack=True)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(1)
plt.hist(p0, bins=401, label="$n_\text{Puls}=0$")
plt.xlabel(r"$B^z$[$J_q$]")
plt.ylabel("Anzahl")
plt.savefig("overDist.pdf")
