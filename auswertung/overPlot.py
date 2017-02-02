import numpy as np
import matplotlib.pyplot as plt

p0, p50, p200, p1000 = np.genfromtxt("./3_noB.overhauserDist.dat", unpack=True)
p_2_0, p_2_50, p_2_200, p_2_1000 = np.genfromtxt("./6.overhauserDist.dat", unpack=True)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(1)
plt.hist(p0, bins=401, label="$n_\text{Puls}=0$")
plt.xlabel(r"$B^z$[$J_q$]")
plt.ylabel("Anzahl")
plt.savefig("./plots/10_10_1.pdf")

plt.figure(2)
plt.hist(p50, bins=401, label="$n_\text{Puls}=50$")
plt.xlabel(r"$B^z$[$J_q$]")
plt.ylabel("Anzahl")
plt.savefig("./plots/10_10_50.pdf")

plt.figure(3)
plt.hist(p200, bins=401, label="$n_\text{Puls}=200$")
plt.xlabel(r"$B^z$[$J_q$]")
plt.ylabel("Anzahl")
plt.savefig("./plots/10_10_200.pdf")

plt.figure(4)
plt.hist(p1000, bins=401, label="$n_\text{Puls}=1000$")
plt.xlabel(r"$B^z$[$J_q$]")
plt.ylabel("Anzahl")
plt.savefig("./plots/10_10_1000.pdf")
