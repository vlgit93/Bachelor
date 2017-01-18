# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Use global stylize function and global fontsize
from Plots_global import stylize, fontsize, lin_func, pot_func


def plot_S0zS0z(fontsize):
	stylize(fontsize)
	plt.grid(color='grey')
	
	#Autokorrelationsfunktionen:
	#t, SxSx = np.loadtxt('../data/lanczos_gamma0.01_lvl=8_INFTYspins_test.txt', usecols=(0,1), unpack=True)
	#plt.plot(t, SxSx, '-', color='red', linewidth=2, label=r'$\langle S_0^x(t) S_0^x(0) \rangle$')
	#t, SySy = np.loadtxt('../data/lanczos_gamma0.01_lvl=8_INFTYspins_test.txt', usecols=(0,2), unpack=True)
	#plt.plot(t, SySy, '-', color='green', linewidth=2, label=r'$\langle S_0^y(t) S_0^y(0) \rangle$')
	#t, SzSz = np.loadtxt('../data/lanczos_gamma0.01_lvl=8_INFTYspins_test.txt', usecols=(0,3), unpack=True)
	#plt.plot(t, SzSz, '-', color='blue', linewidth=2, label=r'$\langle S_0^z(t) S_0^z(0) \rangle$')
	t, SxSx = np.loadtxt('../data/lanczos_gamma0.01_lvl=8_INFTYspins_pulse_test.txt', usecols=(0,1), unpack=True)
	plt.plot(t, SxSx, '-', color='red', linewidth=2, label=r'$\langle S_0^x(t) S_0^x(0) \rangle$')
	t, SySy = np.loadtxt('../data/lanczos_gamma0.01_lvl=8_INFTYspins_pulse_test.txt', usecols=(0,2), unpack=True)
	plt.plot(t, SySy, '-', color='green', linewidth=2, label=r'$\langle S_0^y(t) S_0^y(0) \rangle$')
	t, SzSz = np.loadtxt('../data/lanczos_gamma0.01_lvl=8_INFTYspins_pulse_test.txt', usecols=(0,3), unpack=True)
	plt.plot(t, SzSz, '-', color='blue', linewidth=2, label=r'$\langle S_0^z(t) S_0^z(0) \rangle$')
	
	plt.xlabel(r'$t\,\left[J_\mathrm{Q}^{-1}\right]$')
	plt.ylabel(r'$\langle S_0^z(t) S_0^z(0) \rangle$')


	plt.legend(loc='best', prop={'size':fontsize}, fancybox=True, numpoints=1)

 
	plt.savefig('./Abbildungen/lanczos_S0zS0z.pdf', bbox_inches='tight')
	plt.show()
	
	plt.close()


def plot_Senv(fontsize):
	stylize(fontsize)
	plt.grid(color='grey')
		
	t, Sx = np.loadtxt('../data/lanczos_gamma0.1_lvl=8_INFTYspins_pulse_test.txt', usecols=(0,4), unpack=True)
	t, Sy = np.loadtxt('../data/lanczos_gamma0.1_lvl=8_INFTYspins_pulse_test.txt', usecols=(0,5), unpack=True)
	
	S_env_prime = np.sqrt(Sx*Sx + Sy*Sy)
	plt.plot(t, S_env_prime, '-', color='red', linewidth=2, label=r'$S_\mathrm{env}^\prime$')
	
	t, SxSx = np.loadtxt('../data/lanczos_gamma0.1_lvl=8_INFTYspins_pulse_test.txt', usecols=(0,1), unpack=True)
	t, SxSy = np.loadtxt('../data/lanczos_gamma0.1_lvl=8_INFTYspins_pulse_test.txt', usecols=(0,7), unpack=True)
	
	S_env = np.sqrt(SxSx*SxSx + SxSy*SxSy)
	
	plt.plot(t, S_env/3., '-', color='blue', linewidth=2, label=r'$S_\mathrm{env}$')
	
	
	plt.xlabel(r'$t\,\left[J_\mathrm{Q}^{-1}\right]$')
	plt.ylabel(r'$S_\mathrm{env}^\prime (t)$')

	

	plt.legend(loc='best', prop={'size':fontsize}, fancybox=True, numpoints=1)

 
	plt.savefig('./Abbildungen/lanczos_Senv.pdf', bbox_inches='tight')
	plt.show()
	
	plt.close()
	





#plot_S0zS0z(fontsize())
plot_Senv(fontsize())


