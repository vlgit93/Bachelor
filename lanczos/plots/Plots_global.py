# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

def lin_func(x, a, b):
	return (a*x+b)
	
def pot_func(x, a, b):
	return (a*x**b)

def stylize(fontsize):
	plt.rc('text', usetex=True)
	plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],'size':fontsize})
	plt.rcParams['ytick.major.pad']='8'
	plt.rcParams['xtick.major.pad']='8'
	plt.minorticks_on()

def fontsize():
	return 18
