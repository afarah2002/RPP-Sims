#----------PYTHON imports----------#
import collections
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import random
# import thread
import time
import scipy.stats as ss
from scipy import optimize, stats
import sys

#----------code starts here!----------#
pi = np.pi
SEE_param_dict = {
				   "Be" : (7.14, 1.803, 15.36, 1.0, 0.55, 0.20),  
				   "Mg" : (4.16, 1.764, 12.97, 1.3, 0.80, 0.24), 
				   "Al" : (3.80, 1.783, 7.970, 1.7, 2.00, 0.40), 
				   "Si" : (4.05, 1.625, 13.6,  2.7, 0.89, 0.45),
				   "Ti" : (3.13, 1.891, 2.84,  0.5, 1.21, 0.25), 
				   "Fe" : (3.23, 1.732, 2.761, 0.6, 1.15, 0.35), 
				   "Ni" : (3.07, 1.724, 2.414, 1.0, 1.19, 0.50), 
				   "Cu" : (3.02, 1.719, 2.544, 0.6, 1.53, 0.40), 
				   "Ge" : (2.80, 1.708, 4.480, 1.0, 1.00, 0.40), 
				   "Pd" : (2.46, 1.660, 2.251, 1.0, 1.41, 0.55), 
				   "Ag" : (2.64, 1.659, 2.560, 1.0, 1.43, 0.60), 
				   "W"  : (2.29, 1.575, 1.975, 0.2, 1.06, 0.25), 
				   "Pt" : (2.29, 1.564, 1.874, 0.5, 1.69, 0.55), 
				   "Au" : (2.35, 1.559, 2.11,  0.5, 1.28, 0.50), 
				   "Pb" : (2.38, 1.553, 3.70,  0.8, 1.06, 0.50), 
				   "Bi" : (2.39, 1.549, 4.84,  2.0, 0.98, 0.70), 
				 }
class PrimaryCurrentPlotter(object):

	def __init__(self, transProb, energyLB, energyUB, solRad, solTurns):
		self.transProb = transProb
		self.energyUnit = "(keV)"
		self.currentUnit = "(A)"
		self.cluster_constant = 4.644e-12
		self.E_PE_list = np.arange(energyLB, energyUB, (energyUB-energyLB)/100000)
		self.solRad = solRad
		self.solTurns = solTurns
		self.clusterArea = pi*(0.188/2)**2
		# print self.E_PE_list

	def computePopulation(self, element, params):
		vac_perm_constant = 4*pi*10**-7
		m, n, B, escape_depth, delta_max, E_PE_max = params

		# find the energy that matches the given transmission probability

		self.T_r_list = []
		# print element
		E_opt_list = []
		for e in self.E_PE_list:
			T_r = np.round(np.e**(-(escape_depth/(B*e**n))**m), 5)
			# print element, T_r, e
			self.T_r_list.append(T_r)
			if np.isclose(T_r, self.transProb,1e-6,1e-4) == True:
				E_opt_list.append((T_r, e))
				# print element, T_r, e

		self.transmitted_percent, self.E_PE = E_opt_list[0] # takes the option with the least energy
		# print element, E_opt_list[0]

		SEE_yield = delta_max*((1.28*(self.E_PE/E_PE_max)**-0.67)*(1-np.e**(-1.614*(self.E_PE/E_PE_max)**1.67)))
		self.I_PE = self.solRad*np.sqrt(125*e*self.cluster_constant)/ \
			   (4*vac_perm_constant*self.solTurns * SEE_yield)

		# number of indident e+ through cluster area to get the incident current
		self.population = self.I_PE/(1.6e-19*self.clusterArea*3e8*np.sqrt(2*self.E_PE/511))
		self.pop_trans = self.population*self.transmitted_percent
		print element, "\n", \
			  "beam population: ", np.format_float_scientific(self.population, precision=3), " e+", "\n", \
			  "beam energy: ", np.round(self.E_PE,3), " keV", "\n", \

		return (self.E_PE, self.transmitted_percent), self.population

	


	def plotter(self):
		# fig, ax = plt.subplots()
		plt.figure(0)
		x = self.E_PE_list
		# y = self.I_PE_list
		y = self.T_r_list
		plt.xlim(0,max(x)+max(x)/50)
		plt.ylim(0,1.1)
		x_label = "Incident e+ energy " + self.energyUnit
		# y_label = "Incident e+ current " + self.currentUnit
		y_label = "e+ transmission probability"
		plt.plot(x,y,"--")
		plt.xlabel(x_label)
		plt.ylabel(y_label)
		plt.title(y_label + " vs " + x_label)
		plt.grid(1)
		pass


transmission_Probability = 0.9999
PCP = PrimaryCurrentPlotter(transmission_Probability, 0.001, 8, 1, 500)

if __name__ == '__main__':
	PCP
	energies = []
	populations = []
	# fig1, ax1 = plt.subplots()
	for element, params in SEE_param_dict.items():
		# print element
		opt_coor, population = PCP.computePopulation(element, params)
		energies.append(opt_coor[0])
		populations.append(population)

		PCP.plotter()
		plt.scatter(opt_coor[0], opt_coor[1])
		plt.text(opt_coor[0], opt_coor[1], element, fontsize=12)
	
		plt.figure(1)
		plt.scatter(opt_coor[0], population)
		plt.xlabel("energies (keV)")
		plt.ylabel("population")
		plt.xlim(0,8)
		plt.ylim(0,1.5e12)
		plt.grid(1)
		plt.text(opt_coor[0], population, element, fontsize=12)
	plt.show()
	# for i, element in enumerate(SEE_param_dict.keys()):
	# 	coor = mins_list[i]
	# 	ax.scatter(coor[0],coor[1])
	# 	ax.text(coor[0], coor[1], element,fontsize=12)

	# plt.title("Minimum e+ currents and energies")
	# x_label = "Incident e+ energy (keV)" 
	# y_label = "Incident e+ current (A)" 
