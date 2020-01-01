#----------GEANT4 imports----------#
from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam

#----------PYTHON imports----------#
import collections
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import random
import thread
import time
import scipy.stats as ss
from scipy import optimize, stats
import sys

#-----------FILE imports------------#
from beam2 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction,MySteppingAction, \
				  Plotter, WipeData
from electron_emission import SecondaryElectronEmissionProcess
from geom_constructor import GeomConstructor 
from visualizer import Visualizer


PLT = Plotter()
WIPE = WipeData()
GC = GeomConstructor()

spherical_coor_LIST = []
# spherical_coor_LIST = [[0,0]]
pi = np.pi
step = 6
for phi in np.arange(0, pi, pi/step): # smaller steps means more clusters, range goes to pi since clusters are double sided
	for theta in np.arange(0, pi, pi/step):
		sph_coor = [theta, phi] # phi, theta
		spherical_coor_LIST.append(sph_coor)

cluster_coor_LIST = []
global cluster_width
cluster_width = 160
spacing = 40./7
interval = cluster_width + spacing
interval *= 3
global edge
edge = 500
for x in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval, interval):
	for y in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval, interval):
		for z in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval, interval):
			cluster_coor = [x,y,z]
			# cluster_coor_LIST.append(cluster_coor)
			# print cluster_coor
			# time.sleep(1)
			opposite = np.multiply(cluster_coor, -1)
			for i in cluster_coor:
				if np.abs(i) == edge-cluster_width/2 - spacing and cluster_coor not in cluster_coor_LIST:
					cluster_coor_LIST.append(cluster_coor)
					break


# print cluster_coor_LIST
'''
# USE TO DOUBLE CHECK CLUSTER COOR LIST #

fig = plt.figure()
# Axes3D.scatter(self.px, self.py, self.pz)

# first subplot: a 3D scatter plot of positions
ax = fig.add_subplot(111, projection='3d')
axmin = -(500-cluster_width/2 - spacing)
axmax = 500 - (cluster_width/2 + spacing)	
axes = plt.gca()
axes.set_xlim([axmin,axmax])
axes.set_ylim([axmin,axmax])
axes.set_zlim([axmin,axmax])

ax.set_xlabel('mm')
ax.set_ylabel('mm')
ax.set_zlabel('mm')

for cluster_coor in cluster_coor_LIST:
	print cluster_coor
	ax.scatter(cluster_coor[0], cluster_coor[1], cluster_coor[2])

# plt.draw() 

plt.show()
'''

# energy_LIST = [0.002]
# energy_LIST = list(np.arange(0.002, 0.003, 0.0002))
energy_LIST = list(np.arange(1001, 2000, 200))
global energyUnit
energyUnit = eV

cluster_time_median_LIST = []

global opt_be_right_LIST
global opt_be_left_LIST
opt_be_right_LIST = []
opt_be_left_LIST = []

class Constructor(object):
	def __init__(self):
		# from EZsim.EZgeom import G4EzVolume
		g4py.NISTmaterials.Construct()
		# set DetectorConstruction to the RunManager
		# G4VUserDetectorConstruction.Construct()
		g4py.ezgeom.Construct()
		exN03PL = g4py.EMSTDpl.PhysicsListEMstd()
		gRunManager.SetUserInitialization(exN03PL)
		# reset world material
		air= gNistManager.FindOrBuildMaterial("G4_AIR")
		# g4py.ezgeom.SetWorldMaterial(air)

		# myDC= MyDetectorConstruction()
		# gRunManager.SetUserInitialization(myDC)



	def construct(self):
		# dummy box
		# GC = GeomConstructor()

		# calorimeter placed inside the box
		# cal= G4EzVolume("Calorimeter") #initialize volume
		NaI= gNistManager.FindOrBuildMaterial("G4_SODIUM_IODIDE")
		# GRAPHITE = gNistManager.FindOrBuildMaterial("G4_GRAPHITE")
		# GC.ConstructBox("Detector Box", GRAPHITE, [0., 0., 0.], cm, [95., 95., 95.])

		# GC.ConstructOrb("Orb", NaI, [10., 10., 10.], cm, 10.)

		# GC.ConstructTube("Tube", NaI, [-1., -1., -1.], cm, 10., 30., 30., 0, 300.)
		scale_factor = 5
		# GC.ConstructSphere("Sphere", material1 , [0., 0., 0.], cm, 0, 1., 0., 360., 0., 180)
		# GC.ConstructSphere("Sphere", GRAPHITE , [0., 0., 0.], m, (.10-30e-9)*scale_factor, .10*scale_factor, 0., 360., 0., 90)


		# GC.ConstructCone("Cone", NaI, [-20., -20., -20.], cm, 0., 20., 0., 0., 25., 0., 180) # dphi = 359.9999 is basically 360, but we can still see it
		# gRunManager.Initialize()

###############################################################################################################################################
# Functions for curve fitting by scipy.optimize.curve_fit()
class BaseFunctions(object):
	
	def rational(self, x, p, q):
	    """
	    The general rational function description.
	    p is a list with the polynomial coefficients in the numerator
	    q is a list with the polynomial coefficients (except the first one)
	    in the denominator
	    The zeroth order coefficient of the denominator polynomial is fixed at 1.
	    Numpy stores coefficients in [x**2 + x + 1] order, so the fixed
	    zeroth order denominator coefficent must comes last. (Edited.)
	    """
	    return np.polyval(p, x) / np.polyval(q + [1.0], x)

BF = BaseFunctions()

def rational3_3(x, p0, p1, p2, q1, q2):
    return BF.rational(x, [p0, p1, p2], [q1, q2])

class CurveFitter(object):

	def fit(self, test_function, x_data, y_data):
		popt, pcov = optimize.curve_fit(test_function, x_data, y_data, p0=None)
		print popt, "\n"
		return popt
CF = CurveFitter()

uniqueClusters = [[None, None, None]]
class FieldDesign(object):

	def cartesianfieldParam(self, energy, b, x, y, z, cluster_width):
		
		radius = np.sqrt(x**2 + y**2 + z**2)

		mag0 = x*b/radius
		mag1 = y*b/radius
		mag2 = z*b/radius

		# mag0 = np.sqrt(3*b**2)
		# mag1 = np.sqrt(3*b**2)
		# mag2 = np.sqrt(3*b**2)

		magVec = [mag0, mag1, mag2] 
		# print magVec
		
		# get index of "max" value in magVec, most significant axis
		maxIndex = list(np.abs(magVec)).index(max(np.abs(magVec))) 

		#get scaling factor
		scale = edge / np.abs(magVec[maxIndex])
		magVecScaled = np.multiply(magVec, scale)

		receiverDimScaled = [cluster_width, cluster_width, cluster_width]
		receiverDimScaled[maxIndex] = 1

		# print "Receiver DIMENSIONS  ", receiverDimScaled

		material1 = G4Material.GetMaterial("G4_W")
		GC.ConstructBox("Receiver", material1, [x,y,z], mm, receiverDimScaled)
		GC.ConstructBox("Receiver", material1, np.multiply([x,y,z], -1), mm, receiverDimScaled) # receiver for opposite cluster

		# if the magVecScaled is not in the uniqueClusters list, append to it
		flag = 0
		for elem in uniqueClusters:
			if collections.Counter(elem) == collections.Counter(magVecScaled):
				# print True
				flag = 1

		if flag == 0:
			uniqueClusters.append(magVecScaled)
		else: 
			pass


		return magVec, magVecScaled

	def spherefieldParam(self, energy, b, phi, theta):

		# vectorList = [list(np.multiply([energy, energy, energy], be))]
		# radius = np.sqrt(np.square(vectorList[0][0]) + np.square(vectorList[0][1]) + np.square(vectorList[0][2]))

		mag0 = b*(np.sin(phi)*np.cos(theta))
		mag1 = b*(np.sin(phi)*np.sin(theta))
		mag2 = b*np.cos(phi)

		magVec = [mag0, mag1, mag2] 
		
		# get index of "max" value in magVec, most significant axis
		maxIndex = list(np.abs(magVec)).index(max(np.abs(magVec))) 

		#get scaling factor
		scale = 500 / np.abs(magVec[maxIndex])
		magVecScaled = np.multiply(magVec, scale)

		receiverDimScaled = [cluster_width, cluster_width, cluster_width]
		receiverDimScaled[maxIndex] = 1

		# print "Receiver DIMENSIONS  ", receiverDimScaled

		material1 = G4Material.GetMaterial("G4_W")
		# GC.ConstructBox("Receiver", material1, magVecScaled, mm, receiverDimScaled)
		# GC.ConstructBox("Receiver", material1, np.multiply(magVecScaled, -1), mm, receiverDimScaled) # receiver for opposite cluster

		# if the magVecScaled is not in the uniqueClusters list, append to it
		flag = 0
		for elem in uniqueClusters:
			if collections.Counter(elem) == collections.Counter(magVecScaled):
				flag = 1

		if flag == 0:
			uniqueClusters.append(magVecScaled)
		else: 
			pass


		return magVec, magVecScaled

# -----------CLASS ASSIGNMENTS----------- #
FD = FieldDesign()
Constructor = Constructor()
VIS = Visualizer()
SEEP = SecondaryElectronEmissionProcess()
# --------------------------------------- #

initialMomenta = []
finalMomenta = []

viz_theta = 35
viz_phi = 35
zoom = 1.5


global ALL_clusters_positions 
global ALL_clusters_momenta 
ALL_clusters_positions = []
ALL_clusters_momenta = []

class ClusterClass(object):

	def run(self, energy, location_range):
		# print(energy_LIST)
		# print(cluster_coor_LIST)
		# time.sleep(1)
		data_right = open("data_right.txt", "a")
		data_left = open("data_left.txt", "a")

		# print spherical_coor_LIST

		# time.sleep(1)

		e = energy

		ALL_clusters_positions[:] = []
		ALL_clusters_momenta[:] = []

		if energyUnit == MeV:
			constant = 4.644e-9
		if energyUnit == keV:
			constant = 4.644e-12
		if energyUnit == eV:
			constant = 4.644e-15

		b = np.sqrt(e*3*constant) 

		for location in location_range:
			# Constructor
			# Constructor.construct()
			if len(location) == 3: # this is cartesian
				x = location[0]
				y = location[1]
				z = location[2]
				magVec, magVecScaled = FD.cartesianfieldParam(e, b, x, y ,z, cluster_width)
			if len(location) == 2: # this is spherical
				phi = location[0]
				theta = location[1]
				magVec, magVecScaled = FD.spherefieldParam(e, b, phi, theta)

			# set user actions ...
			# PGA_1 = MyPrimaryGeneratorAction(e, energyUnit,cart_coor)
			PGA_1 = MyPrimaryGeneratorAction(e, energyUnit, magVecScaled)
			gRunManager.SetUserAction(PGA_1)
			myEA = MyEventAction()
			gRunManager.SetUserAction(myEA)
			mySA = MySteppingAction()
			gRunManager.SetUserAction(mySA)
			# VIS.visualizer(viz_theta, viz_phi)
			fieldMgr = gTransportationManager.GetFieldManager()
			myField = G4UniformMagField(G4ThreeVector(magVec[0], magVec[1], magVec[2]))
			# myField = MyField(1)
			fieldMgr.SetDetectorField(myField)
			fieldMgr.CreateChordFinder(myField)
			# print "|B-field| = ", vectorList[0][0]
			# CI = ClusterIsolation(magVecScaled)
			# CI
			# CI.getClusterWidth()
			myRA = MyRunAction()
			gRunManager.SetUserAction(myRA)
			gRunManager.Initialize()
			gRunManager.BeamOn(1)

			cluster_positions, cluster_momenta = PLT.clusterDataReturner()
			print len(cluster_positions), "\n", len(cluster_momenta)

			# saving ALL clusters for this energy to be SEEPed
			ALL_clusters_positions.append(cluster_positions) 
			print "appended", ALL_clusters_positions
			# time.sleep(12)
			ALL_clusters_momenta.append(cluster_momenta) 
			SEEP.runSEE(e, cluster_positions, cluster_momenta)
			# # std_devs_LIST, means_LIST, n_LIST, n_sd_LIST = PLT.dataReturner()
			std_devs_LIST, n_LIST, n_sd_LIST, cluster_time_LIST, cluster_size_LIST = PLT.dataReturner() # for 3D positions

				# WIPE.wipeCluster()
				# PLT.wipeData() #clean lists before starting another run
			'''
			for num, n_sd in enumerate(n_sd_LIST):
				if num == 0: # right cluster
					print "NSD RIGHT"
					# time.sleep(1)
					max_n_sd = max(n_sd) # gets the peak n/sd that shows the greatest clustering efficiency
					opt_be = n_sd.index(max_n_sd) * be_step # find the be_ratio that produces that max n/sd 

					opt_be_right_LIST.append(opt_be)
					data_right.write(str(opt_be)+"\n")
				else: # left cluster
					print "NSD LEFT"
					# time.sleep(1)
					max_n_sd = max(n_sd) # gets the peak n/sd that shows the greatest clustering efficiency
					opt_be = n_sd.index(max_n_sd) * be_step # find the be_ratio that produces that max n/sd 

					opt_be_left_LIST.append(opt_be)
					data_left.write(str(opt_be)+"\n")
			'''
			data  = {
					"SD" : std_devs_LIST, \
					# "means" : means_LIST, \
					"n/sd" : n_sd_LIST, \
					"cluster_size" : n_LIST \
					}
			median = np.median(cluster_time_LIST)
			cluster_time_median_LIST.append(median)
							
			# graph the cluster time distribution
			# n_bins = 50
			# # mu, sigma = np.mean(cluster_time_LIST), np.std(cluster_time_LIST)

			# n, bins, patches = plt.hist(cluster_time_LIST, n_bins, normed=0, facecolor='red', alpha=0.5)
			# plt.xlabel("Clustering time (ns)")
			# plt.ylabel("Frequency")
			# plt.title("Clustering Time distribution (ns)")
			# cluster_time_LIST = [round(t, 2) for t in cluster_time_LIST]
			# print cluster_time_LIST
			# print "time median (ns): ", median, "\n"
			   #    "time mean (ns): ", np.mean(cluster_time_LIST), "\n", \
				  # "time mode (ns): ", stats.mode(cluster_time_LIST).mode[0], "\n" 
			'''SUMMARY - MODE IS THE BEST ESTIMATE OF THE CENTER BECAUSE 
				  	 MEAN IS HIGHER THAN MOST OF THE DATA AND MODE IS LOWER THAN
				  	 MOST OF THE DATA'''
			# plt.show()
			
			# WIPE.wipeTime()
		time.sleep(1)
		for i in ALL_clusters_positions:
			print i[-1] # <------------------- WHAT?!?!?!? HOW, THESE SHOULD BE DIFFERENT!!!! APPENDING RIGHT????
		time.sleep(1.5)
		return ALL_clusters_positions, ALL_clusters_momenta

	# 		'''
	# 		fig, (sd, n_sd, n) = plt.subplots(3, sharex=True, sharey=False)
	# 		# plt.tight_layout()
	# 		plt.xlabel("Ratio of B-field to Particle Beam Energy (T/keV) ", fontsize=18)

	# 		fontdict = {'fontsize': 18,
	# 					'fontweight': 5,
	# 					}

	# 		for dep_var_name, dep_var_LIST in data.items():
	# 			for dep_var in dep_var_LIST:

	# 				if dep_var_LIST.index(dep_var) == 0:
	# 					label = 'right'
	# 				if dep_var_LIST.index(dep_var) == 1:
	# 					label = 'left'

	# 				if dep_var_name == 'SD':
	# 					sd.plot(be_ratio, dep_var, label=label)	
	# 					sd.set(ylabel=dep_var_name)
	# 					title = dep_var_name + " vs be_ratio (T/eV)"
	# 					sd.set_title(title, fontdict=fontdict)
	# 					# popt = CF.fit(function, be_ratio, dep_var)
	# 					# sd.plot(be_ratio, function(be_ratio, *popt), label=label)
	# 				if dep_var_name == "n/sd":
	# 					n_sd.plot(be_ratio, dep_var, label=label)	
	# 					n_sd.set(ylabel=dep_var_name)
	# 					title = dep_var_name + " vs be_ratio (T/eV)"
	# 					n_sd.set_title(title, fontdict=fontdict)
	# 				if dep_var_name == "cluster_size":				
	# 					n.plot(be_ratio, dep_var, label=label)
	# 					n.set(ylabel=dep_var_name)
	# 					title = dep_var_name + " vs be_ratio (T/eV)"
	# 					n.set_title(title, fontdict=fontdict)

	# 		plt.xticks(tickMarks)
	# 		plt.legend()
	# 	# plt.show()
	# 	'''
	# 	#### plotting E vs optimal B/E ####
	# 	# plt.figure()
	# 	# plt.xlabel("Energy (MeV)", fontsize=18)
	# 	# plt.ylabel("Optimal B/E ratio (mT/MeV)", fontsize=18)
	# 	# plt.ylim(0,0.0002*1000)
	# 	# plt.title("Optimal B/E ratio (mT/MeV) vs. e+ Energy (MeV)", fontsize=24)
	# 	# plt.plot(dummy_x, np.multiply(gathered_data_right, 1000), label='right')
	# 	# plt.plot(dummy_x, np.multiply(gathered_data_left, 1000), label='left')
	# 	# # plt.ticklabel_format(axis='both', style='sci', scilimits=(-7,0))
	# 	# plt.legend()
	# 	# plt.show()

	# # fit median data
	# # x = np.arange(0, len(cluster_time_median_LIST))
	# # y = cluster_time_median_LIST
	# # function = rational3_3
	# # popt = CF.fit(function, x, y)

	# '''
	# PLT.grapher()

	# fig, ax = plt.subplots(1, sharey=True, sharex=False, tight_layout=False)
	# n_bins = 10
	# ax.hist(cluster_size_LIST, n_bins)
	# plt.title("Cluster sizes (mm)")
	# # plt.figure()
	# # plt.ylabel("Median cluster time (ns)", fontsize=18)
	# # plt.xlabel("Run number", fontsize=18)
	# # plt.plot(x,y)
	# # plt.plot(np.arange(-60, 60, 1), function(np.arange(-60, 60, 1), *popt))
	# cluster_count = 2*len(cluster_coor_LIST)
	# true_cluster_count = (1000/cluster_width)**3
	# avg_cluster_size = np.mean(cluster_size_LIST)
	# # cluster_count = 2 * (len(uniqueClusters)-1)
	# print "There are ", cluster_count, " clusters"
	# print "Avg cluster size: ", avg_cluster_size, " mm"
	# plt.show()
	# '''

CC = ClusterClass()
if __name__ == '__main__':

	for e in energy_LIST:
		positions, momenta = CC.run(e, spherical_coor_LIST)
		# print "# of clusters = ", len(positions)
		# time.sleep(3)

		# print len(positions)
		# time.sleep(1)
		# counter = 0
		# for pos in positions:
		# 	mom = momenta[positions.index(pos)]
		# 	# print mom
		# 	SEEP.runSEE(e, pos, mom)
		# 	# time.sleep(3)
