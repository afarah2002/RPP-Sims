#----------imports----------#
from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam

import scipy.stats as ss
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages
import random
import time
import thread
import numpy as np
from scipy import optimize, stats
#----file imports--------#
from geom_constructor import GeomConstructor 
# from beam import BeamInitializer
from beam2 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction, MySteppingAction, MyField, Plotter, WipeData


PLT = Plotter()
WIPE = WipeData()

spherical_coor_LIST = []
pi = np.pi
step = 8
for phi in np.arange(0, pi, pi/step): # smaller steps means more clusters, range goes to pi since clusters are double sided
	for theta in np.arange(0, pi, pi/step):
		sph_coor = [theta, phi] # phi, theta
		spherical_coor_LIST.append(sph_coor)


# energy_LIST = lists(np.arange(2., 9., 1.)) # MeV
# energy_LIST = list(np.arange(1., 50., 1.)) # MeV
energy_LIST = [2.5]
dummy_x  = list(np.arange(1., 50., 1.)) # MeV
dummy_y = [0.0001]*49

cluster_time_median_LIST = []

gathered_data_right = [6e-05, \
8e-05, \
0.00013, \
8e-05, \
4e-05, \
5e-05, \
5e-05, \
0.00011, \
5e-05, \
4e-05, \
6e-05, \
8e-05, \
6e-05, \
6e-05, \
7e-05, \
4e-05, \
4e-05, \
5e-05, \
9e-05, \
5e-05, \
7e-05, \
6e-05, \
5e-05, \
3e-05, \
5e-05, \
3e-05, \
4e-05, \
9e-05, \
8e-05, \
5e-05, \
9e-05, \
7e-05, \
7e-05, \
5e-05, \
5e-05, \
4e-05, \
3e-05, \
6e-05, \
6e-05, \
4e-05, \
3e-05, \
5e-05, \
8e-05, \
5e-05, \
8e-05, \
7e-05, \
2e-05, \
5e-05, \
2e-05, \
]
gathered_data_left = [5e-05, \
6e-05, \
5e-05, \
5e-05, \
7e-05, \
0.00011, \
6e-05, \
6e-05, \
0.00016, \
7e-05, \
9e-05, \
6e-05, \
6e-05, \
6e-05, \
3e-05, \
8e-05, \
3e-05, \
4e-05, \
4e-05, \
7e-05, \
5e-05, \
0.0001, \
8e-05, \
0.00011, \
4e-05, \
3e-05, \
4e-05, \
3e-05, \
6e-05, \
7e-05, \
4e-05, \
0.00012, \
3e-05, \
9e-05, \
3e-05, \
3e-05, \
4e-05, \
4e-05, \
5e-05, \
4e-05, \
6e-05, \
6e-05, \
9e-05, \
7e-05, \
5e-05, \
7e-05, \
3e-05, \
3e-05, \
0.00015, \
]

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
		GC = GeomConstructor()
		# GC.ConstructBox("Detector Box", air, [0., 0., 0.], cm, [20., 20., 40.])

		# calorimeter placed inside the box
		# cal= G4EzVolume("Calorimeter") #initialize volume
		nai= gNistManager.FindOrBuildMaterial("G4_SODIUM_IODIDE")
		material1 = G4Material.GetMaterial("G4_W")

		# GC.ConstructOrb("Orb", nai, [10., 10., 10.], cm, 10.)

		# GC.ConstructTube("Tube", au, [-1., -1., -1.], cm, 10., 30., 30., 0, 300.)
		scale_factor = 5
		# GC.ConstructSphere("Sphere", material1 , [0., 0., 0.], cm, 0, 1., 0., 360., 0., 180)
		# GC.ConstructSphere("Sphere", material1 , [0., 0., 0.], cm, 9.5*scale_factor, 10.*scale_factor, 0., 360., 0., 180)

	# GC.ConstructCone("Cone", nai, [-20., -20., -20.], cm, 0., 20., 0., 0., 25., 0., 180) # dphi = 359.9999 is basically 360, but we can still see it
		# gRunManager.Initialize()


class Visualizer(object):

	def visualizer(self, view_angle):
		# time.sleep(.1)

		gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
		gApplyUICommand("/vis/viewer/create OGLSX oglsxviewer")
		gApplyUICommand("/vis/drawVolume")
		# gApplyUICommand("/globalField/verbose level")
		# gApplyUICommand("/MagneticFieldModels/UniformField/SetBVec 3. 1. 1. Tesla")

		gApplyUICommand("/vis/viewer/select oglsxviewer")
		gApplyUICommand("/vis/ogl/set/displayListLimit 100000")
		gApplyUICommand("/vis/scene/add/trajectories")

		gApplyUICommand("/tracking/storeTrajectory 1")
		gApplyUICommand("/vis/scene/endOfEventAction accumulate")
		gApplyUICommand("/vis/scene/endOfRunAction accumulate")

		gApplyUICommand("/vis/viewer/set/viewpointThetaPhi 0" + str(angle))
		gApplyUICommand("/vis/viewer/zoom 1.00001")
		
		gApplyUICommand("/vis/viewer/update")

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


	# def fit(self):
	# 	pass
	# def grapher(self):
	# 	pass
Constructor = Constructor()
Constructor.construct()
VIS = Visualizer()

initialMomenta = []
finalMomenta = []

angle = 35
zoom = 1.5

if __name__ == '__main__':
	# print(energy_LIST)
	print(spherical_coor_LIST)
	# time.sleep(1)
	data_right = open("data_right.txt", "a")
	data_left = open("data_left.txt", "a")
	for sph_coor in spherical_coor_LIST:
		phi = sph_coor[0]
		theta = sph_coor[1]

		print "coors = (", phi, ",", theta, ")"
		for e in energy_LIST:
			WIPE.wipeComps()
			energy = e
			# energy = 2.5
			be_step = 1.e-5

			tickMarks = np.arange(1e-7, 2.5e-4, be_step*5.)
			be_ratio = [1.25e-4]
			# be_ratio = np.arange(1e-7, 2.5e-4, be_step) # ratio between magnetic field (varied) and particle energy (fixed @ 2.5 MeV)
			# be_ratio = [1.e-4, 2.e-4, be_step] # 2.5 MeV after all tests, used to verify best be_ratio, should display 3D position plot
			# be_ratio = [6e-5] # 0.5 MeV after all tests, used to verify best be_ratio, should display 3D position plot
			print energy, "\n"

			print("be len: ", len(be_ratio))
			# time.sleep(1)

			for be in be_ratio:
				# angle += 10
				# angle += 0.075 # +0.075 is a recommended delta theta

				# set user actions ...
				PGA_1 = MyPrimaryGeneratorAction(energy)
				gRunManager.SetUserAction(PGA_1)

				myEA = MyEventAction()
				gRunManager.SetUserAction(myEA)

				mySA = MySteppingAction()
				gRunManager.SetUserAction(mySA)


				vectorList = [list(np.multiply([energy, energy, energy], be))]
				radius = np.sqrt(np.square(vectorList[0][0]) + np.square(vectorList[0][1]) + np.square(vectorList[0][2]))

				fieldMgr = gTransportationManager.GetFieldManager()
				myField = G4UniformMagField(G4ThreeVector(vectorList[0][0]*(np.sin(phi)*np.cos(theta)), \
														  vectorList[0][1]*(np.sin(phi)*np.sin(theta)), \
														  vectorList[0][2]*np.cos(phi)))
				# myField = MyField(1)
				fieldMgr.SetDetectorField(myField)
				fieldMgr.CreateChordFinder(myField)
				# print "|B-field| = ", vectorList[0][0]

				myRA = MyRunAction()
				gRunManager.SetUserAction(myRA)

				gRunManager.Initialize()

				gRunManager.BeamOn(1)

				VIS.visualizer(angle)

				# std_devs_LIST, means_LIST, n_LIST, n_sd_LIST = PLT.dataReturner()
				std_devs_LIST, n_LIST, n_sd_LIST, cluster_time_LIST = PLT.dataReturner() # for 3D positions

				# PLT.wipeData() #clean lists before starting another run

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

			

			data  = {
					"SD" : std_devs_LIST, \
					# "means" : means_LIST, \
					"n/sd" : n_sd_LIST, \
					"cluster_size" : n_LIST \
					}
			median = np.median(cluster_time_LIST)
			cluster_time_median_LIST.append(median)
			'''					
			# graph the cluster time distribution
			n_bins = 50
			# mu, sigma = np.mean(cluster_time_LIST), np.std(cluster_time_LIST)

			n, bins, patches = plt.hist(cluster_time_LIST, n_bins, normed=0, facecolor='red', alpha=0.5)
			plt.xlabel("Clustering time (ns)")
			plt.ylabel("Frequency")
			plt.title("Clustering Time distribution (ns)")


			# cluster_time_LIST = [round(t, 2) for t in cluster_time_LIST]
			# print cluster_time_LIST
			print "time median (ns): ", median, "\n"
			   #    "time mean (ns): ", np.mean(cluster_time_LIST), "\n", \
				  # "time mode (ns): ", stats.mode(cluster_time_LIST).mode[0], "\n" 
			SUMMARY - MODE IS THE BEST ESTIMATE OF THE CENTER BECAUSE 
				  	 MEAN IS HIGHER THAN MOST OF THE DATA AND MODE IS LOWER THAN
				  	 MOST OF THE DATA
			plt.show()
			'''
			# WIPE.wipeTime()


			# fig, (sd, n_sd, n) = plt.subplots(3, sharex=True, sharey=False)
			# # plt.tight_layout()
			# plt.xlabel("Ratio of B-field to Particle Beam Energy (T/MeV) ", fontsize=18)

			# fontdict = {'fontsize': 18,
			# 			'fontweight': 5,
			# 			}

		# 	for dep_var_name, dep_var_LIST in data.items():

		# 		for dep_var in dep_var_LIST:

		# 			if dep_var_LIST.index(dep_var) == 0:
		# 				label = 'right'
		# 			if dep_var_LIST.index(dep_var) == 1:
		# 				label = 'left'

		# 			if dep_var_name == 'SD':
		# 				sd.plot(be_ratio, dep_var, label=label)	
		# 				sd.set(ylabel=dep_var_name)
		# 				title = dep_var_name + " vs be_ratio (T/MeV)"
		# 				sd.set_title(title, fontdict=fontdict)
		# 				# popt = CF.fit(function, be_ratio, dep_var)
		# 				# sd.plot(be_ratio, function(be_ratio, *popt), label=label)

		# 			if dep_var_name == "n/sd":
		# 				n_sd.plot(be_ratio, dep_var, label=label)	
		# 				n_sd.set(ylabel=dep_var_name)
		# 				title = dep_var_name + " vs be_ratio (T/MeV)"
		# 				n_sd.set_title(title, fontdict=fontdict)


		# 			if dep_var_name == "cluster_size":				
		# 				n.plot(be_ratio, dep_var, label=label)
		# 				n.set(ylabel=dep_var_name)
		# 				title = dep_var_name + " vs be_ratio (T/MeV)"
		# 				n.set_title(title, fontdict=fontdict)

		# 	plt.xticks(tickMarks)
		# 	plt.legend()
		# plt.show()


		#### plotting E vs optimal B/E ####
		# plt.figure()
		# plt.xlabel("Energy (MeV)", fontsize=18)
		# plt.ylabel("Optimal B/E ratio (mT/MeV)", fontsize=18)
		# plt.ylim(0,0.0002*1000)
		# plt.title("Optimal B/E ratio (mT/MeV) vs. e+ Energy (MeV)", fontsize=24)
		# plt.plot(dummy_x, np.multiply(gathered_data_right, 1000), label='right')
		# plt.plot(dummy_x, np.multiply(gathered_data_left, 1000), label='left')
		# # plt.ticklabel_format(axis='both', style='sci', scilimits=(-7,0))
		# plt.legend()
		# plt.show()

	# fit median data
	x = np.arange(0, len(cluster_time_median_LIST))
	y = cluster_time_median_LIST
	function = rational3_3
	popt = CF.fit(function, x, y)

	plt.figure()
	plt.ylabel("Median cluster time (ns)", fontsize=18)
	plt.xlabel("Run number", fontsize=18)
	# plt.plot(x,y)
	plt.plot(np.arange(-60, 60, 1), function(np.arange(-60, 60, 1), *popt))
	plt.show()







