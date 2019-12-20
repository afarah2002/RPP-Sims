#----------imports----------#
from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import time
import thread
import numpy as np
from scipy import optimize
from scipy.stats import linregress
import math

#----file imports--------#
from geom_constructor import GeomConstructor 
# from beam import BeamInitializer
from beam2 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction, MySteppingAction, MyField, Plotter, WipeData


PLT = Plotter()
WIPE = WipeData()

# energy_LIST = list(np.arange(2., 9., 1.)) # MeV
# energy_LIST = list(np.arange(1., 50., 1.)) # MeV
# energy_LIST = list(np.logspace(-6., 3., num=500, endpoint=True, base=10)) # eV
energy_LIST = list(np.arange(1.901e-3, 1e3, 1e-5))
# energy_LIST = [1000.]


# energy_LIST = []

global opt_b_right_LIST
global opt_b_left_LIST
opt_b_right_LIST = []
opt_b_left_LIST = []

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
		# gApplyUICommand("/vis/viewer/zoom 1.00001")
		
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

def reject_outliers(data, m=2):

    return np.array(data)[abs(data - np.mean(data)) < m * np.std(data)]

if __name__ == '__main__':
	print(len(energy_LIST))
	# time.sleep(1)
	data_right = open("data_right.txt", "a")
	data_left = open("data_left.txt", "a")
	while True:
		for e in energy_LIST:
			WIPE.wipeComps()
			energy = e
			# energy = 2.5
			print energy, "\n"

			# time.sleep(1)
			b_step = 1.e-9
			B = np.arange(0, 5.e-2, b_step)
			# B = [0.05]
			# print("B len: ", len(B))
			B = [np.sqrt(e*5.16e-16)]

			for b in B:
				# angle += 10
				# angle += 0.075 # +0.075 is a recommended delta theta
				# set user actions ...
				PGA_1 = MyPrimaryGeneratorAction(energy)
				gRunManager.SetUserAction(PGA_1)

				myEA = MyEventAction()
				gRunManager.SetUserAction(myEA)

				mySA = MySteppingAction()
				gRunManager.SetUserAction(mySA)

				# b = np.sqrt(3*(energy/100)**2)
				# b = e/(np.log(2*e))

				# b = abs(10**np.log(e)*e)
				print "energy: ", e*10e6, "eV", "\n", "B: ", b,"T"

				# b = e*np.e**(np.sqrt(e))/np.log(e)
				# b = (np.log(e)-1)/(np.log(e)**2)
				# vectorList = [list(np.multiply([energy, energy, energy], be))]
				# vectorList = [[b, b, b]]
				# for v in vectorList: 
				fieldMgr = gTransportationManager.GetFieldManager()
				myField = G4UniformMagField(G4ThreeVector(np.sqrt(3*b**2),np.sqrt(3*b**2),np.sqrt(3*b**2)))
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
				std_devs_LIST, n_LIST, n_sd_LIST, cluster_size_LIST = PLT.dataReturner() # for 3D positions
				# print "CLUSTER SIZES", "\n", cluster_size_LIST
				# # print len(n_sd_LIST[0])
				gradient_length = 100
				if e > 1e-3:
					gradient_length = 100
				if e > 1:
					gradient_length = 400
				if e > 1e2:
					gradient_length = 500

				if len(n_sd_LIST[0]) > gradient_length:
					# print "HIIIIII"
					avg_change = np.mean(np.gradient(n_sd_LIST[0][-(gradient_length + 1):-1]))
					print "avg change: ", avg_change, 
					if avg_change < -0.035:
						print "should stop now"
						time.sleep(.5)
						break 
						
						
					

				# PLT.wipeData() #clean lists before starting another run

			for num, n_sd in enumerate(n_sd_LIST):
				if num == 0: # right cluster
					print "NSD RIGHT"
					# time.sleep(1)
					max_n_sd = max(n_sd) # gets the peak n/sd that shows the greatest clustering efficiency
					opt_b = n_sd.index(max_n_sd) * b_step # find the be_ratio that produces that max n/sd 

					opt_b_right_LIST.append(opt_b)
					data_right.write(str(opt_b)+"\n")

				else: # left cluster
					print "NSD LEFT"
					# time.sleep(1)
					max_n_sd = max(n_sd) # gets the peak n/sd that shows the greatest clustering efficiency
					opt_b = n_sd.index(max_n_sd) * b_step # find the be_ratio that produces that max n/sd 

					opt_b_left_LIST.append(opt_b)
					data_left.write(str(opt_b)+"\n")

			function = rational3_3

			data  = {
					"SD" : std_devs_LIST, \
					# "means" : means_LIST, \
					"n/sd" : n_sd_LIST, \
					"cluster_size" : n_LIST \
					}

		# fig, (sd, n_sd, n) = plt.subplots(3, sharex=True, sharey=False)
		# # plt.tight_layout()
		# plt.xlabel("Ratio of B-field to Particle Beam Energy (T/MeV) ", fontsize=18)

		# fontdict = {'fontsize': 18,
		# 			'fontweight': 5,
		# 			}

		# for dep_var_name, dep_var_LIST in data.items():

		# 	for dep_var in dep_var_LIST:

		# 		if dep_var_LIST.index(dep_var) == 0:
		# 			label = 'right'
		# 		if dep_var_LIST.index(dep_var) == 1:
		# 			label = 'left'

		# 		if dep_var_name == 'SD':
		# 			sd.plot(be_ratio, dep_var, label=label)	
		# 			sd.set(ylabel=dep_var_name)
		# 			title = dep_var_name + " vs be_ratio (T/MeV)"
		# 			sd.set_title(title, fontdict=fontdict)
		# 			# popt = CF.fit(function, be_ratio, dep_var)
		# 			# sd.plot(be_ratio, function(be_ratio, *popt), label=label)

		# 		if dep_var_name == "n/sd":
		# 			n_sd.plot(be_ratio, dep_var, label=label)	
		# 			n_sd.set(ylabel=dep_var_name)
		# 			title = dep_var_name + " vs be_ratio (T/MeV)"
		# 			n_sd.set_title(title, fontdict=fontdict)


		# 		if dep_var_name == "cluster_size":				
		# 			n.plot(be_ratio, dep_var, label=label)
		# 			n.set(ylabel=dep_var_name)
		# 			title = dep_var_name + " vs be_ratio (T/MeV)"
		# 			n.set_title(title, fontdict=fontdict)

		# plt.xticks(tickMarks)
		# plt.legend()
		# plt.show()
	# plt.figure()
	# plt.xlabel("Energy (MeV)", fontsize=18)
	# plt.ylabel("Optimal B/E ratio (mT/MeV)", fontsize=18)
	# plt.ylim(0,0.0002*1000)
	# plt.title("Optimal B/E ratio (mT/MeV) vs. e+ Energy (MeV)", fontsize=24)
	# # right = reject_outliers(gathered_data_right)
	# # for i in right:
	# # 	data_right.write(str(i) + "\n")
	# # left = reject_outliers(gathered_data_left)
	# # for i in left:
	# # 	data_left.write(str(i) + "\n")
	# plt.plot(dummy_x, np.multiply(gathered_data_right, 1000), label='right')
	# plt.plot(dummy_x, np.multiply(gathered_data_left, 1000), label='left')
	# # plt.ticklabel_format(axis='both', style='sci', scilimits=(-7,0))
	# plt.legend()
	# plt.show()









