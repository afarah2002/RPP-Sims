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
#----file imports--------#
from geom_constructor import GeomConstructor 
# from beam import BeamInitializer
from beam2 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction, MySteppingAction, MyField, Plotter, WipeData


## setting up lists for std devs x,y,z, pos/neg
global std_devs_LIST
std_devs_LIST  = []

global std_dev_pos_x_right_LIST
global std_dev_pos_x_left_LIST
global std_dev_pos_y_right_LIST
global std_dev_pos_y_left_LIST
global std_dev_pos_z_right_LIST
global std_dev_pos_z_left_LIST

std_dev_pos_x_right_LIST = []
std_dev_pos_x_left_LIST = []
std_dev_pos_y_right_LIST = []
std_dev_pos_y_left_LIST = []
std_dev_pos_z_right_LIST = []
std_dev_pos_z_left_LIST = []

## setting up lists for means x,y,z, pos/neg
global mean_pos_x_right_LIST
global mean_pos_x_left_LIST
global mean_pos_y_right_LIST
global mean_pos_y_left_LIST
global mean_pos_z_right_LIST
global mean_pos_z_left_LIST

mean_pos_x_right_LIST = []
mean_pos_x_left_LIST = []
mean_pos_y_right_LIST = []
mean_pos_y_left_LIST = []
mean_pos_z_right_LIST = []
mean_pos_z_left_LIST = []

## setting up lists for cluster sizes x,y,z, pos/neg
global n_pos_x_right_LIST
global n_pos_x_left_LIST
global n_pos_y_right_LIST
global n_pos_y_left_LIST
global n_pos_z_right_LIST
global n_pos_z_left_LIST

n_pos_x_right_LIST = []
n_pos_x_left_LIST = []
n_pos_y_right_LIST = []
n_pos_y_left_LIST = []
n_pos_z_right_LIST = []
n_pos_z_left_LIST = []

PLT = Plotter()
WIPE = WipeData()

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

if __name__ == '__main__':
	Constructor = Constructor()
	Constructor.construct()
	VIS = Visualizer()

	initialMomenta = []
	finalMomenta = []

	angle = 35
	zoom = 1.5
	# energy = 2.5
	tickMarks = np.arange(1e-7, 5e-1, 0.05)
	be_ratio = np.arange(1e-7, 5e-1, .005) # ratio between magnetic field (varied) and particle energy (fixed @ 2.5 MeV)
	# be_ratio = [0.0001] # after all tests, used to verify best be_ratio, should display 3D position plot
	bound_lower, bound_upper, vectorCount, energy = PLT.paramReturner()
	print energy, "\n"

	print(len(be_ratio))
	time.sleep(1)

	for be in be_ratio:
		# angle += 0.075 # +0.075 is a recommended delta theta
		# for be in be_ratio:
		print "\n", "NEW BE_RATIO: ", be, "\n"
		# time.sleep(1)

		PLT # initialize the x y z lists for the plotter
		# set user actions ...
		PGA_1 = MyPrimaryGeneratorAction()
		gRunManager.SetUserAction(PGA_1)

		myEA = MyEventAction()
		gRunManager.SetUserAction(myEA)

		mySA = MySteppingAction()
		gRunManager.SetUserAction(mySA)


		vectorList = [
				# [1., 1., 1.], 
			 	# [10., 10., 10.]
			 	list(np.multiply([energy, energy, energy], be))
			 	# list(np.multiply([0, 0.1, 0.1], be))
			 	# list(np.multiply([0.1, 0.1, 0.1], be))
			 	# list(np.multiply([0.1, 0.1, 0.1], be))
			 	# list(np.multiply([0.1, 0.1, 0.1], be))
			 	# list(np.multiply([0.1, 0.1, 0.1], be))

			 	# [0., 0., 1]
			 	# [0,0,0]
			 ]
		for v in vectorList: 
			fieldMgr = gTransportationManager.GetFieldManager()
			myField = G4UniformMagField(G4ThreeVector(v[0],v[1],v[2]))
			# myField = MyField(1)
			fieldMgr.SetDetectorField(myField)
			fieldMgr.CreateChordFinder(myField)
			# print "|B-field| = ", vectorList[0][0]

		myRA = MyRunAction()
		gRunManager.SetUserAction(myRA)

		gRunManager.Initialize()

		gRunManager.BeamOn(1)

		VIS.visualizer(angle)

		std_devs_LIST, means_LIST, n_LIST, n_sd_LIST = PLT.dataReturner()
		# PLT.wipeData() #clean lists before starting another run


	function = rational3_3

	data  = {
			"SD" : std_devs_LIST, \
			# "means" : means_LIST, \
			"n/sd" : n_sd_LIST, \
			"cluster_size" : n_LIST \
			}

	fig, (sd, n_sd, n) = plt.subplots(3, sharex=True, sharey=False)
	# plt.tight_layout()
	plt.xlabel("Ratio of B-field to Particle Beam Energy (T/MeV) ", fontsize=9)

	for dep_var_name, dep_var_LIST in data.items():
		for dep_var in dep_var_LIST:
			if dep_var_LIST.index(dep_var) == 0:
				label = 'pos_x_right'
			if dep_var_LIST.index(dep_var) == 1:
				label = 'pos_x_left'
			if dep_var_LIST.index(dep_var) == 2:
				label = 'pos_y_right'
			if dep_var_LIST.index(dep_var) == 3:
				label = 'pos_y_left'
			if dep_var_LIST.index(dep_var) == 4:
				label = 'pos_z_right'
			if dep_var_LIST.index(dep_var) == 5:
				label = 'pos_z_left'

			if dep_var_name == 'SD':
				sd.plot(be_ratio, dep_var, label=label)	
				sd.set(ylabel=dep_var_name)
				title = dep_var_name + " vs be_ratio (T/MeV)"
				sd.set_title(title)
				# popt = CF.fit(function, be_ratio, dep_var)
				# sd.plot(be_ratio, function(be_ratio, *popt), label=label)

			if dep_var_name == "n/sd":
				n_sd.plot(be_ratio, dep_var, label=label)	
				n_sd.set(ylabel=dep_var_name)
				title = dep_var_name + " vs be_ratio (T/MeV)"
				n_sd.set_title(title)

			if dep_var_name == "cluster_size":				
				n.plot(be_ratio, dep_var, label=label)
				n.set(ylabel=dep_var_name)
				title = dep_var_name + " vs be_ratio (T/MeV)"
				n.set_title(title)

	plt.xticks(tickMarks)
	plt.legend()

		# csfont = {'fontname':'Times New Roman'}
		# matplotlib.rcParams["font.family"] = "Times New Roman"

	# pp = PdfPages("RESULTS/" + str(bound_lower) + "_" + str(bound_upper) + ".pdf")
	# pp.savefig(fig)
	# pp.close()

	plt.show()











