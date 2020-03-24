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
plt.rc('font',family='Times New Roman')
import numpy as np
import random
import thread
import time
import scipy.stats as ss
from scipy import optimize, stats
import sys
#-----------FILE imports------------#
from beam2_1 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction,MySteppingAction, \
				    DataAnalysis
from electron_emission import SecondaryElectronEmissionProcess
from geom_constructor import GeomConstructor 
from visualizer import Visualizer
from field_designer import FieldDesign
from cluster_generator import ClusterGenerator
#----------code starts here!----------#

# energy_LIST = [100]
energy_LIST = list(np.arange(0, 1e6, 1000))
# energy_LIST = list(np.arange(1001., 2000., 200.))
global energyUnit
energyUnit = eV

GC = GeomConstructor()
class SpaceConstructor(object):
	def __init__(self):
		g4py.NISTmaterials.Construct()
		g4py.ezgeom.Construct()
		exN03PL = g4py.EMSTDpl.PhysicsListEMstd()
		gRunManager.SetUserInitialization(exN03PL)
		air= gNistManager.FindOrBuildMaterial("G4_AIR")

###############################################################################################################################################
# -----------CLASS ASSIGNMENTS----------- #
FD = FieldDesign()
SpaceConst = SpaceConstructor()
VIS = Visualizer()
SEEP = SecondaryElectronEmissionProcess()
DA = DataAnalysis()
ClusGen = ClusterGenerator()
# --------------------------------------- #
class ClusterClass(object):

	def __init__(self):

		self.ALL_clusters_positions = []
		self.ALL_clusters_momenta = []

		self.viz_theta = 90
		self.viz_phi = 0
	def run(self, energy, location_range, particleCount, cluster_width, edge):

		print energy

		if energyUnit == MeV:
			cluster_constant = 4.644e-9
		if energyUnit == keV:
			cluster_constant = 4.644e-12
		if energyUnit == eV:
			cluster_constant = 4.644e-15

		b = np.sqrt(energy*cluster_constant) 

		for location in location_range:

			SpaceConst
			# gRunManager.GetRunManager().GeometryHasBeenModified()

			if len(location) == 3: # this is cartesian
				x = location[0]
				y = location[1]
				z = location[2]
				magVec, magVecScaled = FD.cartesianfieldParam(energy, b, x, y ,z, cluster_width, edge)
			if len(location) == 2: # this is spherical
				phi = location[0]
				theta = location[1]
				magVec, magVecScaled = FD.spherefieldParam(energy, b, phi, theta, cluster_width, edge)

			print "energy = ", energy, " eV"
			# set user actions ...
			PGA_1 = MyPrimaryGeneratorAction(e, energyUnit, magVecScaled, particleCount)
			gRunManager.SetUserAction(PGA_1)
			myEA = MyEventAction()
			gRunManager.SetUserAction(myEA)
			mySA = MySteppingAction()
			gRunManager.SetUserAction(mySA)
			VIS.visualizer(self.viz_theta, self.viz_phi, "cluster_gen")
			fieldMgr = gTransportationManager.GetFieldManager()
			myField = G4UniformMagField(G4ThreeVector(magVec[0], magVec[1], magVec[2]))
			fieldMgr.SetDetectorField(myField)
			fieldMgr.CreateChordFinder(myField)
			myRA = MyRunAction()
			gRunManager.SetUserAction(myRA)
			gRunManager.Initialize()
			gRunManager.BeamOn(1)
			# DA.grapher()

			# DA.timeAnalysis()
			# cluster_time_LIST, cluster_sizes_LIST = DA.dataReturner()

			# saving ALL clusters for this energy to be SEEPed
			# cluster_time, cluster_positions, cluster_momenta = DA.clusterDataReturner()
			# cluster_time, cluster_positions, cluster_momenta = 0,0,0
			# print "time to cluster = ", cluster_time 

			# self.ALL_clusters_positions.append(cluster_positions) 
			# self.ALL_clusters_momenta.append(cluster_momenta) 

		return self.ALL_clusters_positions, self.ALL_clusters_momenta, location

CC = ClusterClass()
if __name__ == '__main__':
	step = 2
	particleCount = 500
	edge = 500
	cluster_width = 188
	# location_range = ClusGen.sphericalClusters(step)
	location_range = ClusGen.cartesianClusters(cluster_width, edge)

	# time testing
	location_range = [[0,0]]
	cluster_results = {}

	for e in energy_LIST:
		CC
		positions, momenta, location = CC.run(e, location_range, particleCount, cluster_width, edge)

		cluster_results[e] = [positions, momenta, location]
	# plt.show()
	gApplyUICommand("/vis/viewer/zoom 1.5")
	# for e, data in cluster_results.items():
	# 	positions, momenta, location = data[0], data[1], data[2]
	# 	for pos in positions:
	# 		mom = momenta[positions.index(pos)]
	# 		SEEP.runSEE(e, pos, mom, location)

