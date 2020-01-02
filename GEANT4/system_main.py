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
from beam2_1 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction,MySteppingAction, \
				    DataAnalysis
from electron_emission import SecondaryElectronEmissionProcess
from geom_constructor import GeomConstructor 
from visualizer import Visualizer
from field_designer import FieldDesign
from cluster_generator import ClusterGenerator
#----------code starts here!----------#



global cluster_width
cluster_width = 160
spacing = 40./7
interval = cluster_width + spacing
interval *= 3
global edge
edge = 500
# energy_LIST = [0.002]
# energy_LIST = list(np.arange(1, 4, 0.5))
energy_LIST = list(np.arange(1001, 2000, 200))
global energyUnit
energyUnit = eV

cluster_time_median_LIST = []


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
GC = GeomConstructor()
ClusGen = ClusterGenerator()
# --------------------------------------- #

viz_theta = 35
viz_phi = 35
zoom = 1.5

class ClusterClass(object):

	def __init__(self):

		self.ALL_clusters_positions = []
		self.ALL_clusters_momenta = []


	def run(self, energy, location_range, particleCount):

		if energyUnit == MeV:
			constant = 4.644e-9
		if energyUnit == keV:
			constant = 4.644e-12
		if energyUnit == eV:
			constant = 4.644e-15

		b = np.sqrt(energy*constant) 

		for location in location_range:
			SpaceConst
			if len(location) == 3: # this is cartesian
				x = location[0]
				y = location[1]
				z = location[2]
				magVec, magVecScaled = FD.cartesianfieldParam(energy, b, x, y ,z, cluster_width)
			if len(location) == 2: # this is spherical
				phi = location[0]
				theta = location[1]
				magVec, magVecScaled = FD.spherefieldParam(energy, b, phi, theta, cluster_width)

			# set user actions ...
			PGA_1 = MyPrimaryGeneratorAction(e, energyUnit, magVecScaled, particleCount)
			gRunManager.SetUserAction(PGA_1)
			myEA = MyEventAction()
			gRunManager.SetUserAction(myEA)
			mySA = MySteppingAction()
			gRunManager.SetUserAction(mySA)
			VIS.visualizer(viz_theta, viz_phi)
			fieldMgr = gTransportationManager.GetFieldManager()
			myField = G4UniformMagField(G4ThreeVector(magVec[0], magVec[1], magVec[2]))
			fieldMgr.SetDetectorField(myField)
			fieldMgr.CreateChordFinder(myField)
			myRA = MyRunAction()
			gRunManager.SetUserAction(myRA)
			gRunManager.Initialize()
			gRunManager.BeamOn(1)

			# print len(cluster_positions), "\n", len(cluster_momenta)
			# SEEP.runSEE(e, cluster_positions, cluster_momenta)
			# saving ALL clusters for this energy to be SEEPed
			cluster_positions, cluster_momenta = DA.clusterDataReturner()
			self.ALL_clusters_positions.append(cluster_positions) 
			# print cluster_positions
			# print "appended", ALL_clusters_positions
			# time.sleep(12)
			self.ALL_clusters_momenta.append(cluster_momenta) 
			# cluster_time_LIST, cluster_size_LIST = DA.dataReturner() # for 3D positions

			# median = np.median(cluster_time_LIST)
			# cluster_time_median_LIST.append(median)
							
		for i in self.ALL_clusters_positions:
			print i[-1]
			time.sleep(2)

		return self.ALL_clusters_positions, self.ALL_clusters_momenta

CC = ClusterClass()
if __name__ == '__main__':
	step = 6
	particleCount = 500
	location_range = ClusGen.sphericalClusters(step)
	for e in energy_LIST:
		CC
		positions, momenta = CC.run(e, location_range, particleCount)
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
