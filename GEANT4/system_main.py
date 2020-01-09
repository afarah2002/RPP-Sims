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

# energy_LIST = [1001]
# energy_LIST = list(np.arange(1, 4, 0.5))
energy_LIST = list(np.arange(1001., 2000., 200.))
global energyUnit
energyUnit = eV

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
class ClusterClass(object):

	def __init__(self):

		self.ALL_clusters_positions = []
		self.ALL_clusters_momenta = []

		self.viz_theta = 35
		self.viz_phi = 35
	def run(self, energy, location_range, particleCount, cluster_width, edge):

		if energyUnit == MeV:
			constant = 4.644e-9
		if energyUnit == keV:
			constant = 4.644e-12
		if energyUnit == eV:
			constant = 4.644e-15

		b = np.sqrt(energy*constant) 

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
				magVec, magVecScaled = FD.spherefieldParam(energy, b, phi, theta, cluster_width)

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

			# saving ALL clusters for this energy to be SEEPed
			cluster_positions, cluster_momenta = DA.clusterDataReturner()
			self.ALL_clusters_positions.append(cluster_positions) 
			self.ALL_clusters_momenta.append(cluster_momenta) 

		return self.ALL_clusters_positions, self.ALL_clusters_momenta

CC = ClusterClass()
if __name__ == '__main__':
	step = 6
	particleCount = 500
	edge = 500
	cluster_width = 190
	# location_range = ClusGen.sphericalClusters(step)
	cluster_results = {}
	location_range = ClusGen.cartesianClusters(cluster_width, edge)
	for e in energy_LIST:
		CC
		positions, momenta = CC.run(e, location_range, particleCount, cluster_width, edge)

		cluster_results[e] = [positions, momenta]
	gApplyUICommand("/vis/viewer/zoom 1.5")
	for e, data in cluster_results.items():
		positions, momenta = data[0], data[1]
		for pos in positions:
			mom = momenta[positions.index(pos)]
			# print mom
			SEEP.runSEE(e, pos, mom)
			# time.sleep(3)
		# gRunManager.GeometryHasBeenModified()

