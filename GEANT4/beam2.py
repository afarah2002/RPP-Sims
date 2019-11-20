#----------imports----------#
from Geant4 import * 
import random
import numpy as np
import scipy.stats as ss
import pandas as pd
import seaborn as sns  # for nicer graphics

## matplotlib stuff
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import sleep
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import time


global bound_lower
global bound_upper

bound_lower = 400 # only consider the cluster within these bounds, range of 100
bound_upper = bound_lower + 200
# bound_upper = 550

global vectorCount
vectorCount = 500 # number of scattered e+ per run

# global energy
# energy = SP.sendEnergy()

global pos_3D_right
pos_3D_right = []
global n_sd_pos_3D_right_LIST
sd_pos_3D_right_LIST = []
global n_pos_3D_right_LIST
n_pos_3D_right_LIST = []
global n_sd_pos_3D_right_LIST
n_sd_pos_3D_right_LIST = []

global pos_3D_left
pos_3D_left = []
global n_sd_pos_3D_left_LIST
sd_pos_3D_left_LIST = []
global n_pos_3D_left_LIST
n_pos_3D_left_LIST = []
global n_sd_pos_3D_left_LIST
n_sd_pos_3D_left_LIST = []

global cluster_time_LIST
cluster_time_LIST = []

global px
global py
global pz

px = []
py = []
pz = []



#----------code starts here!----------#


class WipeData(object):
	# wipe lists for next data collection
	def wipe(self):
		pos_3D_right[:] = []
		pos_3D_left[:] = []
		px[:] = []
		py[:] = []
		pz[:] = []

	def wipeTime(self):
		cluster_time_LIST[:] = []

	def wipeComps(self):

		sd_pos_3D_right_LIST[:] = []
		n_pos_3D_right_LIST[:] = []
		n_sd_pos_3D_right_LIST[:] = []
		sd_pos_3D_left_LIST[:] = []
		n_pos_3D_left_LIST[:] = []
		n_sd_pos_3D_left_LIST[:] = []



WIPE = WipeData()

class Plotter(object):
	"graphs 3D positions"

	def __init__(self):
		pass

	def dataCollection(self, posf, momf, tcluster):

		# time.sleep(1)
		# print px, "\n", py, "\n", pz, "\n"
		# time.sleep(1)
		px.append(posf[0])
		py.append(posf[1])
		pz.append(posf[2])


		radius = np.sqrt(np.square(posf[0]) + np.square(posf[1]) + np.square(posf[2]))
		# if tcluster > 25: # weird set of outlier <---- INTERESTING PHENOMENON
		# 	print radius

		# if radius < bound_lower:
		# 	print tcluster

		# if posf[0] < 0 and posf[1] < 0 and posf[2] < 0 and -radius < -bound_lower and -radius > -bound_upper:
		# 	pos_3D_left.append(radius)
		# 	cluster_time_LIST.append(tcluster)
		# if posf[0] > 0 and posf[1] > 0 and posf[2] > 0 and radius > bound_lower and radius < bound_upper:
		# 	pos_3D_right.append(radius)
		# 	cluster_time_LIST.append(tcluster)

		radius_lower = np.sqrt(3 * np.square(bound_lower))
		radius_upper = np.sqrt(3 * np.square(bound_upper))


		if radius > radius_lower and radius < radius_upper:
			cluster_time_LIST.append(tcluster)
			print tcluster


		# print("DATA STORED")
		# print len(self.px), "\n", len(self.py), "\n", len(self.pz), "\n"

	def dataAnalysis(self):
		# results = open("RESULTS/results_10212019_1.txt", "a")

		n_pos_3D_right = len(pos_3D_right)
		mean_pos_3D_right = np.mean(pos_3D_right)
		self.std_dev_pos_3D_right = np.std(pos_3D_right)
		if self.std_dev_pos_3D_right == 0:
			self.std_dev_pos_3D_right = 0.0001
		median_pos_3D_right = np.median(pos_3D_right)
		n_sd_pos_3D_right = float(n_pos_3D_right / self.std_dev_pos_3D_right)

		n_pos_3D_left = len(pos_3D_left)
		mean_pos_3D_left = np.mean(pos_3D_left)
		self.std_dev_pos_3D_left = np.std(pos_3D_left)
		if self.std_dev_pos_3D_left == 0:
			self.std_dev_pos_3D_left = 0.0001
		median_pos_3D_left = np.median(pos_3D_left)
		n_sd_pos_3D_left = float(n_pos_3D_left / self.std_dev_pos_3D_left)

		sd_pos_3D_right_LIST.append(self.std_dev_pos_3D_right)
		n_sd_pos_3D_right_LIST.append(n_sd_pos_3D_right)
		n_pos_3D_right_LIST.append(n_pos_3D_right)

		sd_pos_3D_left_LIST.append(self.std_dev_pos_3D_left)
		n_sd_pos_3D_left_LIST.append(n_sd_pos_3D_left)
		n_pos_3D_left_LIST.append(n_pos_3D_left)



	def dataReturner(self):

		return [sd_pos_3D_right_LIST, sd_pos_3D_left_LIST], \
			   [n_pos_3D_right_LIST, n_pos_3D_left_LIST], \
			   [n_sd_pos_3D_right_LIST, n_sd_pos_3D_left_LIST], \
			   cluster_time_LIST



		pass

	def grapher(self):
		fig = plt.figure()
		# Axes3D.scatter(self.px, self.py, self.pz)

		# first subplot: a 3D scatter plot of positions
		ax = fig.add_subplot(111, projection='3d')
		axmin = -600 # what units are these???????
		axmax = 600
		axes = plt.gca()
		axes.set_xlim([axmin,axmax])
		axes.set_ylim([axmin,axmax])
		axes.set_zlim([axmin,axmax])

		ax.set_xlabel('X position units')
		ax.set_ylabel('Y position units')
		ax.set_zlabel('Z position units')

		ax.scatter(px, py, pz)
		plt.title("3D Positions of Randomly Scattered e+")
		plt.show()

	def paramReturner(self):
		return bound_lower, bound_upper, vectorCount, energy


		 


PLT = Plotter()

class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self,energy):
		G4VUserPrimaryGeneratorAction.__init__(self)
		self.particleGun = G4ParticleGun(1)
		# print("\n Particle gun defined \n")
		self.energy = energy
	def GeneratePrimaries(self, event):


		# Particle param
		#################################################
		locationArray = [0, 0, 0]

		particle = "e+"
		# energy_2 = 2.5
		energyUnit = MeV 
		dimensionUnit = cm

		energy = self.energy
		self.particleGun.SetParticleByName(particle) # define particle
		self.particleGun.SetParticleEnergy(energy*energyUnit) # define particle energy 

		for i in range(0, vectorCount): # creates random momentum vectors originating from [0, 0, 0]
			mx = random.uniform(-1,1)
			my = random.uniform(-1,1)
			mz = random.uniform(-1,1)
			momentumArray = [mx, my, mz]
			self.particleGun.SetParticlePosition(G4ThreeVector(locationArray[0], locationArray[1], locationArray[2])*dimensionUnit) # define first particle generator location
			self.particleGun.SetParticleMomentumDirection(G4ThreeVector(momentumArray[0], momentumArray[1], momentumArray[2])*dimensionUnit) # define first particle generator momentum
			self.particleGun.GeneratePrimaryVertex(event)
		#################################################

#-------------------------------------------------------------------
class MyRunAction(G4UserRunAction):
	"My Run Action"

	def EndOfRunAction(self, run):
		PLT.grapher()
		PLT.dataAnalysis()
		WIPE.wipe()
		# print "*** End of Run"
		# print "- Run sammary : (id= %d, #events= %d)" \
		# % (run.GetRunID(), run.GetNumberOfEventToBeProcessed())

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
	"My Event Action"

	def EndOfEventAction(self, event):
		pass

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
	"My Stepping Action"

	def UserSteppingAction(self, step):
		preStepPoint = step.GetPreStepPoint()
		postStepPoint = step.GetPostStepPoint()

		track = step.GetTrack()
		touchable = track.GetTouchable()
		KE = track.GetKineticEnergy()

		# kinetic energy in MeV - PRE
		# initialKE = preStepPoint.GetKineticEnergy() 
		# kinetic energy in MeV - POST
		# finalKE = postStepPoint.GetKineticEnergy()

		m = [track.GetMomentum().x, track.GetMomentum().y, track.GetMomentum().z] # equal to the postStepPoint momentum
		p = [postStepPoint.GetPosition().x, postStepPoint.GetPosition().y, postStepPoint.GetPosition().z] # (mm)
		mm = np.sqrt((m[0])**2 + (m[1])**2 + (m[2])**2)
		t = track.GetGlobalTime() # (ns)

		# momenta - PRE
		initialMomentum = [preStepPoint.GetMomentum().x, preStepPoint.GetMomentum().y, preStepPoint.GetMomentum().z]
		# momenta - POST
		finalMomentum = [postStepPoint.GetMomentum().x, postStepPoint.GetMomentum().y, postStepPoint.GetMomentum().z]

		# print KE, "\n", p, "\n", initialMomentum, "\n", finalMomentum, "\n\n" 
		# energy = step.GetTotalEnergyDeposit()


		PLT.dataCollection(p, m, t) # calls data collection and analysis on final positions and momenta
		# return initialMomentum, finalMomentum 

class MyField(G4MagneticField): ### used when mag field NOT parameterized in main filed
	"My Magnetic Field"

	def __init__(self, eb_ratio):
		self.eb_ratio = eb_ratio

	def GetFieldValue(self, pos, time):
		self.eb_ratio = 1
		vectorList = [
						# [1., 1., 1.], 
					 	# [10., 10., 10.]
					 	list(np.multiply([0.1, 0.1, 0.1], self.eb_ratio))
					 	# [0., 0., 1]
					 	# [0,0,0]
					 ]
		for v in vectorList:

			bfield = G4ThreeVector()
			bfield.x = v[0]*tesla
			bfield.y = v[1]*tesla
			bfield.z = v[2]*tesla
			# print "\n", "B-field activated", "\n" ### gets rid of other prints for some reason
			return bfield


