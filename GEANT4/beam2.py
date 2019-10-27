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

global bound_lower
global bound_upper

bound_lower = 0 # only consider the cluster within these bounds, range of 100
bound_upper = bound_lower + 100

global vectorCount
vectorCount = 100 # number of scattered e+ per run

global energy
energy_1 = 2.5

## setting up lists for std devs x,y,z, pos/neg
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

#----------code starts here!----------#
class Plotter(object):
	"graphs 3D positions"

	def __init__(self):
		self.px = []
		self.py = []
		self.pz = []

		self.px_pos = []
		self.px_neg = []
		self.py_pos = []
		self.py_neg = []
		self.pz_pos = []
		self.pz_neg = []
	
	def dataCollection(self, posf, momf):
		self.px.append(posf[0])
		self.py.append(posf[1])
		self.pz.append(posf[2])
		self.pos3D = np.sqrt(np.square(self.px) + np.square(self.py) + np.square(self.pz)) # 3D position

		# cutoff = 400 # above/below this value (+ or -) the cluster is analyzed
		# bound_upper = 300
		# bound_lower = 200 
		#isolate clusters

		# if posf[0] > cutoff:
		if posf[0] > bound_lower and posf[0] < bound_upper:
			self.px_pos.append(posf[0])
		# if posf[0] < -cutoff:
		if posf[0] < -bound_lower and posf[0] > -bound_upper:
			self.px_neg.append(posf[0])
		# if posf[1] > cutoff:
		if posf[1] > bound_lower and posf[1] < bound_upper:
			self.py_pos.append(posf[1])
		# if posf[1] < -cutoff:
		if posf[1] < -bound_lower and posf[1] > -bound_upper:
			self.py_neg.append(posf[1])
		# if posf[2] > cutoff:
		if posf[2] > bound_lower and posf[2] < bound_upper:
			self.pz_pos.append(posf[2])
		# if posf[2] < -cutoff:
		if posf[2] < -bound_lower and posf[2] > -bound_upper:
			self.pz_neg.append(posf[2])

		print("DATA STORED")
		# print len(self.px), "\n", len(self.py), "\n", len(self.pz), "\n"

	def dataAnalysis(self):
		# results = open("RESULTS/results_10212019_1.txt", "a")

		# POSITION X - ALL
		# mean_pos_ = 
		# std_dev_pos_ = 
		# median_pos_ = 		
		# # POSITION Y - ALL
		# mean_pos_ = 
		# std_dev_pos_ = 
		# median_pos_ = 
		# # POSITION Z - ALL
		# mean_pos_ = 
		# std_dev_pos_ = 
		# median_pos_ = 
		# POSITION X - POS(RIGHT)
		n_pos_x_right = len(self.px_pos)
		mean_pos_x_right = np.mean(self.px_pos)
		self.std_dev_pos_x_right = np.std(self.px_pos)
		median_pos_x_right = np.median(self.px_pos)
		# results.write(str(mean_pos_x_right) +  "	" +  str(self.std_dev_pos_x_right) + "	" + str(median_pos_x_right) + "\n")
		# POSITION X - NEG(LEFT)
		n_pos_x_left = len(self.px_neg)
		mean_pos_x_left = np.mean(self.px_neg)
		self.std_dev_pos_x_left = np.std(self.px_neg)
		median_pos_x_left = np.median(self.px_neg)
		# results.write(str(mean_pos_x_left) +  " 	" +  str(self.std_dev_pos_x_left) + "	" + str(median_pos_x_left) + "\n")
		# POSITION Y - POS(RIGHT)
		n_pos_y_right = len(self.py_pos)
		mean_pos_y_right = np.mean(self.py_pos)
		self.std_dev_pos_y_right = np.std(self.py_pos)
		median_pos_y_right = np.median(self.py_pos)
		# results.write(str(mean_pos_y_right) +  "	" +  str(self.std_dev_pos_y_right) + "	" + str(median_pos_y_right) + "\n")
		# POSITION Y - NEG(LEFT)
		n_pos_y_left = len(self.py_neg)
		mean_pos_y_left = np.mean(self.py_neg)
		self.std_dev_pos_y_left = np.std(self.py_neg)
		median_pos_y_left = np.median(self.py_neg)
		# results.write(str(mean_pos_y_left) +  "    " +  str(self.std_dev_pos_y_left) + "	" + str(median_pos_y_left) + "\n")
		# POSITION Z - POS(RIGHT)
		n_pos_z_right = len(self.pz_pos)
		mean_pos_z_right = np.mean(self.pz_pos)
		self.std_dev_pos_z_right = np.std(self.pz_pos)
		median_pos_z_right = np.median(self.pz_pos)
		# results.write(str(mean_pos_z_right) +  "	" +  str(self.std_dev_pos_z_right) + "	" + str(median_pos_z_right) + "\n")
		# POSITION Z - NEG(LEFT)
		n_pos_z_left = len(self.pz_neg)
		mean_pos_z_left = np.mean(self.pz_neg)
		self.std_dev_pos_z_left = np.std(self.pz_neg)
		median_pos_z_left = np.median(self.pz_neg)
		# results.write(str(mean_pos_z_left) +  "	  " +  str(self.std_dev_pos_z_left) + "	" + str(median_pos_z_left) + "\n\n")
		# print self.std_dev_pos_x_right, self.std_dev_pos_x_left, self.std_dev_pos_y_right, self.std_dev_pos_y_left, self.std_dev_pos_z_right, self.std_dev_pos_z_left

		std_dev_pos_x_right_LIST.append(self.std_dev_pos_x_right)
		std_dev_pos_x_left_LIST.append(self.std_dev_pos_x_left)
		std_dev_pos_y_right_LIST.append(self.std_dev_pos_y_right)
		std_dev_pos_y_left_LIST.append(self.std_dev_pos_y_left)
		std_dev_pos_z_right_LIST.append(self.std_dev_pos_z_right)
		std_dev_pos_z_left_LIST.append(self.std_dev_pos_z_left)

		mean_pos_x_right_LIST.append(mean_pos_x_right)
		mean_pos_x_left_LIST.append(mean_pos_x_left)
		mean_pos_y_right_LIST.append(mean_pos_y_right)
		mean_pos_y_left_LIST.append(mean_pos_y_left)
		mean_pos_z_right_LIST.append(mean_pos_z_right)
		mean_pos_z_left_LIST.append(mean_pos_z_left)

		n_pos_x_right_LIST.append(n_pos_x_right)
		n_pos_x_left_LIST.append(n_pos_x_left)
		n_pos_y_right_LIST.append(n_pos_y_right)
		n_pos_y_left_LIST.append(n_pos_y_left)
		n_pos_z_right_LIST.append(n_pos_z_right)
		n_pos_z_left_LIST.append(n_pos_z_left)

	def dataReturner(self):
		return [std_dev_pos_x_right_LIST, std_dev_pos_x_left_LIST, std_dev_pos_y_right_LIST, std_dev_pos_y_left_LIST, std_dev_pos_z_right_LIST, std_dev_pos_z_left_LIST], \
			   [mean_pos_x_right_LIST, mean_pos_x_left_LIST, mean_pos_y_right_LIST, mean_pos_y_left_LIST, mean_pos_z_right_LIST, mean_pos_z_left_LIST], \
			   [n_pos_x_right_LIST, n_pos_x_left_LIST, n_pos_y_right_LIST, n_pos_y_left_LIST, n_pos_z_right_LIST, n_pos_z_left_LIST]
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

		ax.scatter(self.px, self.py, self.pz)

		# second subplot: a histogram of positions
		fig, axs = plt.subplots(11, sharey=False, tight_layout=False)
		n_bins = 25
		axs[0].hist(self.px, bins=n_bins) # histogram of 3D position x
		axs[1].hist(self.py, bins=n_bins) # histogram of 3D position x
		axs[2].hist(self.pz, bins=n_bins) # histogram of 3D position z
		axs[3].hist(self.pos3D, bins=n_bins) # histogram of 3D position magnitude

		axs[5].hist(self.px_pos, bins=n_bins) # histogram of 3D position x
		axs[6].hist(self.px_neg, bins=n_bins) # histogram of 3D position x
		axs[7].hist(self.py_pos, bins=n_bins) # histogram of 3D position x
		axs[8].hist(self.py_neg, bins=n_bins) # histogram of 3D position x
		axs[9].hist(self.pz_pos, bins=n_bins) # histogram of 3D position z
		axs[10].hist(self.pz_neg, bins=n_bins) # histogram of 3D position z
		# plt.figure()
		# n_bins = 1000
		# plt.hist([self.px, self.py, self.pz], bins=n_bins, histtype='barstacked', normed=True)
		# plt.show()
	def paramReturner(self):
		return bound_lower, bound_upper, vectorCount, energy_1

	def wipeData(self):
		# wipe lists for next data collection
		self.px = []
		self.py = []
		self.pz = []

		self.px_pos = []
		self.px_neg = []
		self.py_pos = []
		self.py_neg = []
		self.pz_pos = []
		self.pz_neg = []
 


PLT = Plotter()

class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self):
		G4VUserPrimaryGeneratorAction.__init__(self)
		self.particleGun = G4ParticleGun(1)
		# print("\n Particle gun defined \n")

	def GeneratePrimaries(self, event):


		# Particle param
		#################################################
		locationArray = [0, 0, 0]

		particle = "e+"

		# energy_2 = 2.5
		energyUnit = MeV 
		dimensionUnit = cm

		self.particleGun.SetParticleByName(particle) # define particle
		self.particleGun.SetParticleEnergy(energy_1*energyUnit) # define particle energy 

		for i in range(0,vectorCount): # creates random momentum vectors originating from [0, 0, 0]
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
		# PLT.grapher()
		PLT.dataAnalysis()
		print "*** End of Run"
		print "- Run sammary : (id= %d, #events= %d)" \
		% (run.GetRunID(), run.GetNumberOfEventToBeProcessed())

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
	"My Event Action"

	def EndOfEventAction(self, event):
		#print "*** dE/dx in current step=", step.GetTotalEnergyDeposit()
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
		p = [track.GetPosition().x, track.GetPosition().y, track.GetPosition().z]
		mm = np.sqrt((m[0])**2 + (m[1])**2 + (m[2])**2)

		# momenta - PRE
		initialMomentum = [preStepPoint.GetMomentum().x, preStepPoint.GetMomentum().y, preStepPoint.GetMomentum().z]
		# momenta - POST
		finalMomentum = [postStepPoint.GetMomentum().x, postStepPoint.GetMomentum().y, postStepPoint.GetMomentum().z]

		# print KE, "\n", p, "\n", initialMomentum, "\n", finalMomentum, "\n\n" 
		# energy = step.GetTotalEnergyDeposit()



		PLT.dataCollection(p, m) # calls data collection and analysis on final positions and momenta
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


