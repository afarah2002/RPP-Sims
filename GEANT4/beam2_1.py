#----------GEANT4 IMPORTS----------#
from Geant4 import * 

# --------PYTHON IMPORTS ----------#
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy as np
import pandas as pd
import random
from scipy.signal import find_peaks
import scipy.stats as ss
import seaborn as sns  # for nicer graphics
import time
# -------- FILE IMPORTS ------- #
from arrow_generator import Arrow3D

#----------code starts here!----------#
global times
global cluster_time
times = []

time_file = open("time_analysis_TIMES.txt", "w")
time_means_file = open("time_analysis_TIMES_MEANS", "w")

class DataAnalysis(object):
	"Performs data collection, analysis and visualization for the pre-SEE analysis"

	# for a specific cluster
	avg_cluster_time = 0

	def __init__(self):

		self.p3D = []
		self.px = []
		self.py = []
		self.pz = []

		self.m3D = []
		self.mx = []
		self.my = []
		self.mz = []

		self.cluster_time_LIST = []
		self.cluster_sizes_LIST = []

		times = []


	def wipeData(self):

		self.p3D[:] = []
		self.px[:] = []
		self.py[:] = []
		self.pz[:] = []

		self.m3D[:] = []
		self.mx[:] = []
		self.my[:] = []
		self.mz[:] = []

	def dataCollection(self, posf, momf, tcluster):

		radius = np.sqrt(np.square(posf[0]) + np.square(posf[1]) + np.square(posf[2]))

		# p3D.append(posf)
		self.p3D.append(posf)
		# print len(self.p3D)	
		self.m3D.append(momf)

		self.px.append(posf[0])
		self.py.append(posf[1])
		self.pz.append(posf[2])

		self.mx.append(momf[0])
		self.my.append(momf[1])
		self.mz.append(momf[2])

		bound_lower = 490 # only consider the cluster within these bounds, range of 100
		bound_upper = bound_lower + 400

		radius_lower = np.sqrt(3 * np.square(bound_lower))
		radius_upper = np.sqrt(3 * np.square(bound_upper))

		# if radius > bound_lower and radius < bound_upper:
		self.cluster_time_LIST.append(tcluster)

	def dataReturner(self):
		# print self.cluster_time_LIST
		return self.cluster_time_LIST, self.cluster_sizes_LIST

	def grapher(self):

		fig = plt.figure()
		# Axes3D.scatter(self.self.px, self.self.py, self.self.pz)

		# first subplot: a 3D scatter plot of positions
		ax = fig.add_subplot(111, projection='3d')
		axmin = -600 
		axmax = 600
		axes = plt.gca()
		axes.set_xlim([axmin,axmax])
		axes.set_ylim([axmin,axmax])
		axes.set_zlim([axmin,axmax])

		ax.set_xlabel('mm')
		ax.set_ylabel('mm')
		ax.set_zlabel('mm')

		ax.scatter(self.px, self.py, self.pz)
		scale_fac = 50
		# for i in np.arange(0, len(self.px)):
		# 	a = Arrow3D([self.px[i], self.px[i] + mx[i]*scale_fac], [self.py[i], self.py[i] + my[i]*scale_fac], [self.pz[i], self.pz[i] + mz[i]*scale_fac], mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
		# 	ax.add_artist(a)


		plt.title("3D Positions of Randomly Scattered e+")

		# plt.draw() 
		
		plt.show()

	def computeClusterMomentum(self):

		clusterx = []
		clustery = []
		clusterz = []

		momx = []		
		momy = []
		momz = []

		# print "length of p3D = ", len(p3D)

		for pos in self.p3D:
			difference = np.sqrt((clusterCenter[0] - pos[0])**2+ \
								 (clusterCenter[1] - pos[1])**2+ \
								 (clusterCenter[2] - pos[2])**2)

			pos_index = self.p3D.index(pos)

			if difference < 160:
				clusterx.append(pos[0] - clusterCenter[0])
				clustery.append(pos[1] - clusterCenter[1])
				clusterz.append(pos[2] - clusterCenter[2])

				momx.append(self.m3D[pos_index][0])
				momy.append(self.m3D[pos_index][1])
				momz.append(self.m3D[pos_index][2])

		# print "number of clustered positrons = ", len(clusterx)

		# saving data

		global C_positions_LIST 
		global C_momenta_LIST 
		C_positions_LIST = []
		C_momenta_LIST = []

		for index in range(len(clusterx)):
			position = [clusterx[index], clustery[index], clusterz[index]]
			momentum = [momx[index], momy[index], momz[index]]
			
			# positions.write("["+ str(position[0])+ ", "+ str(position[1])+ ", "+ str(position[2])+ "]"+"\n")
			# momenta.write("["+ str(momentum[0])+ ", "+ str(momentum[1])+ ", "+ str(momentum[2])+ "]"+"\n")

			C_positions_LIST.append(position)
			C_momenta_LIST.append(momentum)


		# print len(C_positions_LIST), "\n", len(C_momenta_LIST)
		# fig = plt.figure()
		# # Axes3D.scatter(self.self.px, self.self.py, self.self.pz)

		# # first subplot: a 3D scatter plot of positions
		# ax = fig.add_subplot(111, projection='3d')
		# axmin = -160 
		# axmax = 160
		# axes = plt.gca()
		# axes.set_xlim([axmin,axmax])
		# axes.set_ylim([axmin,axmax])
		# axes.set_zlim([axmin,axmax])

		# ax.set_xlabel('X (mm)')
		# ax.set_ylabel('Y (mm)')
		# ax.set_zlabel('Z (mm)')

		# ax.scatter(clusterx, clustery, clusterz)
		# for i in np.arange(0, len(clusterx)):
		# 	a = Arrow3D([clusterx[i], clusterx[i] + 5000*momx[i]], [clustery[i], clustery[i] + 1000*momy[i]], [clusterz[i], clusterz[i] + 1000*momz[i]], mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
		# 	ax.add_artist(a)


		# plt.title("A single cluster")

		# # plt.draw() 
		
		# plt.show()


	def computeClusterSize(self):
		n_bins = 25
		position_LIST = [self.px, self.py, self.pz]
		for i in position_LIST:
			# fig2, ax2 = plt.subplots(2, sharey=True, sharex=False, tight_layout=False)
			range_LIST = []
			positive = []
			negative = []
			# counts, bins, bars = axs[position_LIST.index(i)].hist(i, bins=n_bins) 
			# print counts, "\n\n", bins

			# separate x, y, z into positive and negative
			for pos in i: 
				# print pos
				if pos > 0:
					positive.append(pos)
				if pos < 0:
					negative.append(pos)
			frame = 0
			for i in [positive, negative]:
				# print "sign change" , i, "\n"
				# plt.hist(i, bins=n_bins)
				counts, bins = np.histogram(i, bins=n_bins)
				counts = list(counts)
				# print counts, "\n\n", "bins", "\n", bins
				for freq in counts:
					# gets rid of any outliers
					if freq not in np.arange(5, 150): 
						index = counts.index(freq)
						np.delete(bins, index)
				# takes the range of each pos/neg without outliers
				rng = np.max(bins) - np.min(bins)
				range_LIST.append(rng)
				# print range_LIST
				frame += 1
			#combines the pos/neg ranges to get the real range
			true_range = np.sum(range_LIST)
			# print true_range
			if true_range < 200:
				self.cluster_sizes_LIST.append(true_range)
		print "avg cluster size = ", np.mean(self.cluster_sizes_LIST)
				# print "\n", "Range = ", true_range, "\n"

	def timeAnalysis(self):
		times = self.cluster_time_LIST
		# global avg_cluster_time
		median_cluster_time = np.median(times)
		time_file.write(str(median_cluster_time)+"\n")

		mean_cluster_time = np.mean(times)
		time_means_file.write(str(mean_cluster_time)+"\n")




		self.avg_cluster_time = mean_cluster_time
		n_bins = 50
		print "time to cluster = ", self.avg_cluster_time 
		# plt.hist(times, n_bins)
		# plt.xlabel("Time to cluster (ns)")
		# plt.ylabel("Frequency")
		# plt.title("Distribution of cluster times for a specific energy")
		# times[:] = []
		# plt.show()

	def clusterDataReturner(self):
		positions = C_positions_LIST
		momenta = C_momenta_LIST

		return avg_cluster_time, positions, momenta


DA = DataAnalysis()

'''																									     											   '''
   #                         								   ____________   _  ____________ 								                           #
   ########################################################## / ___/ __/ _ | / |/ /_  __/ / / ##########################################################
   ##########################################################/ (_ / _// __ |/    / / / /_  _/ ##########################################################
   ##########################################################\___/___/_/ |_/_/|_/ /_/   /_/   ##########################################################
   #																																				   #
'''																									    											   '''

class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self, energy, energyUnit, center, particleCount):
		G4VUserPrimaryGeneratorAction.__init__(self)
		self.particleGun = G4ParticleGun(1)
		# print("\n Particle gun defined \n")
		self.energy = energy
		self.energyUnit = energyUnit
		self.particleCount = particleCount
		global clusterCenter
		clusterCenter = center

	def GeneratePrimaries(self, event):
		# Particle param
		#################################################
		locationArray = [0, 0, 0]

		particle = "e+"
		# energy_2 = 2.5
		energyUnit = self.energyUnit 
		dimensionUnit = cm

		energy = self.energy
		self.particleGun.SetParticleByName(particle) # define particle
		self.particleGun.SetParticleEnergy(energy*energyUnit) # define particle energy 

		for i in range(0, self.particleCount): # creates random momentum vectors originating from [0, 0, 0]
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
		# DA
		# DA.grapher()
		DA.timeAnalysis()
		# DA.computeClusterMomentum()
		# DA.computeClusterSize()
		# DA.wipeData()
		# WIPE.wipe()
		# print "*** End of Run"
		# print "- Run sammary : (id= %d, #events= %d)" \
		# % (run.GetRunID(), run.GetNumberOfEventToBeProcessed())
		pass

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
	"My Event Action"

	def EndOfEventAction(self, event):
		# print "See count", SEE_count

		pass

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
	"My Stepping Action"

	def UserSteppingAction(self, step):
		preStepPoint = step.GetPreStepPoint()
		postStepPoint = step.GetPostStepPoint()

		# change = step.particleChange()
		track = step.GetTrack()
		parentId = track.GetParentID()
		particleName = track.GetDefinition().GetParticleName() 
		# print particleName
		touchable = track.GetTouchable()
		KE = track.GetKineticEnergy()

		# kinetic energy in MeV - PRE
		# initialKE = preStepPoint.GetKineticEnergy() 
		# kinetic energy in MeV - POST
		# finalKE = postStepPoint.GetKineticEnergy()

		if particleName == 'e+':
			# print "positron"
			p_test = [step.GetDeltaPosition().x,step.GetDeltaPosition().y,step.GetDeltaPosition().z]
			p = [postStepPoint.GetPosition().x, postStepPoint.GetPosition().y, postStepPoint.GetPosition().z] # (mm)
			# p and p_test are the SAME 
			t = step.GetDeltaTime()
			# t = track.GetGlobalTime() # (ns)
			# t and t_test are the SAME

			m = [postStepPoint.GetMomentum().x, postStepPoint.GetMomentum().y, postStepPoint.GetMomentum().z]
			# m = [step.GetDeltaMomentum().x, step.GetDeltaMomentum().y, step.GetDeltaMomentum().z] # equal to the postStepPoint momentum
			mm = np.sqrt((m[0])**2 + (m[1])**2 + (m[2])**2)
			# momenta - PRE
			initialMomentum = [preStepPoint.GetMomentum().x, preStepPoint.GetMomentum().y, preStepPoint.GetMomentum().z]
			# momenta - POST
			# print KE, "\n", p, "\n", initialMomentum, "\n", finalMomentum, "\n\n" 
			# energy = step.GetTotalEnergyDeposit()
			# print p 
			DA.dataCollection(p, m, t) # calls data collection and analysis on final positions and momenta
			# return initialMomentum, finalMomentum 


############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
