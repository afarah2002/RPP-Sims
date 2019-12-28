#----------imports----------#
from Geant4 import * 
import random
import numpy as np
import scipy.stats as ss
from scipy.signal import find_peaks
import pandas as pd
import seaborn as sns  # for nicer graphics

## matplotlib stuff
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from time import sleep
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import time

global vectorCount
vectorCount = 500 # number of scattered e+ per run

global SEE_count
SEE_count = 0


#----------code starts here!----------#
class WipeData(object):
	# wipe lists for next data collection
	def wipe(self):
		pass

WIPE = WipeData()

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class Plotter(object):
	"graphs 3D positions"

	def __init__(self):
		pass

	def dataCollection(self, posf, momf, tcluster):

		pass


	def dataAnalysis(self):

		pass



	def dataReturner(self):

		pass

	def grapher(self):

		pass


		# plt.show()

PLT = Plotter()


class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self, energy, p_LIST, m_LIST):
		G4VUserPrimaryGeneratorAction.__init__(self)
		self.particleGun = G4ParticleGun(1)
		# print("\n Particle gun defined \n")
		self.energy = energy
		self.positions_LIST = p_LIST
		self.momenta_LIST = m_LIST

	def GeneratePrimaries(self, event):


		# Particle param
		#################################################
		particle = "e+"
		# energy_2 = 2.5
		energyUnit = MeV 
		dimensionUnit = cm

		energy = self.energy
		self.particleGun.SetParticleByName(particle) # define particle
		self.particleGun.SetParticleEnergy(energy*energyUnit) # define particle energy 

		
		for position in self.positions_LIST: # creates random momentum vectors originating from [0, 0, 0]
			momentumArray = self.momenta_LIST[self.positions_LIST.index(position)]
			self.particleGun.SetParticlePosition(G4ThreeVector(position[0], position[1], position[2])*dimensionUnit) # define first particle generator location
			self.particleGun.SetParticleMomentumDirection(G4ThreeVector(momentumArray[0]/1000, momentumArray[1]/1000, momentumArray[2]/1000)*dimensionUnit) # define first particle generator momentum
			self.particleGun.GeneratePrimaryVertex(event)
		#################################################

#-------------------------------------------------------------------
class MyRunAction(G4UserRunAction):
	"My Run Action"

	def EndOfRunAction(self, run):
		global SEE_count
		SEE_count = 0
		# time.sleep(1)
		# pass

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

		# change = step.particleChange()
		track = step.GetTrack()
		touchable = track.GetTouchable()
		KE = track.GetKineticEnergy()
		parentId = track.GetParentID()
		particleName = track.GetDefinition().GetParticleName() 

		# kinetic energy in MeV - PRE
		# initialKE = preStepPoint.GetKineticEnergy() 
		# kinetic energy in MeV - POST
		# finalKE = postStepPoint.GetKineticEnergy()


		p_test = [step.GetDeltaPosition().x,step.GetDeltaPosition().y,step.GetDeltaPosition().z]
		p = [postStepPoint.GetPosition().x, postStepPoint.GetPosition().y, postStepPoint.GetPosition().z] # (mm)
		# p and p_test are the SAME 
		t_test = step.GetDeltaTime()
		t = track.GetGlobalTime() # (ns)
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
		# PLT.dataCollection(p, m, t) # calls data collection and analysis on final positions and momenta
		# return initialMomentum, finalMomentum 


		if particleName == 'e-' and parentId != 0:
			global SEE_count
			SEE_count += 1
			print SEE_count
