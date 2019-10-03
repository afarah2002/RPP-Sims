from Geant4 import *
import g4py.Qmaterials, g4py.NISTmaterials
import g4py.ExN03geom
from g4py.ezgeom import G4EzVolume
#import g4py.ExN03pl
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam
import sys
from time import *
from subprocess import *
import os

class BeamInitializer(G4VUserPrimaryGeneratorAction):

	def __init__(self):
		G4VUserPrimaryGeneratorAction.__init__(self)
		# pgPGA= g4py.ParticleGun.ParticleGunAction()
		# gRunManager.SetUserAction(pgPGA)
		# self.beam = pgPGA.GetParticleGun()
		self.beam = G4ParticleGun(1)

	def parameters(self, particle, energy, energy_unit, dimension_unit, locationArray, directionArray):
		self.particle = particle
		self.energy = energy
		self.energy_unit = energy_unit
		self.dimension_unit = dimension_unit
		self.locationArray = locationArray
		self.directionArray = directionArray
		# self.beam = G4ParticleGun(2)
		# self.beam = g4py.ParticleGun.Construct()

	def GeneratePrimaries(self, event): # set beam parameters
		# beam.GeneratePrimaryVertex(beam)
		self.beam.SetParticleByName(self.particle)
		self.beam.SetParticleEnergy(self.energy*self.energy_unit)
		self.beam.SetParticlePosition(G4ThreeVector(self.locationArray[0], self.locationArray[1], self.locationArray[2])*self.dimension_unit)
		self.beam.SetParticleMomentumDirection(G4ThreeVector(self.directionArray[0], self.directionArray[1], self.directionArray[2])*self.dimension_unit)
		# print(type(self.beam))

		self.beam.GeneratePrimaryVertex(event) # runs beam

	# def generatePrimaries(self, event):
	# 	# sleep(1)
	# 	# gRunManager.SetUserInitialization(self.beam)
	# 	# self.beam.G4PrimaryVertex(event)
	# 	pass



