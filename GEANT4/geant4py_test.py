#----------imports----------#
from Geant4 import *
import g4py.Qmaterials, g4py.NISTmaterials
import g4py.ExN03geom
#import g4py.ExN03pl
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam
import sys
from time import *
from subprocess import *
import os

#----------code starts here!----------#

class BeamInitializer(object):
	'''Description: Creates particle guns for specified particles, energies, positions, and directions
	'''
	def __init__(self, gun_name, particle, energy, E_unit, pos_vect, pos_unit, dir_vect):
		#open datafile


		self.gun_name = gun_name
		rand_engine = Ranlux64Engine()
		HepRandom.setTheEngine(rand_engine)
		HepRandom.setTheSeed(20050830L)

		self.exN03geom = g4py.ExN03geom.ExN03DetectorConstruction()
		gRunManager.SetUserInitialization(self.exN03geom)

		self.exN03PL = g4py.EMSTDpl.PhysicsListEMstd()
		gRunManager.SetUserInitialization(self.exN03PL)

		self.gun_name = g4py.ParticleGun.Construct()

		# set parameters of particle gun
		self.gun_name.SetParticleByName(particle)
		self.gun_name.SetParticleEnergy(energy*E_unit)
		self.gun_name.SetParticlePosition(G4ThreeVector(pos_vect[0], pos_vect[1], pos_vect[2])*pos_unit)
		self.gun_name.SetParticleMomentumDirection(G4ThreeVector(dir_vect[0], dir_vect[1] , dir_vect[2]))

		gRunManager.Initialize()

	def cmd_beamOn(self):

		#set moderator parameters
		#set number of mod layers
		self.layerVar = 10
		self.exN03geom.SetNbOfLayers(self.layerVar)
		#set mod material
		self.absorbermaterialVar = "Lead"
		self.exN03geom.SetAbsorberMaterial(self.absorbermaterialVar)
		#set mod layer thickness
		self.absorberthickVar = 10.0
		self.exN03geom.SetAbsorberThickness(self.absorberthickVar  * mm/2.0)
		#set substance between mod layers
		self.gapmaterialVar = "liquidArgon"
 		self.exN03geom.SetGapMaterial(self.gapmaterialVar)
		#set gap size between layers
		self.gapthickVar = 5.0
		self.exN03geom.SetGapThickness(self.gapthickVar  * mm/2.0)

		# self.exN03geom.SetCalorSizeYZ(self.calorsizeYZVar.get() * mm)
		position = -self.layerVar*(self.absorberthickVar + self.gapthickVar)*1.2

		self.exN03geom.UpdateGeometry()
		# self.exN03PL.SetDefaultCutValue(self.cutVar.get() * mm)
		# self.exN03PL.SetCutsWithDefault()
		self.magVar = 0.
		self.exN03geom.SetMagField(self.magVar * tesla)

		# print("Now geometry updated")


		# self.cmd_particle(self.particleVar.get())
		# self.cmd_energy(self.energyVar.get())

		# print position

		eventNum = 10
		for i in range(eventNum):

			self.gun_name.SetParticlePosition(G4ThreeVector(position*mm, (i-eventNum/2)*5.*mm, 0.*cm))
			gRunManager.BeamOn(1)
			# sleep(0.1)
		# gApplyUICommand("/vis/viewer/update")


	def writeOutput(self):
		self.f.close()
		#close datafile



if __name__ == '__main__':
	gun_name = "pg1"
	f = open("DATAFILES/" + gun_name + "data.txt", 'w+')
	with sys.stdout as f:
		BI = BeamInitializer(gun_name, "e+", 50., MeV, (-40.,0.,0.), cm, (1.,0.,0.))
		BI
		BI.cmd_beamOn()
	BI.writeOutput()

