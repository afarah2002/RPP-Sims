#----------imports----------#
from Geant4 import * 
import random

#----------code starts here!----------#
class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self):
		G4VUserPrimaryGeneratorAction.__init__(self)
		self.particleGun = G4ParticleGun(1)
		print("\n Particle gun defined \n")

	def GeneratePrimaries(self, event):


		# Particle param
		#################################################
		spaceParamDict = {}
		vectorCount = 100
		locationArray = [0, 0, 0]

		particle = "e+"

		energy_1 = 2.5	
		energy_2 = 2.5

		energyUnit = MeV
		dimensionUnit = cm

		self.particleGun.SetParticleByName(particle) # define particle
		self.particleGun.SetParticleEnergy(energy_1*energyUnit) # define particle energy 

		for i in range(0,vectorCount): # creates random momentum vectors originating from [0, 0, 0]
			px = random.uniform(-1,1)
			py = random.uniform(-1,1)
			pz = random.uniform(-1,1)
			momentumArray = [px, py, pz]
			self.particleGun.SetParticlePosition(G4ThreeVector(locationArray[0], locationArray[1], locationArray[2])*dimensionUnit) # define first particle generator location
			self.particleGun.SetParticleMomentumDirection(G4ThreeVector(momentumArray[0], momentumArray[1], momentumArray[2])*dimensionUnit) # define first particle generator momentum
			self.particleGun.GeneratePrimaryVertex(event)
		#################################################

#-------------------------------------------------------------------
class MyRunAction(G4UserRunAction):
	"My Run Action"

	def EndOfRunAction(self, run):
		print "*** End of Run"
		print "- Run sammary : (id= %d, #events= %d)" \
		% (run.GetRunID(), run.GetNumberOfEventToBeProcessed())

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
	"My Event Action"

	def EndOfEventAction(self, event):
		pass

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
	"My Stepping Action"

	def UserSteppingAction(self, step):
		#print "*** dE/dx in current step=", step.GetTotalEnergyDeposit()
		preStepPoint = step.GetPreStepPoint()
		postStepPoint = step.GetPostStepPoint()
		energy = step.GetTotalEnergyDeposit()
		print(
			postStepPoint.GetCharge(), 
			postStepPoint.GetKineticEnergy(), # MeV
			[preStepPoint.GetMomentum().x, postStepPoint.GetMomentum().x] 
			# postStepPoint.GetMomentum().y, postStepPoint.GetMomentum().z]
			)
		track = step.GetTrack()
		touchable = track.GetTouchable()
		#print " *** vid= ", touchable.GetReplicaNumber()
		pass

class MyField(G4MagneticField):
	"My Magnetic Field"

	def GetFieldValue(self, pos, time):

		vectorList = [
						# [0, 1, 0], 
					 	[-10., -10., -10.]
					 ]
		for v in vectorList:

			bfield = G4ThreeVector()
			bfield.x = v[0]*tesla
			bfield.y = v[1]*tesla
			bfield.z = v[2]*tesla
			return bfield

# class ScoreSD(G4VSensitiveDetector):
# 	"SD for score voxels"

# 	def __init__(self):
# 		G4VSensitiveDetector.__init__(self, "ScoreVoxel")

# 	def ProcessHits(self, step, rohist):
# 		preStepPoint = step.GetPreStepPoint()
# 		if(preStepPoint.GetCharge() == 0):
# 			return

# 		track = step.GetTrack()
# 		touchable = track.GetTouchable()
# 		voxel_id = touchable.GetReplicaNumber()
# 		dedx = step.GetTotalEnergyDeposit()
# 		xz = posXZ(voxel_id)
# 		# hist_dosezse2d.Fill(xz[1], xz[0], dedx/MeV)
# 		if( abs(xz[0]) <= 100 ):
# 			# hist_dosez.Fill(xz[1],  dedx/MeV)
# 			pass



		# # Magnetic field param
		# #################################################
		# magUnit = tesla
		# magVectorList = [
		# 				[.5,.5,.5],
		# 				# [0,0,0]
		# 				]

		# for vector in magVectorList:
		# 	fieldMgr = gTransportationManager.GetFieldManager()
		# 	vOrigin = G4ThreeVector(0., 0., 0.)
		# 	myField = G4UniformMagField(G4ThreeVector(vector[0]*magUnit, vector[1]*magUnit, vector[2]*magUnit))

		# return myField
		# #################################################
