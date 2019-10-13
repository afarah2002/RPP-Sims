#----------imports----------#
from Geant4 import * 
import random

#----------code starts here!----------#
class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self):
		G4VUserPrimaryGeneratorAction.__init__(self)
		self.particleGun = G4ParticleGun(1)
		print("Particle gun defined \n")

	def GeneratePrimaries(self, event):
		# Magnetic field param
		#################################################
		magUnit = tesla
		magVectorList = [[10,10,10],[-1,-1,-1]]
		for vector in magVectorList:
			fieldMgr = gTransportationManager.GetFieldManager()
			myField = G4UniformMagField(G4ThreeVector(vector[0]*magUnit, vector[1]*magUnit, vector[2]*magUnit))
			#myField= MyField()
			fieldMgr.SetDetectorField(myField)
			fieldMgr.CreateChordFinder(myField)
		#################################################

		# Particle param
		#################################################
		spaceParamDict = {}
		vectorCount = 100
		locationArray = [0, 0, 0]

		particle = "e+"

		energy_1 = 50	
		energy_2 = 50

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
		pass
		#print "*** dE/dx in current step=", step.GetTotalEnergyDeposit()
		preStepPoint= step.GetPreStepPoint()
		track= step.GetTrack()
		touchable= track.GetTouchable()
		#print " *** vid= ", touchable.GetReplicaNumber()


