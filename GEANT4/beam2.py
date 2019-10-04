#----------imports----------#
from Geant4 import * 

#----------code starts here!----------#
class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
	"My Primary Generator Action"

	def __init__(self):
		G4VUserPrimaryGeneratorAction.__init__(self)
		self.particleGun = G4ParticleGun(1)
		print("Particle gun defined \n")

	def GeneratePrimaries(self, event):

		# locationArray_TEST = [0., 0., 0.]
		# momentumArray_TEST = [10., 10., 10.]
		# dimensionUnit = cm

		# ---- particle paramteters ---- #
		particle = "e+"

		energy_1 = 50	
		energy_2 = 50

		energyUnit = MeV
		dimensionUnit = cm

		locationArray_1 = [-9.5*5+.1, 0., 0.]
		momentumArray_1 = [1., 0., 0.]

		locationArray_2 = [9.5*5-.1, 0., 0.]
		momentumArray_2 = [-1., 0., 0.]
		# ------------------------------ #

		self.particleGun.SetParticleByName(particle) # define particle
		self.particleGun.SetParticleEnergy(energy_1*energyUnit) # define particle energy 

		self.particleGun.SetParticlePosition(G4ThreeVector(locationArray_1[0], locationArray_1[1], locationArray_1[2])*dimensionUnit) # define first particle generator location
		self.particleGun.SetParticleMomentumDirection(G4ThreeVector(momentumArray_1[0], momentumArray_1[1], momentumArray_1[2])*dimensionUnit) # define first particle generator momentum
		self.particleGun.GeneratePrimaryVertex(event)

		self.particleGun.SetParticlePosition(G4ThreeVector(locationArray_2[0], locationArray_2[1], locationArray_2[2])*dimensionUnit) # define second particle generator location
		self.particleGun.SetParticleMomentumDirection(G4ThreeVector(momentumArray_2[0], momentumArray_2[1], momentumArray_2[2])*dimensionUnit) # define second particle generator momentum
		self.particleGun.GeneratePrimaryVertex(event)
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
