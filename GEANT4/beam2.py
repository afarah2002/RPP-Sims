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
		#print "*** vid= ", touchable.GetReplicaNumber()
