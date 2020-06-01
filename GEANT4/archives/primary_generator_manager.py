#----------imports----------#
from Geant4 import *

#----------code starts here!----------#
class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
  "My Primary Generator Action"

  def __init__(self):
    G4VUserPrimaryGeneratorAction.__init__(self)
    self.particleGun= G4ParticleGun(1)

  def GeneratePrimaries(self, event):
    #dx= random.gauss(0., 0.1)
    dx=0.
    self.particleGun.SetParticleMomentumDirection(G4ThreeVector(dx, 0., 1.))
    self.particleGun.GeneratePrimaryVertex(event)

# ------------------------------------------------------------------
class MyRunAction(G4UserRunAction):
  "My Run Action"

  def BeginOfRunAction(self, run):
    print "*** #event to be processed (BRA)=",
    run.GetNumberOfEventToBeProcessed()

  def EndOfRunAction(self, run):
    print "*** run end run(ERA)=", run.GetRunID()

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
  "My Event Action"

  #def BeginOfEventAction(self, event):
    #print "*** current event (BEA)=", event.GetEventID()
