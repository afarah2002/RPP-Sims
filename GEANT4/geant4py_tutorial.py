from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import random
import time



# def ConstructGeom():
# 	print "* Constructing geometry..."
# 	g4py.ezgeom.Construct()
# 	g4py.NISTmaterials.Construct()
# 	# reset world material
# 	air= G4Material.GetMaterial("G4_AIR")
# 	g4py.ezgeom.SetWorldMaterial(air)

# 	# a target box is placed
# 	global target
# 	target= G4EzVolume("Target")
# 	au= G4Material.GetMaterial("G4_Au")
# 	target.CreateTubeVolume(au, 0., 1.*cm, 1.*mm)
# 	target.PlaceIt(G4ThreeVector(0.,0.,-10.*cm))

# ConstructGeom()

# from EZsim.EZgeom import G4EzVolume
g4py.NISTmaterials.Construct()
# set DetectorConstruction to the RunManager
g4py.ezgeom.Construct()
exN03PL = g4py.EMSTDpl.PhysicsListEMstd()
gRunManager.SetUserInitialization(exN03PL)
# reset world material
air= gNistManager.FindOrBuildMaterial("G4_AIR")
g4py.ezgeom.SetWorldMaterial(air)
# dummy box
detector_box = G4EzVolume("DetectorBox") #i nitialize volume
detector_box.CreateBoxVolume(air, 20.*cm, 20.*cm, 40.*cm)# create volume
detector_boxPos = G4ThreeVector(0.,0.,20.*cm) # specify xyz location of volume
detector_box.PlaceIt(detector_boxPos) # add volume to space
# calorimeter placed inside the box
cal= G4EzVolume("Calorimeter") #initialize volume
nai= gNistManager.FindOrBuildMaterial("G4_SODIUM_IODIDE")
au= G4Material.GetMaterial("G4_Au")

orb = G4EzVolume("Orb") # initialize volume
orbPos = G4ThreeVector(1., 1., 1.) # specify xyz location of volume
orb.CreateOrbVolume(au, 10.*cm)# create volume
orb.PlaceIt(orbPos) # add volume to space

tube = G4EzVolume("Tube") # initialize volume
tubePos = G4ThreeVector(-1., -1., -1.) # specify xyz location of volume
tube.CreateTubeVolume(au, 30., 10.*cm, 30.*cm) # create volume
tube.PlaceIt(tubePos) # add volume to space

gRunManager.Initialize()

angle = 0
while True:
	angle += 5.
	# angle = float(random.randint(1, 101))
	
	gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
	gApplyUICommand("/vis/viewer/create OGLSX oglsxviewer")
	gApplyUICommand("/vis/drawVolume")
	time.sleep(.05)
	# gApplyUICommand("/vis/viewer/zoom " + str(zoom))
	gApplyUICommand("/vis/viewer/set/viewpointThetaPhi 0" + str(angle))
	gApplyUICommand("/vis/viewer/update")

# class GeomBuild(object):

# 	def __init__(self):

# 	def Box(self, material, length, width, height):

# 	def Tube(self, outer_radius, inner_radius, length)
