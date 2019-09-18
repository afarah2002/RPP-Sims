#----------imports----------#
from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import random
import time

#----file imports--------#
from geom_constructor import GeomConstructor 

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
GC = GeomConstructor()
GC.ConstructBox("Detector Box", air, [0., 0., 0.], cm, [20., 20., 40.])

# calorimeter placed inside the box
cal= G4EzVolume("Calorimeter") #initialize volume
nai= gNistManager.FindOrBuildMaterial("G4_SODIUM_IODIDE")
au= G4Material.GetMaterial("G4_Au")

GC.ConstructOrb("Orb", nai, [10., 10., 10.], cm, 10.)

GC.ConstructTube("Tube", au, [-1., -1., -1.], cm, 10., 30., 30., 0, 300.)

GC.ConstructSphere("Sphere", au, [0., 0., 0.], cm, 0., 10., 0., 300., 0., 150)

GC.ConstructCone("Cone", nai, [-20., -20., -20.], cm, 0., 20., 0., 0., 25., 0., 180) # dphi = 359.9999 is basically 360, but we can still see it
# sphere = G4EzVolume("Sphere")
# spherePos = G4ThreeVector(1., 1., 1.)
# sphere.CreateShpereVolume(au, 1.*cm, .5*m) # haha "Shpere" not "Sphere"
# sphere.PlaceIt(spherePos)

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
