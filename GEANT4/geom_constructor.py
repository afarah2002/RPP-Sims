'''
This geometry generator simplifies the construction of sumple 
3D shapes in a GEANT4PY environment. 
Includes: Box/Rectangular Prism, Tube, Cone, Sphere, Orb
'''

#----------imports----------#
from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl

#----------code starts here!----------#


class GeomConstructor(object):

	def ConstructBox(self, nameString, material, locationArray, unit, dimensionsArray):
		box = G4EzVolume(nameString) #initialize volume
		boxPos = G4ThreeVector(locationArray[0]*unit, locationArray[1]*unit, locationArray[2]*unit) # specify xyz location of volume
		box.CreateBoxVolume(material, dimensionsArray[0]*unit, dimensionsArray[1]*unit, dimensionsArray[2]*unit)# create volume
		box.PlaceIt(boxPos) # add volume to space

	def ConstructTube(self, nameString, material, locationArray, unit, minRadius, maxRadius, length, phi0, dphi):
		tube = G4EzVolume(nameString) # initialize volume
		tubePos = G4ThreeVector(locationArray[0]*unit, locationArray[1]*unit, locationArray[2]*unit) # specify xyz location of volume
		tube.CreateTubeVolume(material, minRadius*unit, maxRadius*unit, length*unit, phi0, dphi*deg) # create volume
		tube.PlaceIt(tubePos) # add volume to space

	def ConstructCone(self, nameString, material, locationArray, unit, minRadius1, maxRadius1, minRadius2, maxRadius2, height, phi0, dphi):
		cone = G4EzVolume(nameString) #initialize volume
		conePos = G4ThreeVector(locationArray[0]*unit, locationArray[1]*unit, locationArray[2]*unit) # specify xyz location of volume
		cone.CreateConeVolume(material, minRadius1*unit, maxRadius1*unit, minRadius2*unit, maxRadius2*unit, height*unit, phi0, dphi*deg) # create volume
		cone.PlaceIt(conePos) # add volume to space

	def ConstructSphere(self, nameString, material, locationArray, unit, minRadius, maxRadius, phi0, dphi, theta0, dtheta):
		sphere = G4EzVolume(nameString) #initialize volume
		spherePos = G4ThreeVector(locationArray[0]*unit, locationArray[1]*unit, locationArray[2]*unit) # specify xyz location of volume
		sphere.CreateShpereVolume(material, minRadius*unit, maxRadius*unit, phi0, dphi*deg, theta0, dtheta*deg) # haha "Shpere" not "Sphere"  --> check docs https://nusoft.fnal.gov/larsoft/doxsvn/html/pyEzgeom_8cc_source.html#l00116
		sphere.PlaceIt(spherePos) # add volume to space

	def ConstructOrb(self, nameString, material, locationArray, unit, maxRadius):
		orb = G4EzVolume(nameString) # initialize volume
		orbPos = G4ThreeVector(locationArray[0]*unit, locationArray[1]*unit, locationArray[2]*unit) # specify xyz location of volume
		orb.CreateOrbVolume(material, maxRadius*unit)# create volume
		orb.PlaceIt(orbPos) # add volume to space
