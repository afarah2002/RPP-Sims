'''
secondary electron emission analysis script - 
unfortunately, GEANT4 does not do SE generation :(
'''
#----------imports----------#
from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam

import scipy.stats as ss
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages
import random
import time
import thread
import numpy as np
from scipy import optimize, stats
#----file imports--------#
from geom_constructor import GeomConstructor 
# from beam import BeamInitializer
from beam3 import ClusteredPositronGenerator, MyRunAction, MyEventAction, MySteppingAction
from visualizer import Visualizer

#----------code starts here!----------#
energy_LIST = [.002] 

class SEEConstructor(object):
	def __init__(self):
		# from EZsim.EZgeom import G4EzVolume
		g4py.NISTmaterials.Construct()
		# set DetectorConstruction to the RunManager
		# G4VUserDetectorConstruction.Construct()
		g4py.ezgeom.Construct()
		exN03PL = g4py.EMSTDpl.PhysicsListEMstd()
		gRunManager.SetUserInitialization(exN03PL)
		# reset world material
		air = gNistManager.FindOrBuildMaterial("G4_AIR")
		# g4py.ezgeom.SetWorldMaterial(air)




	def construct(self):
		# materials
		NaI = gNistManager.FindOrBuildMaterial("G4_SODIUM_IODIDE")
		GRAPHITE = gNistManager.FindOrBuildMaterial("G4_GRAPHITE")
		Al = gNistManager.FindOrBuildMaterial("G4_Al")
		W = G4Material.GetMaterial("G4_W")

		scale_factor = 5

		dist_cent = 60 # from 1 to 100
		thickness = 30e-9
		width = 1000

		# GC.ConstructSphere("Sphere", GRAPHITE , [0., 0., 0.], m, (.10-30e-9)*scale_factor, .10*scale_factor, 0., 360., 0., 90)
		# GC.ConstructBox("Box", W, [(dist_cent)*scale_factor, 0., 0.], mm, [ thickness,width,width]) # ON X-AXIS
		# GC.ConstructBox("Box", W, [-(dist_cent)*scale_factor, 0., 0.], mm, [ thickness,width,width]) # ON X-AXIS

		# GC.ConstructBox("Box", Al, [0., (dist_cent)*scale_factor, 0.], mm, [width, thickness, width]) # ON Y-AXIS
		# GC.ConstructBox("Box", Al, [0., -(dist_cent)*scale_factor, 0.], mm, [width, thickness, width]) # ON Y-AXIS

		# GC.ConstructBox("Box", Al, [0., 0., (dist_cent)*scale_factor], mm, [width,width, thickness]) # ON Z-AXIS
		# GC.ConstructBox("Box", Al, [0., 0., -(dist_cent)*scale_factor], mm, [width,width, thickness]) # ON Z-AXIS


GC = GeomConstructor()
SEEConstructor = SEEConstructor()
VIS = Visualizer()
# # CC = ClusterClass()
viz_theta = 0
viz_phi = 0


class SecondaryElectronEmissionProcess(object):


	def runSEE(self, energy, positions, momenta, location):

		print "energy = ", energy, " eV"
		SEEConstructor
		SEEConstructor.construct()
		gRunManager.GeometryHasBeenModified()
		PGA_1 = ClusteredPositronGenerator(energy, positions, momenta)
		gRunManager.SetUserAction(PGA_1)

		myEA = MyEventAction()
		gRunManager.SetUserAction(myEA)

		fieldMgr = gTransportationManager.GetFieldManager()
		myField = G4UniformMagField(G4ThreeVector(0, 0, 0))

		# myField = MyField(1)
		fieldMgr.SetDetectorField(myField)
		fieldMgr.CreateChordFinder(myField)

		mySA = MySteppingAction()
		gRunManager.SetUserAction(mySA)

		myRA = MyRunAction()
		gRunManager.SetUserAction(myRA)

		gRunManager.Initialize()

		gRunManager.BeamOn(1)

		VIS.visualizer(viz_theta, viz_phi, "cluster_gen")
		
		# gApplyUICommand("/vis/viewer/update")

		time.sleep(2)


