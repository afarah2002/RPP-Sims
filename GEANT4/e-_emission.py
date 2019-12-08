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
from beam3 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction, MySteppingAction, Plotter, WipeData
from geant4py_tutorial import Constructor, Visualizer, ClusterClass

#----------code starts here!----------#

Constructor = Constructor
Constructor
Constructor.construct()
VIS = Visualizer()
CC = ClusterClass()

viz_theta = 35
viz_phi = 35
energy = .250 # MeV 

if __name__ == '__main__':
	while True: 
		PGA_1 = MyPrimaryGeneratorAction(energy, [0,0,0], [1,1,1])
		gRunManager.SetUserAction(PGA_1)

		myEA = MyEventAction()
		gRunManager.SetUserAction(myEA)

		mySA = MySteppingAction()
		gRunManager.SetUserAction(mySA)

		myRA = MyRunAction()
		gRunManager.SetUserAction(myRA)

		gRunManager.Initialize()

		gRunManager.BeamOn(1)

		VIS.visualizer(viz_theta, viz_phi)