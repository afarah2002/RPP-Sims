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
from beam2 import MyPrimaryGeneratorAction, MyRunAction, MyEventAction, MySteppingAction, MyField, Plotter, WipeData
from geant4py_tutorial import Constructor, Visualizer

#----------code starts here!----------#

Constructor = Constructor
Constructor
Constructor.construct()
VIS = Visualizer()

angle = 35
energy = 2.5 
if __name__ == '__main__':

	PGA_1 = MyPrimaryGeneratorAction(energy)
	gRunManager.SetUserAction(PGA_1)

	myEA = MyEventAction()
	gRunManager.SetUserAction(myEA)

	mySA = MySteppingAction()
	gRunManager.SetUserAction(mySA)

	myRA = MyRunAction()
	gRunManager.SetUserAction(myRA)

	gRunManager.Initialize()

	gRunManager.BeamOn(1)

	VIS.visualizer(angle)