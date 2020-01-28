#----------GEANT4 imports----------#
from Geant4 import *
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.NISTmaterials
import g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam
#----------PYTHON imports----------#
import collections
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import random
import thread
import time
import scipy.stats as ss
from scipy import optimize, stats
import sys
#-----------FILE imports------------#
from geom_constructor import GeomConstructor 
# -----------CLASS ASSIGNMENTS----------- #
GC = GeomConstructor()
# --------------------------------------- #

#----------code starts here!----------#
uniqueClusters = [[None, None, None]]

class FieldDesign(object):

	def cartesianfieldParam(self, energy, b, x, y, z, cluster_width, edge):
		
		radius = np.sqrt(x**2 + y**2 + z**2)

		mag0 = x*b/radius
		mag1 = y*b/radius
		mag2 = z*b/radius

		# mag0 = np.sqrt(3*b**2)
		# mag1 = np.sqrt(3*b**2)
		# mag2 = np.sqrt(3*b**2)

		magVec = [mag0, mag1, mag2] 
		# print magVec
		
		# get index of "max" value in magVec, most significant axis
		maxIndex = list(np.abs(magVec)).index(max(np.abs(magVec))) 

		#get scaling factor
		scale = edge / np.abs(magVec[maxIndex])
		magVecScaled = np.multiply(magVec, scale)

		receiverDimScaled = [cluster_width, cluster_width, cluster_width]
		receiverDimScaled[maxIndex] = 1

		# print "Receiver DIMENSIONS  ", receiverDimScaled

		material1 = G4Material.GetMaterial("G4_W")
		# GC.ConstructBox("Receiver", material1, [x,y,z], mm, receiverDimScaled)
		# GC.ConstructBox("Receiver", material1, np.multiply([x,y,z], -1), mm, receiverDimScaled) # receiver for opposite cluster

		# if the magVecScaled is not in the uniqueClusters list, append to it
		flag = 0
		for elem in uniqueClusters:
			if collections.Counter(elem) == collections.Counter(magVecScaled):
				# print True
				flag = 1

		if flag == 0:
			uniqueClusters.append(magVecScaled)
		else: 
			pass


		return magVec, magVecScaled

	def spherefieldParam(self, energy, b, phi, theta, cluster_width, edge):

		# vectorList = [list(np.multiply([energy, energy, energy], be))]
		# radius = np.sqrt(np.square(vectorList[0][0]) + np.square(vectorList[0][1]) + np.square(vectorList[0][2]))

		mag0 = b*(np.sin(phi)*np.cos(theta))
		mag1 = b*(np.sin(phi)*np.sin(theta))
		mag2 = b*np.cos(phi)

		magVec = [mag0, mag1, mag2] 
		
		# get index of "max" value in magVec, most significant axis
		maxIndex = list(np.abs(magVec)).index(max(np.abs(magVec))) 

		#get scaling factor
		scale = edge / np.abs(magVec[maxIndex])
		magVecScaled = np.multiply(magVec, scale)

		receiverDimScaled = [cluster_width, cluster_width, cluster_width]
		receiverDimScaled[maxIndex] = 1

		# print "Receiver DIMENSIONS  ", receiverDimScaled

		material1 = G4Material.GetMaterial("G4_W")
		# GC.ConstructBox("Receiver", material1, magVecScaled, mm, receiverDimScaled)
		# GC.ConstructBox("Receiver", material1, np.multiply(magVecScaled, -1), mm, receiverDimScaled) # receiver for opposite cluster

		# if the magVecScaled is not in the uniqueClusters list, append to it
		flag = 0
		for elem in uniqueClusters:
			if collections.Counter(elem) == collections.Counter(magVecScaled):
				flag = 1

		if flag == 0:
			uniqueClusters.append(magVecScaled)
		else: 
			pass


		return magVec, magVecScaled