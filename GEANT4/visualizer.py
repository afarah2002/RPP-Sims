'''
This class provides basic parameters to view the 
GEANT4PY environment. Call it with the viewing 
angles and view name to open a viewer. 
'''

from Geant4 import *
class Visualizer(object):

	def visualizer(self, viz_theta, viz_phi, viewer_name):

		gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
		gApplyUICommand("/vis/viewer/create OGLSX " + viewer_name)
		gApplyUICommand("/vis/drawVolume")

		gApplyUICommand("/vis/viewer/select " + viewer_name)
		gApplyUICommand("/vis/ogl/set/displayListLimit 100000")
		gApplyUICommand("/vis/scene/add/trajectories")

		gApplyUICommand("/tracking/storeTrajectory 1")
		gApplyUICommand("/vis/scene/endOfEventAction accumulate")
		gApplyUICommand("/vis/scene/endOfRunAction accumulate")

		gApplyUICommand("/vis/viewer/set/viewpointThetaPhi " + str(viz_theta) + " " + str(viz_phi))
		# gApplyUICommand("/vis/viewer/zoom 1.00001")
		
		gApplyUICommand("/vis/viewer/update")