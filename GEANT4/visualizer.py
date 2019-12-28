from Geant4 import *
class Visualizer(object):

	def visualizer(self, viz_theta, viz_phi):
		# time.sleep(.1)

		gApplyUICommand("/vis/sceneHandler/create OGLSX OGLSX")
		gApplyUICommand("/vis/viewer/create OGLSX oglsxviewer")
		gApplyUICommand("/vis/drawVolume")
		# gApplyUICommand("/globalField/verbose level")
		# gApplyUICommand("/MagneticFieldModels/UniformField/SetBVec 3. 1. 1. Tesla")

		gApplyUICommand("/vis/viewer/select oglsxviewer")
		gApplyUICommand("/vis/ogl/set/displayListLimit 100000")
		gApplyUICommand("/vis/scene/add/trajectories")

		gApplyUICommand("/tracking/storeTrajectory 1")
		gApplyUICommand("/vis/scene/endOfEventAction accumulate")
		gApplyUICommand("/vis/scene/endOfRunAction accumulate")

		gApplyUICommand("/vis/viewer/set/viewpointThetaPhi " + str(viz_theta) + " " + str(viz_phi))
		# gApplyUICommand("/vis/viewer/zoom 1.00001")
		
		gApplyUICommand("/vis/viewer/update")