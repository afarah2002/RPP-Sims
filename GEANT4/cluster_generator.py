#----------PYTHON IMPORTS-------------#
import numpy as np
#----------code starts here!----------#
class ClusterGenerator(object):

	def __init__(self):
		pass
	def sphericalClusters(self, steps):

		spherical_coor_LIST = []
		pi = np.pi
		step = 6
		for phi in np.arange(0, pi, pi/step): # smaller steps means more clusters, range goes to pi since clusters are double sided
			for theta in np.arange(0, pi, pi/step):
				sph_coor = [theta, phi] # phi, theta
				spherical_coor_LIST.append(sph_coor)

		return spherical_coor_LIST
	def cartestianClusters(self, cluster_width, width, spacing, edge):
		cluster_coor_LIST = []

		for x in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval, interval):
			for y in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval, interval):
				for z in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval, interval):
					cluster_coor = [x,y,z]
					# cluster_coor_LIST.append(cluster_coor)
					# print cluster_coor
					# time.sleep(1)
					opposite = np.multiply(cluster_coor, -1)
					for i in cluster_coor:
						if np.abs(i) == edge-cluster_width/2 - spacing and cluster_coor not in cluster_coor_LIST:
							cluster_coor_LIST.append(cluster_coor)
							break

		return cluster_coor_LIST