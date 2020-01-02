#----------PYTHON IMPORTS-------------#
import numpy as np
#----------code starts here!----------#
class ClusterGenerator(object):

	def __init__(self):
		pass
	def sphericalClusters(self, step):

		spherical_coor_LIST = []
		pi = np.pi
		for phi in np.arange(0, pi, pi/step): # smaller steps means more clusters, range goes to pi since clusters are double sided
			for theta in np.arange(0, pi, pi/step):
				sph_coor = [theta, phi] # phi, theta
				spherical_coor_LIST.append(sph_coor)

		return spherical_coor_LIST
	def cartesianClusters(self, cluster_width, edge):
		cluster_coor_LIST = []

		leftover_space = edge*2 % cluster_width
		n_clus = int(2*edge/cluster_width)

		spacing = leftover_space/n_clus

		interval = cluster_width + spacing
		

		for x in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval+1, interval):
			for y in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval+1, interval):
				for z in np.arange(-(edge-cluster_width/2 - spacing), edge - (cluster_width/2 + spacing) + interval+1, interval):
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