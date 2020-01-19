# RPP Simulations 

## Description: Simulating a new setup for efficient radioisotope positron propulsion

# Contents:
 - GEANT4 Python sims for tracking e+ in scattering, moderation, etc. --> "GEANT4/"
	 - geom_constructor.py : Streamlined volume constructor functions for box, tube, cone, sphere, and orb
	 - visualizer.py - OpenGL viewer class
	 - field_designer.py - cartesian and spherical cluster generator class 
	 - cluster_generator.py - cluster placment based on field_designer.py
	 - beam2_1.py - particle generator, cluster analysis, cluster plotter, event/step moderation
	 - beam3.py - generates individual clusters in space
	 - system_main.py - main cluster simulator 
