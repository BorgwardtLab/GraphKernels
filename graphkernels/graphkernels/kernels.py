"""
Collect all the functions
for computing the graph kernels
"""


###################
## Import packages 
###################

import numpy as np
import GKextCPy as gkCpy

from utilities import GetGKInput

import IPython as ip



#########################
## Alla kernels calculate
#########################


### Edge Histogram Kernel
def CalculateEdgeHistKernel(G, par = -1.0):

	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)

	par = gkCpy.DoubleVector([par])

	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 1)
	return K


###########################
### Vertex Histogram Kernel
def CalculateVertexHistKernel(G, par = -1.0):
	
	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)

	par = gkCpy.DoubleVector([par])
		
	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 2)
	
	return K


################################
### Vertex Edge Histogram Kernel
def CalculateVertexEdgeHistKernel(G, par = -1.0):

	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)
	par = gkCpy.DoubleVector([par])

	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 3)
	return K



#######################################
### Vertex Vertex Edge Histogram Kernel
def CalculateVertexVertexEdgeHistKernel(G, par = 1):
	
	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)

	par = gkCpy.DoubleVector([par])
		
	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 4)
	
	return K


def CalculateEdgeHistGaussKernel(G, par = 1):
	
	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)

	par = gkCpy.DoubleVector([par])
		
	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 5)

	return K



def CalculateVertexHistGaussKernel(G, par = 1):
	
	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)

	par = gkCpy.DoubleVector([par])
		
	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 6)
	
	return K


def CalculateVertexEdgeHistGaussKernel(G, par = 1):
	
	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)

	par = gkCpy.DoubleVector([par])
		
	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 7)
	
	return K



def CalculateGeometricRandomWalkKernel(G, par = 1):
	
	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)

	par = gkCpy.DoubleVector([par])
		
	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 8)
	
	return K



def CalculateExponentialRandomWalkKernel(G, par=1):
	
	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)

	par = gkCpy.DoubleVector([par])
		
	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 9)
	
	return K

def CalculateKStepRandomWalkKernel(G, par=1):
	
	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)
	
	if isinstance(par, (int, long, float, complex)):
		par = gkCpy.DoubleVector([par])

	else:
		par = gkCpy.DoubleVector(par)
		
	K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 10)
	
	return K



def CalculateWLKernel(G, par = 5):

	# Extract graph info
	E, V_label, V_count, E_count, D_max = GetGKInput(G)
	
	#par = nuber of WL iterations
	par = int(par)	
	
	K = gkCpy.WLKernelMatrix(E, V_label, V_count, E_count, D_max, par)
	
	return K


