"""
In this demo we show how to
calculate the kernel matrices
on the MUTAG data
"""

import graphkernels.kernels as gk

import IPython as ip
import numpy as np

# Load data
# comment next line if data.mutag is in a different folder
mutag_list = np.load("graphkernels/data.mutag", allow_pickle=True)

# Uncomment next line if data.mutag is in your current folder
#mutag_list = np.load("data.mutag", allow_pickle=True)

### ALL KERNELS COMPUTE
K1 = gk.CalculateEdgeHistKernel(mutag_list)
K2 = gk.CalculateVertexHistKernel(mutag_list) 
K3 = gk.CalculateVertexEdgeHistKernel(mutag_list)
K4 = gk.CalculateVertexVertexEdgeHistKernel(mutag_list)
K5 = gk.CalculateEdgeHistGaussKernel(mutag_list)
K6 = gk.CalculateVertexHistGaussKernel(mutag_list)
K7 = gk.CalculateVertexEdgeHistGaussKernel(mutag_list)
K8 = gk.CalculateGeometricRandomWalkKernel(mutag_list)
K9 = gk.CalculateExponentialRandomWalkKernel(mutag_list)
K10 = gk.CalculateKStepRandomWalkKernel(mutag_list)
K11 = gk.CalculateWLKernel(mutag_list)
K12 = gk.CalculateConnectedGraphletKernel(mutag_list, 4)
K13 = gk.CalculateGraphletKernel(mutag_list, 4)
K14 = gk.CalculateShortestPathKernel(mutag_list)

