"""
In this demo we show how to
calculate the kernel matrices
on the MUTAG data
"""


from gk_functions import *
import IPython as ip
import numpy as np

# Load data
mutag_list = np.load("/home/eghisu/projects/gk_python_wrapper/gk_python_c/data/mutag_pydata.npy")

### ALL KERNELS COMPUTE

K1 = CalculateEdgeHistKernel(mutag_list)
K2 = CalculateVertexHistKernel(mutag_list)
K3 = CalculateVertexEdgeHistKernel(mutag_list)
K4 = CalculateVertexVertexEdgeHistKernel(mutag_list)
K5 = CalculateEdgeHistGaussKernel(mutag_list)
K6 = CalculateVertexHistGaussKernel(mutag_list)
K7 = CalculateVertexEdgeHistGaussKernel(mutag_list)

#K = CalculateVertexHistKernel(mutag_list)
#par = 5
#K1 = CalculateWLKernel(mutag_list, par)

ip.embed()