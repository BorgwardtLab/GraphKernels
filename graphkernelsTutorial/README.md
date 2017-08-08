# Tutorial: how to compute graph kernels

The file demo_mutag.py, if run when you are in the main directory of the graphkernels package, will compute all the available kernels on the benchmark mutag dataset. 

The data that the graph kernels funcitons require as input, should be a list of igraph objects, as provided in the file data.mutag. This example file is available here in the tutorial, but will also be downloaded together with the package. 

# Kernel computation in Python

1) Import the packages

import numpy as np

import graphkernels.kernels as gk

2) Load the data

mutag_list = np.load('data.mutag')

3) Compute the kernels

K_edge = gk.CalculateEdgeHistKernel(mutag_list)


The matrix K_edge is the kernel matrix, obtained with the edge histogram kernel, and therefore is a square matrix of size equal to the number of samples.

# List of Graph Kernels

The following graph kernels can be computed with our package. 









