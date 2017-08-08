# Tutorial: how to compute graph kernels

The file demo_mutag.py can be executed as 

python demo_mutag.py

assuming 

will compute all the available kernels on the benchmark mutag dataset. 

The data that the graph kernels funcitons require as input, should be a list of igraph objects, as provided in the file data.mutag. This example file is available here in the tutorial, but will also be downloaded together with the package. 

# Kernel computation in Python

1) Import the packages

import numpy as np

import graphkernels.kernels as gk

2) Load the data

mutag_list = np.load('data.mutag')

3) Compute the kernels

K_edge = gk.CalculateWLKernel(mutag_list, par = 3)


The matrix K_edge is the kernel matrix, obtained with the edge histogram kernel, and therefore is a square matrix of size equal to the number of samples.

Note that if no par is provided, a default value is used. 

# List of Graph Kernels

The following graph kernels can be computed with our package. 


| Graph Kernel      | Function           | Par  |
| ------------- |:-------------:| -----:|
| Linear kernel between edge histograms	| CalculateEdgeHistKernel |	None |
| Linear kernel between vertex histograms | CalculateVertexHistKernel|	None |
| Linear kernel between vertex-edge histograms | CalculateVertexEdgeHistKernel |	None |
| Linear kernel combination (V + λVE)	| CalculateVertexVertexEdgeHistKernel |	λ |
| Gaussian RBF kernel between vertex histograms	| CalculateVertexHistGaussKernel |	σ |
| Gaussian RBF kernel between edge histograms | CalculateEdgeHistGaussKernel |	σ |
| Gaussian RBF kernel between vertex-edge histograms | CalculateVertexEdgeHistGaussKernel |	σ |
| Geometric random walk kernel | CalculateGeometricRandomWalkKernel |	λ |
| Exponential random walk kernel | CalculateExponentialRandomWalkKernel	| β |
| k-step random walk kernel | CalculateKStepRandomWalkKernel |	λ0, λ1, ..., λk |
| Weisfeiler-Lehman subtree kernel | CalculateWLKernel | h | 





