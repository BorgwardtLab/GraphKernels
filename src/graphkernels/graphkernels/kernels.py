"""
Collect all the functions
for computing the graph kernels
"""


###################
## Import packages 
###################

import numpy as np
import GKextCPy as gkCpy
from igraph import Graph

from .utilities import GetGKInput, GetAdjMatList




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
    
    if isinstance(par, (int, float, complex)):
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


def CalculateGraphletKernel(G, par = 4):

    # If k<3 then assign k=3
    if par < 3:

        par = 3
        print("Warning: k=3 is used (k = 3 or 4 is supported)")
    
    # If k>4 then assign k=4
    if par > 4:
        
        par = 4
        print("Warning: k=4 is used (k = 3 or 4 is supported)")

    # Extract graph info
    adj_mat, adj_list = GetAdjMatList(G)
    par = int(par)

    K = gkCpy.CalculateGraphletKernelPy(adj_mat, adj_list, par)

    return K



def CalculateConnectedGraphletKernel(G, par = 4):

    # If k<3 then assign k=3
    if par < 3:

        par = 3
        print("Warning: k=3 is used (k = 3, 4 or 5 is supported)")
    
    # If k>5 then assign k=5
    if par > 5:
        
        par = 5
        print("Warning: k=5 is used (k = 3, 4 or 5 is supported)")

    # Extract graph info
    adj_mat, adj_list = GetAdjMatList(G)

    K = gkCpy.CalculateConnectedGraphletKernelPy(adj_mat, adj_list, par)

    return K



def CalculateShortestPathKernel(G):
    
    G_floyd = []
    for i in range(len(G)):

        g_floyd_am = G[i].shortest_paths_dijkstra()
        g_floyd_am = np.asarray(g_floyd_am).reshape(len(g_floyd_am), len(g_floyd_am))
        g = Graph.Adjacency((g_floyd_am > 0).tolist())
        g.es['label'] = g_floyd_am[g_floyd_am.nonzero()]
        g.vs['id']= np.arange(len(G[i].vs['label']))
        g.vs['label'] = G[i].vs['label']
        G_floyd.append(g)


    G_floyd = np.array(G_floyd)

    K = CalculateKStepRandomWalkKernel(G_floyd, par = (0,1))
    return K

