"""
Utility funcitons for the graph kernels package
"""

####################
## Import packages
####################

import numpy as np
import igraph

import GKextCPy as gkCpy



##################################
## Extract graphs information from 
## an igraph object
###################################

def GetGraphInfo(g):
        
    # matrix of edges
    E = np.zeros(shape = (len(g.es),2)) 
    for i in range(len(g.es)):
        E[i,:] = g.es[i].tuple
    ## there are multiple edge attributes
    if len(g.es.attributes()) > 1:
        print("There are multiple edge attributes! The first attribute %s is used" % g.es.attributes()[0])

    ## an edge attribute is missing
    if len(g.es.attributes()) == 0:
        g.es["label"] = 1

    e_attr_name = g.es.attributes()[0]
    e_attr_values = np.asarray(g.es[e_attr_name]).reshape(len(g.es),1)
    E = np.hstack((E,e_attr_values))

    #if len(g.vs.attributes()) > 1:
    #   print("There are multiple vertex attributes! The first attribute %s (or the label) is used" % g.vs.attributes()[0])

    if len(g.vs.attributes()) == 0:
        g.vs["label"] = 1

    # Default to using a 'label' attribute of the graph if it is
    # present. This also accounts for the case where there are 2
    # or more attributes.
    if 'label' in g.vs.attributes():
        v_attr_name = 'label'
    else:
        v_attr_name = gs.vs.attributes()[0]

    v_attr_values = np.asarray(g.vs[v_attr_name]).reshape(len(g.vs),1).astype(int)

    res = dict([('edge', E), ('vlabel', v_attr_values), ('vsize', len(g.vs)), ('esize', len(g.es)), ('maxdegree', g.maxdegree())])

    return res




##############################
## Given a list of graphs
## Extract graph information
## Convert in the desired input 
## for graphkernels package
##############################

def GetGKInput(G):

    E = gkCpy.VecMatrixXi()
    V_label = gkCpy.IntIntVector()
    V_count = gkCpy.IntVector()
    E_count = gkCpy.IntVector()
    D_max = gkCpy.IntVector()

    for i in range(len(G)):

        g_info = GetGraphInfo(G[i])
        E.append(g_info['edge'])
        V_label.append(gkCpy.IntVector(g_info['vlabel'].reshape(-1).tolist()))
        V_count.append(g_info['vsize'])
        E_count.append(g_info['esize'])
        D_max.append(g_info['maxdegree'])

    return E, V_label, V_count, E_count, D_max



def GetAdjMatList(G):

    adj_mat = gkCpy.VecMatrixXi()
    adj_list = gkCpy.IntIntIntVector()

    for i in range(len(G)):

        # get adjacency matrix of graph i
        am_cur = G[i].get_adjacency()
        am_cur = np.array(am_cur.data)
        adj_mat.append(am_cur)

        # get adjacency list of graph i
        al_cur = np.asarray(G[i].get_adjlist())
        
        for j in range(len(al_cur)):
            al_cur[j] = gkCpy.IntVector(al_cur[j])

        adj_list.append(gkCpy.IntIntVector(al_cur))

    return adj_mat, adj_list







