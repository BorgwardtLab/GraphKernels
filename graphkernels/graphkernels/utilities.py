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
	for i in xrange(len(g.es)):
		E[i,:] = g.es[i].tuple
	## there are multiple edge attributes
  	if len(g.es.attributes()) > 1:
  		print "There are multiple edge attributes! The first attribute %s is used" % g.es.attributes()[0]

  	## an edge attribute is missing
 	if len(g.es.attributes()) == 0:
  		g.es["label"] = 1

  	e_attr_name = g.es.attributes()[0]
  	e_attr_values = np.asarray(g.es[e_attr_name]).reshape(len(g.es),1)
  	E = np.hstack((E,e_attr_values))

  	#if len(g.vs.attributes()) > 1:
  	#	print "There are multiple vertex attributes! The first attribute %s is used" % g.vs.attributes()[0]

 	if len(g.vs.attributes()) == 0:
  		g.vs["label"] = 1

  	v_attr_name = g.vs.attributes()[1]
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


	for i in xrange(len(G)):

		g_info = GetGraphInfo(G[i])
		E.append(g_info['edge'])
		V_label.append(gkCpy.IntVector(g_info['vlabel'].reshape(-1).tolist()))
		V_count.append(g_info['vsize'])
		E_count.append(g_info['esize'])
		D_max.append(g_info['maxdegree'])


	return E, V_label, V_count, E_count, D_max



