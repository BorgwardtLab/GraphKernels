
"""
Converting mutag to list format
"""

from igraph import Graph
import numpy as np

mutag_list = []

n_graphs = 188
for i in xrange(n_graphs):

	g_cur = Graph.Read_GraphML("/home/eghisu/projects/gk_python_wrapper/gk_python_c/data/mutag/mutag_%d.graphml" % (i+1))
	mutag_list.append(g_cur)


np.save("/home/eghisu/projects/gk_python_wrapper/gk_python_c/data/mutag_pydata", mutag_list)