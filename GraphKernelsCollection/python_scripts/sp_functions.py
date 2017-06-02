# Author: Elisabetta Ghisu

"""

- Script containing functions for computing the shortest path kernel

- The Floyd Warshall algorithm is first implemented

- Then the SP is calculated

"""


#######################
# - IMPORT PACKAGES - #
#######################


import numpy as np
import numpy.matlib as matlib

"""
### FLOYD WARSHALL ALGORITHM

Input:
- Adjancency matrix A

Output:
- Shortest path matrix S

"""

def floyd_warshall(A):

	# nuber of nodes
	n = A.shape[0]

	# initialize shortes path matrix
	S = np.zeros(shape = (n,n))

	for i in range(n):
		for j in range(n):
			if A[i,j] == 0 and i!=j:
				S[i,j] = float("inf")
			else:
				S[i,j] = A[i,j]

	# Compute the shortest path matrix
	for k in range(n):
		for i in range(n):
			for j in range(n):
				if S[i,j] > S[i,k] + S[k,j]:
					S[i,j] = S[i,k] + S[k,j]

	return S								



"""
SHORTEST PATH KERNEL: This is a fast implementation of the shortest path
kernel algorithm

Inputs
- Adjancency matrix
- List of list of node labels for each graph
- Total number of node labels 

Outputs
- Kernel matrix
- Feature matrix

"""

def sp_kernel_fast(adj_mat, labels, L):

	# Number of graphs
	n = len(adj_mat)
	L = int(L)
	S = []

	# shortest path matrices
	for i in xrange(n):
		if i%1000 == 0 and i !=0:
			print "%d" % i
		S.append(floyd_warshall(adj_mat[i]))
	
	# maximum length of shortest paths in the dataset
	max_path = 0

	# for each graph in dataset
	for i in xrange(n):

		S_cur = np.copy(S[i])
		S_cur[S_cur == np.inf] = 0
		new_max = np.max(S_cur)
		
		if new_max > max_path:
			max_path = new_max

	# maximum length of shortest paths
	max_path = int(max_path)

	# initialize feature matrix
	sp = np.zeros(((max_path + 1) * L * (L+1) /2,n))

	# compute feature map for shortest path
	for i in xrange(n):

		if i % 1000 == 0:
			print "Processed %d graphs" %i

		S_graph = S[i]
		labels_graph = np.asarray(labels[i].reshape((len(labels[i]),1)))
		labels_graph = labels_graph + 1
		
		labels_aux = matlib.repmat(labels_graph, 1, len(labels_graph))
		
		min_lab = np.minimum(labels_aux, labels_aux.T)
		
		max_lab = np.maximum(labels_aux, labels_aux.T)
		sub_path = np.triu(~(np.isinf(S_graph))).T

		min_lab = min_lab[sub_path]
		max_lab = max_lab[sub_path]


		ind = S_graph[sub_path] * L * (L + 1) / 2 + (min_lab - 1) * (2*L + 2 - min_lab) / 2 + max_lab - min_lab
		ind = ind.astype(int)
		accum = np.zeros((max_path + 1) * L * (L + 1) /2)
		accum[:ind.max() + 1] += np.bincount(ind.astype(int))
		sp[ind,i] = accum[ind]
	
	sum_cols = np.sum(sp, axis = 1)
	ind_true = sum_cols != 0
	sp = sp[ind_true,:]
	
	# compute kernel matrix
	K = np.dot(sp.T,sp)
	
	return K, sp



