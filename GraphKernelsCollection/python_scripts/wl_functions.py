# Author: Elisabetta Ghisu

"""

- Script for computing the WL graph kernel

- Inputs: adjencency list, node labels, number of WL iterations

- Output: Fetaures and kernel matrices

"""


###########################
# --- IMPORT PACKAGES --- #
###########################

import numpy as np
from sklearn.feature_extraction import FeatureHasher
from sklearn.preprocessing import scale
import copy
from scipy.sparse import lil_matrix, csr_matrix, vstack


###########################
# --- WL GRAPH KERNEL --- #
###########################

"""

INPUTS
- ad_list: the graphs are stored in lists; each element is a list itself containing the adjancency lists
- node_label: list of lists, where each entry is a list of node labels for the graph
- h: iterations of WL

OUTPUTS
- K: List of kernel matrices at each iteration
- phi_list: list of features matrices

"""

def WL_compute(ad_list, node_label, h):

    # Total number of graphs in the dataset
    n = len(ad_list)
    
    # Total number of nodes in dataset: initialized as zero
    tot_nodes = 0


    # list of kernel matrices
    K = [0]*(h+1)
    # list of feature mtrices
    phi_list = [0] * (h+1)

    #total number of nodes in the dataset
    for i in range(n):
        tot_nodes = tot_nodes + int(len(ad_list[i]))

    
    #each column of phi will be the explicit feature representation for the graph j
    phi = lil_matrix((tot_nodes, n), dtype = np.uint32)

    # labels will be used to store the new labels
    labels = [0] * n

    #label lookup is a dictionary which will contain the mapping
    # from multiset labels (strings) to short labels (integers)
    label_lookup = {}

    # counter to create possibly new labels in the update step
    label_counter = 0

    # Note: here we are just renaming the node labels from 0,..,num_labels
    # for each graph
    for i in range(n):

        # copy the original labels
        l_aux = np.copy(node_label[i])

        # will be used to store the new labels
        labels[i] = np.zeros(len(l_aux), dtype = np.int32)

        # for each label in graph
        for j in range(len(l_aux)):
            l_aux_str = str(l_aux[j])

            # If the string do not already exist
            # then create a new short label
            if not label_lookup.has_key(l_aux_str):
                label_lookup[l_aux_str] = label_counter
                labels[i][j] = label_counter                                                              
                label_counter += 1
            else:
                labels[i][j] = label_lookup[l_aux_str]

            # node histograph of the new labels
            phi[labels[i][j],i] += 1

    L = label_counter
    print 'Number of original labels %d' %L

    #####################
    # --- Main code --- #
    #####################

    # Now we are starting with the first iteration of WL

    # features obtained from the original node (renamed) labels
    phi_list[0] = phi

    # Kernel matrix based on original features
    K[0] = phi.transpose().dot(phi).toarray().astype(np.float32)
    
    print "K original is computed"
    
    # Initialize iterations to 0
    it = 0

    # copy of the original labels: will stored the new labels
    new_labels = np.copy(labels)
    
    # until the number of iterations is less than h
    while it < h:

        # Initialize dictionary and counter 
        # (same meaning as before)        
        label_lookup = {}
        label_counter = 0

        # Initialize phi as a sparse matrix
        phi = lil_matrix((tot_nodes, n), dtype = np.int32)
        # convert it to array
        phi = phi.toarray()

        print "Iteration %d: phi is computed" % it

        # for each graph in the dataset
        for i in range(n):

            # will store the multilabel string
            l_aux_long = np.copy(labels[i])

            # for each node in graph
            for v in range(len(ad_list[i])):

                # the new labels convert to tuple
                new_node_label = tuple([l_aux_long[v]]) 

                # form a multiset label of the node neighbors 
                new_ad = np.zeros(len(ad_list[i][v]))
                for j in range(len(ad_list[i][v])):
                    new_ad[j] = ad_list[i][v][j]
            
                ad_aux = tuple([l_aux_long[j] for j in new_ad])

                # long labels: original node plus sorted neughbors
                long_label = tuple(tuple(new_node_label)+tuple(sorted(ad_aux)))
            
                # if the multiset label has not yet occurred , add
                # it to the lookup table and assign a number to it
                if not label_lookup.has_key(long_label):
                    label_lookup[long_label] = str(label_counter)
                    new_labels[i][v] = str(label_counter)
                    label_counter += 1

                # else assign it the already existing number
                else:
                    new_labels[i][v] = label_lookup[long_label]

            # count the node label frequencies
            aux = np.bincount(new_labels[i]) 
            phi[new_labels[i],i] += aux[new_labels[i]]
        
        L = label_counter
        print 'Number of compressed labels %d' %L

        # create phi for iteration it+1
        phi_sparse = lil_matrix(phi)
        phi_list[it+1] = phi_sparse

        print "Itaration %d: phi sparse saved" % it

        # create K at iteration it+1
        K[it+1] = K[it] + phi_sparse.transpose().dot(phi_sparse).toarray().astype(np.float32)
       
        print "Iteration %d: K is computed" % it

        # Initialize labels for the next iteration as the new just computed
        labels = copy.deepcopy(new_labels)

        # increment the iteration
        it = it + 1 

    return K, phi_list



