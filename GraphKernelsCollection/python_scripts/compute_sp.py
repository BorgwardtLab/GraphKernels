# Author: Elisabetta Ghisu

"""

- Script for computing the kernel matrix and features map 
using shortest path kernel

"""

###########################
# --- IMPORT PACKAGES --- #
###########################

import numpy as np
import argparse
import os
import pickle

from numpy import genfromtxt

from sp_functions import *
from parse_graphs import *



##############################
### Command Line Arguments ###
##############################

parser = argparse.ArgumentParser(description = "Compute kernel and features matrices via shortest path kernel")
parser.add_argument("--dataset", required = True, help = "Name of the dataset")
args = parser.parse_args()


#####################
### LOAD THE DATA ###
#####################

"""

- Here we load the data input and targets
- The data are assumed to be in graph formats
- They should be in graphml format 

"""

# path to the list of graphs and dataset
filenames_graphs = "data/%s.list" % (args.dataset)
path_to_dataset = "data/%s" % (args.dataset) 

# Load the targets
filename_labels = "data/%s_label.txt" % (args.dataset)

# load and process graphs
node_label, ad_list, G, Y = load_and_process(filenames_graphs, filename_labels, path_to_dataset)

# output directory
out_path = "kernel_matrices/%s/sp" % args.dataset

# If the output directory does not exist, then create it
if not os.path.exists(out_path):
    os.makedirs(out_path)


#########################
# --- SHORTEST PATH --- #
#########################


# assign labels starting from zero to the nodes
L, labels = rename_nodes(node_label)


# Compute adjancency matrix 
adj_mat = get_adj_mat(G)

# Compute kernel and feature maps using shortest path
K, phi = sp_kernel_fast(adj_mat, labels, L)

# save kernel matrix
file_name = "%s/%s_ker_mat" % (out_path, args.dataset)
np.save(file_name, K)

# save feature map
file_name = "%s/%s_phi_map" % (out_path, args.dataset)
np.save(file_name, phi)


