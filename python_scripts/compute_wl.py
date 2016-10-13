# Author: Elisabetta Ghisu

"""

- This script compute the feature and kernel matrices
obtained with the Weisfeiler Lehman Graph Kernel

- These can be computed for different iterations of WL 

- Choosing a given h, we get all the kernels and feature maps up to h

"""


###########################
# --- IMPORT PACKAGES --- # 
###########################

import numpy as np
import argparse
import os

from numpy import genfromtxt
from parse_graphs import load_and_process
from wl_functions import *


##############################
### COMMAND LINE ARGUMENTS ###
##############################

parser = argparse.ArgumentParser(description = "Classification/regression experiments with SP")
parser.add_argument("--dataset", required = True, help = "Name of the dataset")
parser.add_argument("--h", required = True, help = "WL iterations", type = int)
args = parser.parse_args()

# Number of iterations of WL
h = args.h


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


node_label, ad_list, G, Y = load_and_process(filenames_graphs, filename_labels, path_to_dataset)


# Apply WL graph kernel
# Get a list of h kernel matrices: K
# get a list of h features maps: phi 
K, phi = WL_compute(ad_list, node_label, h)


# For each iteration of WL
for j in xrange(h+1):

	# Path to the output folder
	out_path = "kernel_matrices/%s/wl/h%d" % (args.dataset, j)

	# If the output directory does not exist, then create it
	if not os.path.exists(out_path):
		os.makedirs(out_path)

	# save kernel matrix
	file_name = "%s/%s_ker_mat" % (out_path, args.dataset)
	np.save(file_name,K[j])

# For each iteration of WL
for j in xrange(h+1):
	out_path = "kernel_matrices/%s/wl/h%d" % (args.dataset, j)

	# If the output directory does not exist, then create it
	if not os.path.exists(out_path):
		os.makedirs(out_path)

	# save feature map
	file_name = "%s/%s_phi_map" % (out_path, args.dataset)
	np.save(file_name,phi[j])


