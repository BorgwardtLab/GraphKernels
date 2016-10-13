# Description

This folder contains Python scripts to perform the Shortest Path and Weisfeiler Lehman Graph kernels. The data are assumed to be in Graphml format, as provided in the data folder.  

# Usage

- The files parse_graph, sp_functions and wl_functions contains various functions for analysing graphs and computing shortest path and WL kernels, respectively.

- The files compute_sp.py and compute_wl.py are scripts for computing the kernels and feature maps. The arguments --dataset and eventually number of WL iterations --h, need to be passed as a command line arguments. Example usage: python compute_wl.py --dataset mutag --h 3 will compute kernel and features map on the mutag dataset up to 3 iterations of WL.

- The files compute_perf_gk.py and compute_perf_gk_mlp.py compute the prediction performances using an SVM and MLP respectively. It is assumed that either kernel matrices or feature maps have been previously calculated, for instance using the compute_wl and compute_sp scripts. It is also possible to provide regression datasets, in this case a KernelRidge is used instead of an SVM. The performances are reported in terms of accuracy and root mean squared error, for calssification and regression respectively. Various arguments need to be passed as command line. For instance, you can choose which graph kernel to use (SP or WL) and the number of iterations of WL. A more detailed description is provided in the scripts.
