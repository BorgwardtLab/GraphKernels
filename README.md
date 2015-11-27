# GraphKernels

Code for implementing different Graph Kernels

# Summary

This repository contains Matlab scripts for computing several graph kernels for graphs with unlabeled or 
categorically labeled nodes and scripts for evaluating the SVM classification performance using 
precomputed kernels using cross-validation. In the folder "data" we provide 5 data sets of undirected labeled graphs in Matlab format for graph 
classification which can be used for testing the scripts. More information, code and data can be found on the following website. https://www.bsse.ethz.ch/mlcb/research/machine-learning/graph-kernels.html

# Usage

"help X.m" executed in Matlab will display instructions on the usage of X.m for any script X.m

# Implemented Graph Kernels

For unlabeled graphs the following implementations are available:
- the family of graphlet kernels (subfolders allgraphlets, connectedgraphlets, samplinggraphlets) from [6]
- the random walk kernel (RWkernel.m) from [7]
- the shortest path kernel (SPkernel.m) from [1]

For labeled graphs the following implementations are available:
- 3 kernels from the Weisfeiler-Lehman kernel family: the WL subtree (WL.m), the WL shortest path (WLspdelta.m), 
  and the WL edge (WLedge.m) kernels from [5]
- the labeled 3-node graphlet kernel (l3graphletkernel.m) - extension of an algorithm in [6]
- the labeled random walk kernel (lRWkernel.m) from [7]
- the labeled random walk kernel based on walks up to size p (untilpRWkernel.m) - extension of [2],[3]
- the labeled shortest path kernel (lSPkernel.m) from [1]
- the Ramon and Gaertner subtree kernel (RGkernel.m) from [4]

Note that all WL kernels also support unlabeled graphs.


# References
[1] K. M. Borgwardt and H.-P. Kriegel. 
    Shortest-path kernels on graphs. In Proceedings of the International Conference on Data Mining, 
    pages 74-81, 2005.

[2] T. Gaertner, P. A. Flach, and S. Wrobel. 
    On graph kernels: Hardness results and efficient alternatives. In Proceedings of the 16th Annual 
    Conference on Computational Learning Theory and 7th Kernel Workshop, pages 129-143. 
    ISBN 3-540-40720-0.

[3] H. Kashima, K. Tsuda, and A. Inokuchi. 
    Marginalized kernels between labeled graphs. In Proceedings of the 20th International Conference 
    on Machine Learning, Washington, DC, United States, 2003.

[4] J. Ramon and T. Gaertner. 
    Expressivity versus efficiency of graph kernels. In First International Workshop on Mining Graphs, 
    Trees and Sequences (held with ECML/PKDD'03), 2003.

[5] N. Shervashidze, P. Schweitzer, E. J. van Leeuwen, K. Mehlhorn, and K. M. Borgwardt. 
    Weisfeiler-lehman graph kernels. Journal of Machine Learning Research, 12:2539-2561, 2011.

[6] N. Shervashidze, S. V. N. Vishwanathan, T. Petri, K. Mehlhorn, and K. M. Borgwardt. 
    Efficient graphlet kernels for large graph comparison. In Proceedings of the International 
    Conference on Artificial Intelligence and Statistics, 2009.

[7] S. V. N. Vishwanathan, N. N. Schraudolph, I. R. Kondor, and K. M. Borgwardt. 
    Graph kernels. Journal of Machine Learning Research, 11:1201-1242, 2010.

# Contact
This repository is mantained by

Elisabetta Ghisu

ETH Zürich, D-BSSE, Machine Learning and Computational Biology Lab

elisabetta.ghisu@bsse.ethz.ch

For technical questions on the code you can refer to the corresponding authors

