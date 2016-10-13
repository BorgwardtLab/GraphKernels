#Authors: Elisabetta Ghisu, Felipe Llinares Lopez

"""

- This script includes a list of functions for analyzing 
parsing and formatting graphs

- The graphs are given in graphml format

- It also cntans functions for loading, processing the graphs
and extract graph statistics

"""


import numpy as np
from numpy import genfromtxt

# iGraph imports to handle graphs and for graph I/O
from igraph import Graph


# ---------------------------------GRAPHML I/O FUNCTIONS------------------------------------ #

# INPUT:
# filenames_graphs: list of GraphML files, where each file contains one graph in the dataset
# filename_labels: text file with labels corresponding to each graph in the dataset, in the same order as they are in
#                  filename_graphs
# OUTPUT:
# G: A list containing one iGraph object for each graph in the dataset
# Y: A Numpy array containing the labels corresponding to each graph, in the same order as G
def load_graphml(filenames_graphs, filename_labels):
    G = []
    for fname in filenames_graphs:
        G.append(Graph.Read_GraphML(fname))
    Y = genfromtxt(filename_labels)
    return (G, Y)


# Loads a list of paths to GraphML files from filename_list
def load_file_list(filename_flist):
    f = open(filename_flist, 'r')
    f_graphs = []
    for line in f:
        f_graphs.append(line.strip())
    f.close()
    return f_graphs


# --------------------------------COMPUTE STATISTICS---------------------------------------- #


# Retrieve labels of all vertices belonging to any graph in the list of iGraph objects G and
# returns the entire list, and a list with the alphabet of the vertex labels
def get_all_vertex_labels(G, att_name='label'):
    v_l = []
    for g in G:
        v_l += g.vs[att_name]
    return (v_l, np.unique(v_l))


# Retrieve labels of all edges belonging to any graph in the list of iGraph objects G and
# returns the entire list, and a list with the alphabet of the edge labels
def get_all_edge_labels(G, att_name='label'):
    e_l = []
    for g in G:
        e_l += g.es[att_name]
    return (e_l, np.unique(e_l))


# Returns a list where each element is itself the adjacency list of the corresponding graph
# The adjacency lit of a graph has the following format:
# it is a list where each element is a list containing the id of adjacent nodes
def get_adj_list(G):
    ad_l = []
    for g in G:
        ad_l.append(g.get_adjlist())
    return ad_l

# Returns a list where each element is the adjacency matrix of the graph 
# The adjancency matrix is in iGraph format
def get_adj_mat(G):
    ad_m = []
    for g in G:
        ad_m.append(g.get_adjacency())
    return ad_m

# Returns a list where each element contains the nodes label for a graph
def get_node_labels(G, att_name = 'label'):
    node_l = []
    for g in G:
        node_l.append(g.vs[att_name])
    return node_l



# ----------------- LOAD AND PROCESS THE GRAPHS --------------- #


"""
Inputs:
- list of graphs file
- labels file
- path to the data folder

Outputs:
- List of node labels
- List of adjancency lists
- List of graphs in graphml format
- Targets
- number of classes
- sample size
"""


def load_and_process(filenames_graphs, filename_labels, path_to_dataset):

    # load a list of names to graphml files
    f_graphs = load_file_list(filenames_graphs)
    # sample size
    n = len(f_graphs)

    # create a list of paths to the files
    f_graphs_path =[]

    # for each graph in dataset
    for i in range(n):

        # index the graph
        graph_name = f_graphs[i]

        # path to the data folder
        path = "%s/%s" % (path_to_dataset, graph_name)
        f_graphs_path.append(path)

    # If the data is DD have to delete an element (corrupted file)
    if graph_name == "DD":
        del f_graphs_path[148]
        n = n-1

    # Load the graphs in graphml format
    # G is a llist of graphml graph
    # Y is an array of targets
    G,Y = load_graphml(f_graphs_path, filename_labels)

    # Delete corrupted file in DD
    if graph_name == "DD": 
        Y = np.delete(Y, 148)

    # get adjacency list and matrix for all the graphs in G
    ad_list = get_adj_list(G)
    ad_mat = get_adj_mat(G)

    # get a list containing lists of node labels
    node_label = get_node_labels(G)

    return node_label, ad_list, G, Y



"""

RENAME NODES: function to rename nodes from 0,...,num_nodes

Input
- list of list of node labels in each graph

Output
- L: total number of different labels in the dataset
- node_label: new renamed labels

"""

def rename_nodes(node_label): 
    
    # number of graphs in the dataset
    n = len(node_label)

    # labels will store the new labels
    labels = [0] * n

    # disctionary containing the map from the old to the new labels
    label_lookup = {}

    # counter of unique labels
    label_counter = 0

    # for each graph in dataset
    for i in range(n):


        # number of nodes in graph[i]
        num_nodes = len(node_label[i]) 

        # will be used to store the new labels
        labels[i] = np.zeros(num_nodes, dtype = np.uint64) # positive integers

        # for each node in the graph
        for j in range(num_nodes):

            # the node label to a string
            l_node_str = str(np.copy(node_label[i][j]))
            
            # if the string has not been observed yet
            # the corresponding node is assigned a new label
            # otherwise it will be named with the same label
            # already assigned to an identical string

            if not label_lookup.has_key(l_node_str):
                label_lookup[l_node_str] = label_counter
                labels[i][j] = label_counter                                                              
                label_counter += 1
            else:
                labels[i][j] = label_lookup[l_node_str]

    # total number of labels in the dataset
    L = label_counter
    print 'Number of original labels %d' % L 

    return L, labels

