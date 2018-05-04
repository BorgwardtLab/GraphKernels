# uncompyle6 version 3.1.3
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.6.3 (default, Oct  3 2017, 21:45:48) 
# [GCC 7.2.0]
# Embedded file name: graphkernels/graphkernels.py
# Compiled at: 2017-06-02 08:12:19
"""
Collect all the functions
for computing the graph kernels
"""
import numpy as np, GKextCPy as gkCpy
from utilities import GetGKInput

def CalculateEdgeHistKernel(G, par=-1.0):
    E, V_label, V_count, E_count, D_max = GetGKInput(G)
    par = gkCpy.DoubleVector([par])
    K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 1)
    return K


def CalculateVertexHistKernel(G, par=-1.0):
    E, V_label, V_count, E_count, D_max = GetGKInput(G)
    par = gkCpy.DoubleVector([par])
    K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 2)
    return K


def CalculateVertexEdgeHistKernel(G, par=-1.0):
    E, V_label, V_count, E_count, D_max = GetGKInput(G)
    par = gkCpy.DoubleVector([par])
    K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 3)
    return K


def CalculateVertexVertexEdgeHistKernel(G, par=1):
    E, V_label, V_count, E_count, D_max = GetGKInput(G)
    par = gkCpy.DoubleVector([par])
    K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 4)
    return K


def CalculateEdgeHistGaussKernel(G, par=1):
    E, V_label, V_count, E_count, D_max = GetGKInput(G)
    par = gkCpy.DoubleVector([par])
    K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 5)
    return K


def CalculateVertexHistGaussKernel(G, par=1):
    E, V_label, V_count, E_count, D_max = GetGKInput(G)
    par = gkCpy.DoubleVector([par])
    K = gk.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 6)
    return K


def CalculateVertexEdgeHistGaussKernel(G, par=1):
    E, V_label, V_count, E_count, D_max = GetGKInput(G)
    par = gkCpy.DoubleVector([par])
    K = gkCpy.CalculateKernelPy(E, V_label, V_count, E_count, D_max, par, 7)
    return K


def CalculateWLKernel(G, par=5):
    E, V_label, V_count, E_count, D_max = GetGKInput(G)
    par = int(par)
    K = gkCpy.WLKernelMatrix(E, V_label, V_count, E_count, D_max, par)
    return K
# okay decompiling graphkernels.pyc
