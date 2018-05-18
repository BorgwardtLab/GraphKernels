#!/usr/bin/env python

"""
setup.py file for SWIG 
"""

#from distutils.core import setup, Extension
from setuptools import setup, Extension
import numpy as np

GKextCPy_module = Extension('_GKextCPy', sources=['GKextCPy_wrap.cxx', 'GKextCPy.cpp'],swig_opts=['-c++'], extra_compile_args = ["-std=c++11"])
                           

setup (name = 'GKextCPy',
       version = '0.3.1',
       author      = "Elisabetta Ghisu",
       description = """Graph Kernels: building the extension Python module. This is a wrapper package from C++ to Python.""",
       include_dirs = [np.get_include()],
       ext_modules = [GKextCPy_module],
       py_modules = ["GKextCPy"],
       license = 'ETH Zurich',
       )
