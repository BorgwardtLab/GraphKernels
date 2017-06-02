#!/usr/bin/env python

"""
setup.py file for SWIG 
"""

from distutils.core import setup, Extension

graphkernelsCpy_module = Extension('_graphkernelsCpy', sources=['graphkernelsCpy_wrap.cxx', 'graphkernelsCpy.cpp'],swig_opts=['-c++'], extra_compile_args = ["-std=c++11"])
                           

setup (name = 'graphkernelsCpy',
       version = '0.2',
       author      = "Elisabetta Ghisu",
       description = """Graph Kernels package""",
       ext_modules = [graphkernelsCpy_module],
       py_modules = ["graphkernelsCpy"],
       license = 'ETH Zurich',
       )