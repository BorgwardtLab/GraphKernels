

# graphkernels 

graphkernels is a Python package for computing various graph kernels

# Installation

The users can installl the pakage via pip, by typing in a terminal

pip install graphkernels 

Alternatively, the package can be build from source. After downloading the source code from guthub or pypi (https://pypi.python.org/pypi/graphkernels/0.1.2), users can use the setup.py script to install the package, by typing

python setup.py build
python setup.py install 

# Requirements

Note that graphkernels is a Python library relying on C++ source code. The wrapper is built upon an extension  which need to be installed in order for the package to work. The installation described above will automatically install all the requirements for you and place them in the correct path. If you experience problems with the installation above, please note that the following tools are necessary for the package to run. 

- SWIG
- GKextCPy
- a C++ compiler (e.g. gcc, XCode)
- igraph
- numpy

# Usage

Please have a look at https://github.com/eghisu/GraphKernels/tree/master/graphkernelsTutorial for a short tutorial on the basic usage of our package. Here, we provide an example script for computing graph kernels through our package on a benchmark dataset. 
