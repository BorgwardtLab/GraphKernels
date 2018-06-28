# graphkernels 

`graphkernels` is a Python package for computing various graph kernels. 

For the C++ implementation and R package, please refer to https://github.com/mahito-sugiyama/graph-kernels.

The Python and R packages are described at:

- M. Sugiyama, M.E. Ghisu, F. Llinares-López and K. Borgwardt. **graphkernels: R and Python packages for graph comparison**. Bioinformatics, 2017. 

The paper can be found [here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btx602/4209994) 

# Installation

The users can installl the pakage via pip, by typing in a terminal

    $ pip install graphkernels 

Alternatively, the package can be build from source. After downloading the source code from GitHub or [pypi](https://pypi.python.org/pypi/graphkernels/0.1.2), users can use
the `setup.py` script to install the package, by typing:

    $ python setup.py build
    $ python setup.py install

# Requirements

Note that graphkernels is a Python library relying on C++ source code. The wrapper is built upon an extension  which need to be installed in order for the package to work. The installation described above will automatically install all the requirements for you and place them in the correct path. If you experience problems with the installation above, please note that the following tools are necessary for the package to run. 

- SWIG
- GKextCPy
- a C++ compiler (e.g. gcc, XCode)
- igraph
- numpy

# Usage

We provide a short [tutorial](https://github.com/eghisu/GraphKernels/tree/master/graphkernelsTutorial) for the basic usage of our package; there, you can also find an example script for computing graph kernels through our package on a benchmark dataset. 
