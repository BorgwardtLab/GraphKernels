# GraphKernels

A large collection of source code for computing kernels on graph. The src folder contains the following elements:

- graphkernels is a Python package for graph kernels. The Python interface is created from a C++ source code that is wrapped with SWIG (http://www.swig.org)
- GKextCPy is the package created to build the extension module for the wrapper from C++ to Python

# graphkernels 

`graphkernels` is a Python package for computing various graph kernels. 

For the C++ implementation and R package, please refer to https://github.com/mahito-sugiyama/graph-kernels.

The Python and R packages are described at:

- M. Sugiyama, M.E. Ghisu, F. Llinares-López and K. Borgwardt. **graphkernels: R and Python packages for graph comparison**. Bioinformatics, 2017. 

The paper can be found [here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btx602/4209994) 

# Installation

The users can installl the pakage via pip, by typing in a terminal

    $ pip install graphkernels 

Alternatively, the package can be build from source. After downloading the source code from GitHub, users can use
the `setup.py` script to install the package, by typing:

    $ python setup.py build
    $ python setup.py install

# Requirements

Note that graphkernels is a Python library relying on C++ source code. The wrapper is built upon an extension  which need to be installed in order for the package to work. The installation described above will automatically install all the requirements for you and place them in the correct path. If you experience problems with the installation above, please note that the following tools are necessary for the package to run. 

- a C++ compiler (e.g. gcc, http://gcc.gnu.org/install/, XCode)
- eigen3 (http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)
- pkg-config (https://www.freedesktop.org/wiki/Software/pkg-config/)

Please, note that the dependencies above need to be pre-installed in your machine. The links above provide information and instructions for the installation. The above tools can also be easily installed via homebrew, for instance:

`brew install eigen`

Python dependenies - automatically handled:
- GKextCPy
- igraph
- numpy


# Usage

We provide a short [tutorial](https://github.com/eghisu/GraphKernels/tree/master/Tutorial) for the basic usage of our package; there, you can also find an example script for computing graph kernels through our package on a benchmark dataset. 

# Citation

If you use the `graphkernels` package in your projects please cite our work

```
@article{Sugiyama-2017-Bioinformatics,
author = {Sugiyama, Mahito and Ghisu, M. Elisabetta and Llinares-López, Felipe and Borgwardt, Karsten},
title = {graphkernels: R and Python packages for graph comparison},
journal = {Bioinformatics},
volume = {34},
number = {3},
pages = {530--532},
year = {2017},
doi = {10.1093/bioinformatics/btx602},
URL = {http://dx.doi.org/10.1093/bioinformatics/btx602},
}
```
