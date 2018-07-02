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

Alternatively, the package can be build from source. After downloading the source code from GitHub

    $ git clone https://github.com/BorgwardtLab/GraphKernels.git

users can use the `setup.py` script to install the package

    $ cd GraphKernels/src/graphkernels
    $ python setup.py build
    $ python setup.py install --user


In case of error in any of the above steps, please make sure that all requirements are satisfied. The install requirements section below provide instruction to install all the dependencies, in case you don't have them in your environment.

You should also make sure that you're installing the latest release of our package, in case you've had a previous version installed. To make sure the extension and package are not taken from your cache, you can use the `--no-cache-dir` option and install the package as:

`$ pip --no-cache-dir install graphkernels`

# Installing the requirements

Note that graphkernels is a Python library relying on C++ source code. The wrapper is built upon an extension  which need to be installed in order for the package to work (automatically handled). If you experience problems with the installation above, you might be missing one or more of the dependencies tools, which need to be pre-installed in your environment.  

- a C++ compiler (e.g. gcc, http://gcc.gnu.org/install/, XCode)
- eigen3 (http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)
- pkg-config (https://www.freedesktop.org/wiki/Software/pkg-config/)

We recommend the following steps for installing the dependencies

1) Install eigen3
    
    On Ubuntu:
    
    `$ sudo apt-get install libeigen3-dev`
    
    On MacOSX: 
    
    `$ brew install eigen`
    
2) Install pkg-config

    On Ubuntu:
    
    `$ sudo apt-get install pkg-config`
    
    On MacOSX
    
    `$ brew install pkg-config`
  
Additional Python dependencies are automatically handled when installing the `graphkernels` package: 

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
