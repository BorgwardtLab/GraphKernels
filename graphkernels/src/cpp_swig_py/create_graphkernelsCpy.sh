### Bach script to create the C++ Graphkernels package

swig -c++ -python graphkernelsCpy.i # compile the swig file
python setup.py build_ext --user # to build locally
python setup.py install --user ## to install the package (locally)
