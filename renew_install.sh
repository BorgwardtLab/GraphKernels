pip uninstall graphkernels
pip uninstall GKextCPy

cd GKextCPy
swig -c++ -python GKextCPy.i
python setup.py build_ext --user
pip install . --user 

cd ../graphkernels
pip install . --user
