from setuptools import setup

setup(name = 'graphkernels',
      version = '0.2',
	description = 'Package for computing graph kernels',
	url = 'http://github.com/eghisu/graphkernels',
	author = 'Elisabetta Ghisu',
	author_email = 'elisabetta.ghisu@bsse.ethz.ch',
	license = 'ETH Zurich',
	packages = ['graphkernels'],
	install_requires = ['graphkernelsCpy'],
	zip_safe = False) # maybe to check github url
