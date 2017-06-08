from setuptools import setup

setup(name = 'graphkernels',
      version = '0.1.2',
	description = 'Package for computing graph kernels',
	url = 'https://github.com/eghisu/GraphKernels/tree/master/graphkernels',
	author = 'Elisabetta Ghisu',
	author_email = 'elisabetta.ghisu@bsse.ethz.ch',
	license = 'ETH Zurich',
	packages = ['graphkernels'],

	install_requires = ['GKextCPy'],
	package_data={'graphkernels': ['data.mutag']},
	) # maybe to check github url
