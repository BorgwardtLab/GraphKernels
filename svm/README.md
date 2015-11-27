# Description

The folder contains scripts for evaluating the classification performance of SVM with precomputed kernels 
using cross-validation. You need to install libSVM and its Matlab interface to be able to use scripts in this folder.
The two main functions are:

- runmultiplesvm.m cross-validate over different kernels and the parameter C of SVM and repeat the experiment n times with different folds. 
This is used with the output of all WL (Weisfeiler Lehman) scripts, a cell array of kernel matrices. 
To be able to use it, you have to replace "~/code/libsvm" in runsvm.m by the name of the folder where you put the compiled mex files 
(svmtrain and svmpredict) of the libSVM Matlab interface.

- runntimes.m cross-validate over the parameter C of SVM and repeat the experiment n times with different folds.
This is used with a single kernel matrix. To be able to use it, you have to replace "~/code/libsvm" 
in runIndependent.m by the name of the folder where you put the compiled mex files (svmtrain and svmpredict) of the 
libSVM Matlab interface.

# Usage

To obtain information on input and output, run "help runmultiplesvm" and "help runntimes" respectively.
