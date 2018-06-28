#!/usr/bin/env python3

from sklearn.base import clone
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
from sklearn.model_selection import ParameterGrid
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC

import graphkernels.kernels as gk
import igraph as ig
import numpy as np

import os
import sys

class KernelGridSearchCV:
    """
    A simple class for performing a grid search for kernel matrices with
    a cross-validation strategy. At present, the class interface follows
    the default interface of `scikit-learn`. However, the class is *not*
    yet inheriting from any base class.

    The class only supports scoring based on accuracy so far.
    """

    def __init__(self, clf, param_grid, cv = None, random_state = None, refit = True):
        self.clf_             = clf
        self.grid_            = param_grid
        self.cv_              = cv
        self.random_state_    = random_state
        self.refit_           = refit
        self.best_estimator_  = None
        self.best_score_      = None

    def fit(self, X, y):

        # Use stratified k-folds with a user-specified number
        if self.cv_ is None:
            cv = KFold(
                    n_splits = 3,
                    shuffle = True,
                    random_state = self.random_state_
            )
        elif isinstance(self.cv_, int):
            cv = StratifiedKFold(
                    n_splits = self.cv_,
                    shuffle = True,
                    random_state = self.random_state_
            )
        else:
            cv = self.cv_

        grid = ParameterGrid(self.grid_)
        for parameters in grid:
            clf = self.clf_
            clf.set_params(**parameters)

            scores = []
            for train, test in cv.split(np.zeros(len(y)), y):
                X_train = X[train][:, train]
                y_train = y[train]
                X_test  = X[test][:, train]
                y_test  = y[test]

                clf.fit(X_train, y_train)

                # The class only supports the accuracy score for now.
                ac = accuracy_score(y_test, clf.predict(X_test))
                scores.append(ac)

            score = np.mean(scores)
            if self.best_score_ is None or score > self.best_score_:
                self.best_estimator_ = clone(clf)
                self.best_score_     = score
                self.best_params_    = parameters

if __name__ == '__main__':
    data_directory = 'mutag'
    X = []
    y = np.genfromtxt('mutag.label')
    with open('mutag.list') as f:
        for line in f:
            filename = line.strip()
            full_path = os.path.join(data_directory, filename)
            X.append(ig.Graph.Read_GraphML(full_path))

    # Calculate kernel matrix for a few kernels. We are only using a few
    # kernels with a good runtime performance. Please refer to the other
    # tutorials and documentation for a list of *all* available kernels.

    kernels = {
        'vertex_histogram': gk.CalculateVertexHistKernel,
        'edge_histogram': gk.CalculateEdgeHistKernel,
        'weisfeiler_lehman': gk.CalculateWLKernel
    }

    kernel_matrices = dict()
    for kernel_name, kernel_function in sorted(kernels.items()):
        print('Calculating full kernel matrix for', kernel_name)

        # This delegates the kernel calculation to the kernel function
        # specified above. In a real-world scenario, you should supply
        # more parameters to each kernel, such as `sigma` for Gaussian
        # kernels (histograms), or the number of iterations for the WL
        # kernel. Here, we just use default values, and search for the
        # best parameters in an SVM.
        kernel_matrices[kernel_name] = kernel_function(X)

    for kernel_name, kernel_matrix in sorted(kernel_matrices.items()):

        # Prepare a parameter grid to train an SVM classifier with
        # cross-validation.
        grid = {
                'C': 10. ** np.arange(-2,3)
        }

        # Calling the classifier like this ensures that we can use
        # our kernel matrices from above.
        clf = SVC(kernel='precomputed')

        grid_search = KernelGridSearchCV(
            clf,
            param_grid = grid,
            cv = 10, # 10-fold cross-validation; the interface also
                     # supports your own cross-validator script for
                     # more specific cases

            # Make this tutorial reproducible
            random_state = 42,
        )

        # This runs the parameter search over the specified grid
        # according to the cross-validation that was selected in
        # the above code.
        grid_search.fit(kernel_matrix,y)

        # The grid search stores the best estimator, i.e. the one that
        # yielded the best accuracy. You can use it for subsequent ops
        # as well.
        clf = grid_search.best_estimator_

        print('10-fold cross-validation for {} yields an accuracy of {:2.2f}'.format(kernel_name, grid_search.best_score_ * 100.0))
