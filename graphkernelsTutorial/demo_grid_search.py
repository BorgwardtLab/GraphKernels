#!/usr/bin/env python3

from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold

class KernelGridSearchCV:
    """
    A simple class for performing a grid search for kernel matrices with
    a cross-validation strategy. At present, the class interface follows
    the default interface of `scikit-learn`. However, the class is *not*
    yet inheriting from any base class.

    The class only supports scoring based on accuracy so far.
    """

    def __init__(self, clf, param_grid, n_folds, cv = None, random_state = None, refit = True):
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
                    n _splits = 3,
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
