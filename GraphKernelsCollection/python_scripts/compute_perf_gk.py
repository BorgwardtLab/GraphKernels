# Author: Elisabetta Ghisu

"""

- This script take as input a kernel matrix
and returns the classification or regression performance

- The kernel matrix can be calculated using any of the graph kernels approaches

- The criteria used for prediction are SVM for classification and kernel Ridge regression for regression

- For predition we divide the data in training, validation and test. For each split, we first train on the train data, 
then evaluate the performance on the validation. We choose the optimal parameters for the validation set and finally
provide the corresponding perforance on the test set. If more than one split is performed, the final results 
correspond to the average of the performances on the test sets. 

"""

###########################
# --- IMPORT PACKAGES --- #
###########################

import numpy as np
import pickle
import os
import argparse
import random

from numpy import genfromtxt

from sklearn.kernel_ridge import KernelRidge # 0.17
from sklearn.metrics import accuracy_score, mean_squared_error
from sklearn import svm


##################################
# --- COMMAND LINE ARGUMENTS --- #
##################################

parser = argparse.ArgumentParser(description = "Classification/regression experiments with SP")
parser.add_argument("--graph_kernel", required = True, help = "Which graph kernel: sp or wl")
parser.add_argument("--dataset", required = True, help = "Name of the dataset")
parser.add_argument("--type", required = True, help = "Classification or regression")
parser.add_argument("--trials", required = True, help = "Trials for hyperparameters random search", type = int)
parser.add_argument("--splits", required = True, help = "number of splits", type = int)
parser.add_argument("--h", help = "WL iteration. Only required for WL", type = int)


args = parser.parse_args()

# Number of parameter trials
trials = "%d" % (args.trials) 
trials = int(trials)

# Set the seed for uniform parameter distribution
random.seed(20) 

# Regression or classification problem
model_type = "%s" % args.type 

# Number of splits of the data
splits = args.splits


#########################
# --- LOAD THE DATA --- #
######################### 

graph_kernel = "%s" % args.graph_kernel

# Path to the kernel matrices
if graph_kernel == "sp":
	ker_path = "kernel_matrices/%s/%s" % (args.dataset, args.graph_kernel)
elif graph_kernel == "wl":
	h = args.h
	ker_path = "kernel_matrices/%s/%s/h%d" % (args.dataset, args.graph_kernel, h)

filename = "%s/%s_ker_mat.npy" % (ker_path,args.dataset)

# Load the kernel matrix
print "Loading kernel matrix..."
K = np.load(filename)

# Path to the targets
label_name = "data/%s_label.txt" % (args.dataset)
# Load targets
y = genfromtxt(label_name)

# Write label from 0,..,5
if args.dataset == "enzymes":
	y -=1

# Size of the dataset
n = K.shape[0]


#################################
# --- SET UP THE PARAMETERS --- #
#################################

# You should modify these arguments
# depending on the range of parameters that you want to examine
alpha_grid = np.linspace(0.01, 100, num = trials)
C_grid = np.linspace(0.0001, 10, num = trials)


##############################################################
# --- MAIN CODE: PERMUTE, SPLIT AND EVALUATE PEFORMANCES --- #
##############################################################


"""

-  Here starts the main program

-  First we permute the data, then
for each split we evaluate corresponding performances

-  In the end, the performances are averaged over the test sets

"""

# Initialize the performance of the best parameter trial on validation
# With the corresponding performance on test
val_split = []
test_split = []

# For each split of the data
for j in xrange(10, 10 + splits):

	print "Starting split %d..." % j

	# Set the random set for data permutation
	random_state = int(j)
	np.random.seed(random_state)
	idx_perm = np.random.permutation(n)

	# Set the output path
	if graph_kernel == "sp":
		output_path = "output/kernel/%s/random_state_%d/" % (args.dataset,random_state)
	elif graph_kernel == "wl":
		output_path = "output/kernel/%s/h_%d/random_state_%d" % (args.dataset, h,random_state)

	# if the output directory does not exist, then create it
	if not os.path.exists(output_path):
		os.makedirs(output_path)

		
	# Permute the data
	y_perm = y[idx_perm] #targets permutation
	K_perm = K[:,idx_perm] #inputs permutation
	K_perm = K_perm[idx_perm,:] #inputs permutation

	# Set the training, validation and test
	# Note: the percentage can be set up by the user
	num_train_val = int((n * 90)/100)         #90% (of entire dataset) for training and validation
	num_test = n - num_train_val              #10% (of entire dataset) for test
	num_train = int((num_train_val * 90)/100) #90% (of train + val) for training
	num_val = num_train_val - num_train       # ~10% (of train + val) for validation

	# Split the kernel matrix
	K_train = K_perm[0:num_train,0:num_train]
	K_val = K_perm[num_train:(num_train+num_val),0:num_train]
	K_test = K_perm[(num_train+num_val):n, 0:num_train]

	# Split the targets
	y_train = y_perm[0:num_train]

	# Normalization step (for real valued targets only)
	if model_type == "regression":
		y_train_mean = np.mean(y_train)
		y_train_std = np.std(y_train)
		y_train = (y_train-y_train_mean)/float(y_train_std)

	y_val = y_perm[num_train:(num_train+num_val)]
	y_test = y_perm[(num_train+num_val):n]

	# Record the performance for each parameter trial
	# respectively on validation and test set
	perf_all_val = []
	perf_all_test = []


	#####################################################################
	# --- RUN THE MODEL: FOR A GIVEN SPLIT AND EACH PARAMETER TRIAL --- #
	#####################################################################

	# For each parameter trial
	for i in xrange(trials):

		# For regression use the Kernel Ridge method
		if model_type == "regression":

			print "\n Starting experiment for trial %d and parameter alpha = %3f\n " % (i, alpha_grid[i])

			# Fit the kernel ridge model
			KR = KernelRidge(kernel = 'precomputed', alpha = alpha_grid[i])
			KR.fit(K_train, y_train)

			# predict on the validation and test set
			y_pred = KR.predict(K_val)
			y_pred_test = KR.predict(K_test)
			
			# adjust prediction: needed because the training targets have been normalizaed
			y_pred = y_pred * float(y_train_std) + y_train_mean
			y_pred_test = y_pred_test * float(y_train_std) + y_train_mean

			# root mean squared error on validation
			rmse = np.sqrt(mean_squared_error(y_val, y_pred))
			perf_all_val.append(rmse)

			# root mean squared error in test 
			rmse_test = np.sqrt(mean_squared_error(y_test, y_pred_test))
			perf_all_test.append(rmse_test)

			print "The performance on the validation set is: %3f" % rmse
			print "The performance on the test set is: %3f" % rmse_test

		# For clcassification use SVM
		if model_type == "classification":

			print "\nStarting experiment for trial %d and parameter C = %3f \n\n" % (i, C_grid[i])

			# Fit classifier on training data
			clf = svm.SVC(kernel = 'precomputed', C = C_grid[i])
			clf.fit(K_train, y_train)

			# predict on validation and test
			y_pred = clf.predict(K_val)
			y_pred_test = clf.predict(K_test)

			# accuracy on validation set
			acc = accuracy_score(y_val, y_pred)
			perf_all_val.append(acc)

			# accuracy on test set
			acc_test = accuracy_score(y_test, y_pred_test)
			perf_all_test.append(acc_test)

			print "The performance on the validation set is: %3f" % acc
			print "The performance on the test set is: %3f" % acc_test


	#######################################
	# --- FIND THE OPTIMAL PARAMETERS --- #
	#######################################

	# For regression: minimise the mean squared error
	if model_type == "regression":

		# get optimal parameter on validation (argmin mean squared error)
		min_idx = np.argmin(perf_all_val)
		alpha_opt = alpha_grid[min_idx]

		# performance corresponding to optimal parameter on val
		perf_val_opt = perf_all_val[min_idx]

		# corresponding performance on test for the same parameter
		perf_test_opt = perf_all_test[min_idx]

		print "The best performance is for trial %d with parameter alpha = %3f" % (min_idx, alpha_opt)
		print "The best performance on the validation set is: %3f" % perf_val_opt
		print "The corresponding performance on test set is: %3f" % perf_test_opt

	# For classification: maximise the accuracy
	if model_type == "classification":

		# get optimal parameter on validation (argmax accuracy)
		max_idx = np.argmax(perf_all_val)
		C_opt = C_grid[max_idx]

		# performance corresponsing to the optimal parameter on validation
		perf_val_opt = perf_all_val[max_idx]

		# corresponding performance on the test set for the same parameter
		perf_test_opt = perf_all_test[max_idx]

		print "\n The best performance is for trial %d with parameter C = %3f" % (max_idx, C_opt)
		print "The best performance on the validation set is: %3f" % perf_val_opt
		print "The corresponding performance on test set is: %3f" % perf_test_opt


	#######################
	# --- SAVE RESULTS ---#
	#######################

	# file to save performances
	f_perf = "%s/kernel_perf_%s.txt" % (output_path,args.dataset)
	f_out = open(f_perf,'w')

	# Write on file: random state, iteration of WL, 
	# performances on validation and test set for each parameter trial
	f_out.write("Random state = %d\n" % random_state)

	if graph_kernel == "wl":
		f_out.write("h of WL = %d\n" % h)
	f_out.write("Performance on the validation and test set for different parameter trials \n \n")
	f_out.write("trial \t param \t \t perf val \t perf test\n")

	# For regression: alpha param and min_idx
	if model_type == "regression":

		for j in xrange(trials):
			f_out.write("%d \t %3f \t %3f \t %3f \n" % (j,  alpha_grid[j], perf_all_val[j], perf_all_test[j]))
		
		f_out.write("The best performance is for trial %d with parameter %3f" % (min_idx, alpha_opt))

	# For classification: C param and max_idx
	if model_type == "classification":

		for j in xrange(trials):
			f_out.write("%d \t %3f \t %3f \t %3f \n" % (j, C_grid[j], perf_all_val[j], perf_all_test[j]))
		
		f_out.write("\n\nThe best performance is for trial %d with parameter %3f" % (max_idx, C_opt))

	# Optimal performance on validation
	# and corresponding on test set
	f_out.write("\nThe best performance on the validation set is: %3f" % perf_val_opt)
	f_out.write("\nThe corresponding performance on test set is: %3f" % perf_test_opt)
	f_out.close()

	# save an array of all the performances on validation
	name_val = "%s/perf_val_%s" % (output_path,args.dataset)
	np.save(name_val, perf_all_val)

	# and correspondin performances on test set
	name_test = "%s/perf_test_%s" % (output_path,args.dataset)
	np.save(name_test, perf_all_test)

	# append the best performance on validation
	# at the current split
	val_split.append(perf_val_opt)

	# append the correponding performance on the test set
	test_split.append(perf_test_opt)


###############################
# --- AVERAGE THE RESULTS --- #
###############################

# mean of the validation performances over the splits
val_mean = np.mean(np.asarray(val_split))
# std deviation of validation over the splits
val_std = np.std(np.asarray(val_split))

# mean of the test performances over the splits
test_mean = np.mean(np.asarray(test_split))
# std deviation of the test oer the splits
test_std = np.std(np.asarray(test_split))

print "\n Mean performance on val set: %3f" % val_mean
print "With standard deviation: %3f" % val_std
print "\n Mean performance on test set: %3f" % test_mean
print "With standard deviation: %3f" % test_std

# output path for the averaged results
if graph_kernel == "sp":
	output_path = "output/kernel/%srandom_state_avg" % (args.dataset)
elif graph_kernel == "wl":
	output_path = "output/kernel/%s/h_%d/random_state_avg" % (args.dataset, h)

#if the output directory does not exist, then create it
if not os.path.exists(output_path):
	os.makedirs(output_path)

# filename of the final averaged results
name_all_split = "%s/kernel_perf_avg_%s.txt" % (output_path,args.dataset)
f_out = open(name_all_split,'w')

# Write down: number of averaged splits,
# iteration of WL,
# mean performances on validation and test,
# with respectively standard deviations
f_out.write("Average of performances over %d splits" % splits)

if graph_kernel == "wl":
	f_out.write("\nh of WL = %d\n" % h)
f_out.write("\nThe best mean performance on the validation set is: %3f" % val_mean)
f_out.write("\nThe corresponding standard deviation is: %3f" % val_std)
f_out.write("\nThe best mean performance on the test set is: %3f" % test_mean)

# For classification: report the performance in 1-acc
# i.e. the lower the best
if model_type == "classification":
	f_out.write("\nMean performance (1-acc) on the test set: %3f" % (1-test_mean))

f_out.write("\nThe corresponding standard deviation is: %3f" % test_std)
f_out.close()
	








