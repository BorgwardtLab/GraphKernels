# Author: Elisabetta Ghisu

"""

- Script for predicting performances using and MLP
and the input features representation of the graph_kernel

- The MLP is computed using the tensorflow library

"""


###########################
# --- IMPORT PACKAGES --- #
###########################


import numpy as np
import argparse
import random
import os

import IPython as ip

from sklearn.metrics import mean_squared_error, accuracy_score

import tensorflow as tf


##################################
# --- COMMAND LINE ARGUMENTS --- #
##################################

parser = argparse.ArgumentParser(description = "Perform an MLP with features from WL or SP kernels")
parser.add_argument("--graph_kernel", required = True, help = "Which graph kernel: sp or wl")
parser.add_argument("--dataset", required = True, help = "Name of the dataset")
parser.add_argument("--type", required = True, help = "Classification or regression")
parser.add_argument("--trials", required = True, help = "Trials for hyperparameters random search", type = int)
parser.add_argument("--splits", required = True, help = "number of splits", type = int)
parser.add_argument("--iterations", required = True, help = "Number of neural network iterations", type = int)
parser.add_argument("--h", help = "WL iteration. Only required for WL", type = int)

args = parser.parse_args()

graph_kernel = args.graph_kernel
h = args.h

# number of splits of the data
splits = args.splits

# iterations of the neural network
num_iter = "%d" % (args.iterations)
num_iter = int(num_iter)

# number of parameters trials
trials = "%d" % (args.trials) 
trials = int(trials) 

# set the seed for the distribution of the parameters
random.seed(20)

# classification or regression
model_type = "%s" % args.type 

batch_size = 30

#########################
# --- LOAD FEATURES --- #
#########################


# Function for creating the feature map
# needed for wl only
def create_phi(phi_names):

	phi_old = np.load(phi_names[0]).reshape(1,1)[0,0].toarray()
	phi = phi_old

	if len(phi_names)>0:
		for i in xrange(1,len(phi_names)):

			print "Loading features for h = %d" % i
			phi_cur = np.load(phi_names[i]).reshape(1,1)[0,0].toarray()

			phi = np.concatenate((phi_old,phi_cur), axis = 0)
		
			phi_old = phi

	return phi


print "Loading feature maps..."

# Load feature map
if graph_kernel == "wl":
	# list of paths to feature matrices (up to h)
	phi_paths = []
	for i in xrange(h+1):
		path_name = "kernel_matrices/%s/%s/h%d" % (args.dataset, graph_kernel, i)
		phi_paths.append(path_name)

	phi_names = []
	for i in xrange(h+1):
		phi_cur = "%s/%s_phi_map.npy" % (phi_paths[i], args.dataset)
		phi_names.append(phi_cur)

	phi = create_phi(phi_names)

elif graph_kernel == "sp":
	ker_path = "kernel_matrices/%s/%s" % (args.dataset, graph_kernel)
	filename = "%s/%s_phi_map.npy" % (ker_path,args.dataset)
	phi = np.load(filename)


# print out information on phi
print "Shape of the feature matrix: n_features x n_samples"
print phi.shape

label_name = "/links/groups/borgwardt/Projects/NIPS16_DGK/data/data-graphml/%s_label.txt" % (args.dataset)
if args.dataset == "cep":
	label_name = "/links/groups/borgwardt/Projects/NIPS16_DGK/data/data-graphml/%s_single_feat_0_label.txt" % (args.dataset)
else:
	label_name = "/links/groups/borgwardt/Projects/NIPS16_DGK/data/data-graphml/%s_label.txt" % (args.dataset)

y = np.loadtxt(label_name)

# Write label from 0,..,5
if args.dataset == "enzymes":
	y -=1

#dataset size
n = phi.shape[1]
print "\nDataset size: %d" %n


### PARAMETERS GRID ###

nn_hdim = [random.randint(2,1000) for x in xrange(trials)] # hidden dimensions
learn_rate = [random.uniform(0,0.01) for x in xrange(trials)] # learning rate for gradient descent
reg_lambda = [random.uniform(0,0.01) for x in xrange(trials)] # regularization strength


####################################
# --- NEURAL NETWORK FUNCTIONS --- #
####################################


"""

INITIALIZE THE WEIGHTS

Input
- shape: shape of the WEIGHTS

Output
- weights as a random normal variable

"""

# Weights 
def weight_variable(shape):

	return tf.Variable(tf.random_normal(shape, seed = 5))


# Biases
def bias_variable(shape):
    initial = tf.random_normal(shape, seed = 2)
    return tf.Variable(initial)


"""

MLP model: function to model a basic mlp

"""

def mlp_model(X, w_hid, w_out, b_hid, b_out):

	# hidden layer
	hid = tf.nn.sigmoid(tf.add(tf.matmul(X, w_hid),b_hid))

	# output layer
	out  = tf.add(tf.matmul(hid, w_out), b_out)

	return out 



# Build the prediction model
# We use tensorflow for bilding
# the graph flow

def build_model(train_X, train_y, val_X, val_y, test_X, test_y, nn_hdim, nn_output_dim, learn_rate, reg_lambda, num_iter, batch_size, model_type):

	# Size definitions
	num_examples = train_X.shape[0]
	nn_input_dim = train_X.shape[1]

	# initialize tf graph
	X = tf.placeholder(tf.float32, [None, nn_input_dim])
	y = tf.placeholder(tf.float32, [None, nn_output_dim])


    # Initialize dimensions of the shared variables 
	w_hid = weight_variable([nn_input_dim, nn_hdim])
	w_out = weight_variable([nn_hdim, nn_output_dim])
	b_hid = bias_variable([nn_hdim])
	b_out = bias_variable([nn_output_dim])


	# prediction
	y_hat = mlp_model(X, w_hid, w_out, b_hid, b_out)
	
	# Define loss and optimizer
	w_sum = tf.nn.l2_loss(w_hid) + tf.nn.l2_loss(w_out)
	
	# regularization term
	#loss_reg = (reg_lambda/(float(2*batch_size)))*(w_sum)
	loss_reg = (reg_lambda/(float(2*train_X.shape[0])))*(w_sum)

	if model_type == "classification":
		cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(y_hat, y)) #+ loss_reg
	if model_type == "regression":
		cost = tf.sqrt(tf.reduce_sum(tf.pow(y_hat-y,2))/(2*batch_size) + loss_reg)
        
	optimizer = tf.train.GradientDescentOptimizer(learning_rate=learn_rate).minimize(cost)
    
    # Initialize the variables 
	init = tf.initialize_all_variables()

	val_cost = []
	test_cost = []
	val_perf = []
	test_perf = []

	# Launch the graph in a session
	with tf.Session() as sess:

		# run session with initialized variables
		sess.run(init)

		for i in xrange(num_iter):

			for start, end in zip(range(0,len(train_X), batch_size), range(batch_size, len(train_X)+1, batch_size)):

				_,c = sess.run([optimizer, cost], feed_dict = {X: train_X[start:end], y: train_y[start:end]})

			# display epoch loss
			if i % 10 == 0 or i == (i - 1):

				print "\n Epoch:", '%04d' % (i), "cost=","{:.9f}".format(c)

				va_cost = sess.run(cost, feed_dict = {X:val_X, y:val_y})
				val_cost.append(va_cost)
				print "Validation cost:", va_cost

				te_cost = sess.run(cost, feed_dict = {X:test_X, y:test_y})
				test_cost.append(te_cost)
				print "Test cost:", te_cost

				if model_type == "classification":

					pred = tf.argmax(y_hat, 1)

					va_pred = pred.eval({X:val_X, y:val_y})
					va_acc = accuracy_score(np.argmax(val_y,1), va_pred)
					val_perf.append(va_acc) 

					print "Accuracy on validation set:", va_acc

					te_pred = pred.eval({X:test_X, y:test_y})
					te_acc = accuracy_score(np.argmax(test_y,1), te_pred)
					test_perf.append(te_acc)

					print "Accuracy on test set:", te_acc

				if model_type == "regression":

					pred = y_hat

					va_pred = pred.eval({X:val_X, y:val_y})
					va_rmse = np.sqrt(mean_squared_error(val_y, va_pred))
					val_perf.append(va_rmse) 

					print "RMSE on validation set:", va_rmse

					te_pred = pred.eval({X:test_X, y:test_y})
					te_rmse = np.sqrt(mean_squared_error(test_y, te_pred))
					test_perf.append(te_rmse)  

					print "RMSE on test set:", te_rmse

	return val_cost, test_cost, val_perf, test_perf


split_val_perf = []
split_test_perf = []

# For each slit of the data
for k in xrange(10, 10 + splits):

	print "Starting split %d..." % k


	# --- PERMUTE AND SPLIT THE DATA --- #

	# set the random seed for data shuffling
	random_state = int(k)
	np.random.seed(random_state)

	# indeces permutation
	idx_perm = np.random.permutation(n)

	# permute targets
	y_perm = y[idx_perm]

	# permute features
	phi_perm = phi[:, idx_perm] 

	# Set number of training, validation and test
	# Note: This can be modified
	num_train_val = int((n * 90)/100) # 90% training and validation
	num_test = n - num_train_val # 10% test
	num_train = int((num_train_val * 90)/100) #90% for training
	num_val = num_train_val - num_train # ~10% of 90 for validation

	# Define trainin validation and test set
	phi_train = phi_perm[:,0:num_train].T
	phi_val = phi_perm[:,num_train:(num_train+num_val)].T
	phi_test = phi_perm[:,(num_train+num_val):n].T

	train_X = phi_train.astype(np.float32)
	val_X = phi_val.astype(np.float32)
	test_X = phi_test.astype(np.float32)


	# Output is the vector of target for regression
	# and one hot encoding of classes in classification
	if model_type == "regression":
		nn_output_dim = 1 #output layer dimensionality
	if model_type == "classification":
		nn_output_dim = len(np.unique(y))

	print "\nNumber of output dimensions (classes for classification): %d" % int(nn_output_dim)

	if model_type == "regression":
		y_onehot = y_perm.reshape(n,nn_output_dim).astype(np.float32)
	if model_type == "classification":
		y_onehot = np.zeros(shape = (n,nn_output_dim)).astype(np.float32)
		for i in xrange(len(y_onehot)):
			for j in xrange(nn_output_dim):
				if y_perm[i] == j:
					y_onehot[i,j] = 1

	train_y = y_onehot[0:num_train]
	val_y = y_onehot[num_train:(num_train+num_val)]
	test_y = y_onehot[(num_train+num_val):n]

	if graph_kernel == "wl":
		output_path = "output/mlp/%s/%s/h%d/random_state_%d" % (args.dataset, graph_kernel, h, random_state)
	elif graph_kernel == "sp":
		output_path = "output/mlp/%s/%s/random_state_%d"% (args.dataset, graph_kernel, random_state)
	
	# If the output directory does not exist, then create it
	if not os.path.exists(output_path):
		os.makedirs(output_path)


	## Run the model for different parameters

	trial_val_cost = []
	trial_test_cost = []
	trial_val_perf = []
	trial_test_perf = []

	for i in xrange(trials):

		nn_hdim_cur = np.int(nn_hdim[i])
		learn_rate_cur = np.float32(learn_rate[i])
		reg_lambda_cur = np.float32(reg_lambda[i])

		print "\nStart experiment for parameter trial %d:" % i
		print "hidden dimensions: %d" % nn_hdim_cur 
		print "learning rate: %3f" % learn_rate_cur
		print "reg lambda: %3f\n" % reg_lambda_cur

		# Compute performances with an mlp
		val_cost, test_cost, val_perf, test_perf = build_model(train_X, train_y, val_X, val_y, test_X, test_y, 
			nn_hdim_cur, nn_output_dim, learn_rate_cur, reg_lambda_cur, num_iter, batch_size, model_type)

		file_name_trial = "%s/%s_trial_%d.txt" % (output_path, args.dataset, i)
		f_out_trial = open(file_name_trial, 'w')

		f_out_trial.write("Random state: %d" % random_state)
		f_out_trial.write("\n Trial: %d" % i)
		f_out_trial.write("\nhidden dimensions: %d \nlearn_rate: %3f \nreg lambda: %3f \n\n" % (nn_hdim_cur, learn_rate_cur, reg_lambda_cur))
		f_out_trial.write("iter \t val loss \t test loss \t val perf \t test perf \n")

		# for each epoch write results on txt file
		for j in xrange(len(val_cost)):
			f_out_trial.write("%d \t %.3f \t\t %.3f \t\t %.3f \t\t %.3f \n" % (j*10, val_cost[j], test_cost[j], val_perf[j], test_perf[j]))	


		# Now I choose the best performnce on validation set an report the corresponding result on the test set
		# This would be the equivalent of early stopping

		if model_type == "regression":
			idx_trial = np.argmin(val_perf)
			print "\n Minimum RMSE on validation for trial %d reached after %d epochs: %.3f" % (i, idx_trial * 10, val_perf[idx_trial])
			print " Corresponding RMSE on test set: %.3f" % test_perf[idx_trial]

		elif model_type == "classification":
			idx_trial = np.argmax(val_perf)
			print "\n Maximum accuracy on validation for trial %d reached after %d epoch: %.3f" % (i, idx_trial * 10, val_perf[idx_trial])
			print " Corresponding accuracy on test set: %.3f" % test_perf[idx_trial]

		trial_val_cost.append(val_cost[idx_trial])
		trial_test_cost.append(test_cost[idx_trial])
		trial_val_perf.append(val_perf[idx_trial])
		trial_test_perf.append(test_perf[idx_trial])

		#f_out_trial.write("\n Minimum validation cost for trial %d reached after %d epoch: %.3f" % (i, idx_trial*10, val_cost[idx_trial]))
		f_out_trial.write("\n Optimal performance on validation set: %.3f" % val_perf[idx_trial])	
		f_out_trial.write("\n Corresponding performance on test set: %.3f" % test_perf[idx_trial])
		f_out_trial.close()


	if model_type == "classification":
		idx_opt = np.argmax(trial_val_perf)
	elif model_type == "regression":
		idx_opt = np.argmin(trial_val_perf)


	file_name = "%s/%s_results.txt" % (output_path, args.dataset)
	f_out = open(file_name, 'w')

	f_out.write("Results corresponding to the best performance on validation set for each trial \n")

	for i in xrange(trials): 
		f_out.write("\n Optimal performance for trial %d and parameters:\n" % (i))
		f_out.write("\nhidden dimensions: %d \nlearn_rate: %3f \nreg lambda: %3f \n\n" % (nn_hdim[i], learn_rate[i], reg_lambda[i]))
		f_out.write("val loss \t test loss \t\t val perf \t\t test perf \n")
		f_out.write("%.3f \t\t %.3f \t\t %.3f \t\t %.3f \n " % (trial_val_cost[i], trial_test_cost[i], trial_val_perf[i], trial_test_perf[i]))


	f_out.write("\n \n Optimal performance is for parameters: \n")
	f_out.write("\nhidden dimensions: %d \nlearn_rate: %3f \nreg lambda: %3f \n\n" % (nn_hdim[idx_opt], learn_rate[idx_opt], reg_lambda[idx_opt]))
	f_out.write("\n Loss on validation set: %.3f " % trial_val_cost[idx_opt])
	f_out.write("\n Loss on test set: %.3f" % trial_test_cost[idx_opt])
	f_out.write("\n Performance on validation set: %.3f" % trial_val_perf[idx_opt])
	f_out.write("\n Performance on test set: %.3f" % trial_test_perf[idx_opt])

	if model_type == "classification":
		f_out.write("\n Performance (1-acc) on validation set: %.3f" % (1-trial_val_perf[idx_opt]))
		f_out.write("\n Performance (1-acc) on test set: %.3f" % (1-trial_test_perf[idx_opt]))

	f_out.close()
	print("\n \n Optimal performance is for parameters:\n")
	print("\nhidden dimensions: %d \nlearn_rate: %3f \nreg lambda: %3f \n\n" % (nn_hdim[idx_opt], learn_rate[idx_opt], reg_lambda[idx_opt]))
	print("\n Loss on validation set: %.3f " % trial_val_cost[idx_opt])
	print("\n Loss on test set: %.3f" % trial_test_cost[idx_opt])
	print("\n Performance on validation set: %.3f" % trial_val_perf[idx_opt])
	print("\n Performance on test set: %.3f" % trial_test_perf[idx_opt])

	split_val_perf.append(trial_val_perf[idx_opt])
	split_test_perf.append(trial_test_perf[idx_opt])




###########################################
### --- Average over all the splits --- ###
###########################################

val_mean = np.mean(np.asarray(split_val_perf))
val_std = np.std(np.asarray(split_val_perf))
test_mean = np.mean(np.asarray(split_test_perf))
test_std = np.std(np.asarray(split_test_perf))


print "\n Mean performance on val set: %3f" % val_mean
print "With standard deviation: %3f" % val_std
print "\n Mean performance on test set: %3f" % test_mean
print "With standard deviation: %3f" % test_std

if graph_kernel == "wl":
	output_path = "output/mlp/%s/%s/h%d/random_state_avg" % (args.dataset, graph_kernel, h)
elif graph_kernel == "sp":
	output_path = "output/mlp/%s/%s/random_state_avg"% (args.dataset, graph_kernel)
	

# If the output directory does not exist, then create it
if not os.path.exists(output_path):
    os.makedirs(output_path)

file_name = "%s/%s_avg_results.txt" % (output_path, args.dataset)

f_out = open(file_name,'w')
f_out.write("Average of performances over %d splits" % splits)
f_out.write("\nThe best mean performance on the validation set is: %3f" % val_mean)
f_out.write("\nThe corresponding standard deviation is: %3f" % val_std)
f_out.write("\nThe best mean performance on the test set is: %3f" % test_mean)
if model_type == "classification":
	f_out.write("\nMean performance (1-acc) on the test set: %3f" % (1-test_mean))
f_out.write("\nThe corresponding standard deviation is: %3f" % test_std)
f_out.close()



