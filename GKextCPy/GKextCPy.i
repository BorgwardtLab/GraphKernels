

/* Define the SWIG module: graphkernels */

%module GKextCPy


%{
/* line for specifying that the C file should be built as a python extension */
# include <memory>
#include <vector>
# define SWIG_FILE_WITH_INIT 
#include <Python.h>
# define ITPP_EXPORT
# include "GKextCPy.h"
%}

//include te built-in support for std::vector
%include <typemaps.i>
%include <std_vector.i>
%include <std_string.i>
%include <std_shared_ptr.i>
%include <numpy.i>
%include <eigen.i>
%include <stl.i>

%init %{
import_array();
%}


using namespace std;

namespace std {
    %template(IntVector) vector<int>;
    %template(IntIntVector) vector<vector<int> >;
    %template(FloatVector) vector<float>;
    %template(DoubleVector) vector<double>;
    %template(VecMatrixXi) vector<Eigen::MatrixXi>;
};

double selectLinearGaussian(vector<int>& h1, vector<int>& h2, double sigma);
int productMapping(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, Eigen::MatrixXi& H);

Eigen::MatrixXd productAdjacency(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, Eigen::MatrixXi& H);

void bucketsort(vector<int>& x, vector<int>& index, int label_max);

//All the kernel functions 
double edgeHistogramKernel(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, double sigma);
double vertexHistogramKernel(vector<int>& v1_label, vector<int>& v2_label, double sigma);
double vertexEdgeHistogramKernel(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double sigma);
double vertexVertexEdgeHistogramKernel(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double lambda);

double geometricRandomWalkKernel(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double lambda);
double exponentialRandomWalkKernel(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double beta);
double kstepRandomWalkKernel(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, vector<double>& lambda_list);

Eigen::MatrixXd WLKernelMatrix(vector<Eigen::MatrixXi>& E, vector<vector<int> >& V_label, vector<int>& num_v, vector<int>& num_e, vector<int>& degree_max, int h_max);

double computeKernelValue(Eigen::MatrixXi& e1, Eigen::MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, vector<double>& par, int kernel_type);

Eigen::MatrixXd CalculateKernelPy(vector<Eigen::MatrixXi>& E, vector<vector<int> >& V_label, vector<int>& V_count, vector<int>& E_count, vector<int>& D_max, vector<double>& par, int kernel_type);



/*
void getGraphInfo(List& graph_info_list, vector<MatrixXi>& E, vector<vector<Int>>& V_label, vector<Int>& V_count, vector<Int>& E_count, vector<Int>& D_max);



void CalculateKernel(List& graph_info_list, vector<double>& par, Int kernel_type);
*/
