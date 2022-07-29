/*
Header file for the graphkernels module
*/
#include <stdint.h>
#define Int int32_t
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <stdio.h>
#include <string.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

typedef Eigen::Triplet<double> T;



using namespace Eigen;

using namespace std;

/*template <typename S>
ostream &operator<<(ostream& out, const vector<S>& vec) {
  for (Int i = 0; i < vec.size() - 1; ++i) { cout << vec[i] << " "; }
  cout << vec[vec.size() - 1];
  return out;
}*/

/*
template <class ForwardIterator, class T>
  void iota (ForwardIterator first, ForwardIterator last, T val)
{
  while (first!=last) {
    *first = val;
    ++first;
    ++val;
  }
}

*/

double selectLinearGaussian(vector<int>& h1, vector<int>& h2, double sigma);

int productMapping(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, MatrixXi& H);

MatrixXd productAdjacency(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, MatrixXi& H);

void bucketsort(vector<int>& x, vector<int>& index, int label_max);



//All the kernel functions 
double edgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, double sigma);
double vertexHistogramKernel(vector<int>& v1_label, vector<int>& v2_label, double sigma);
double vertexEdgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double sigma);
double vertexVertexEdgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double lambda);

double geometricRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double lambda);
double exponentialRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double beta);
double kstepRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, vector<double>& lambda_list);

MatrixXd WLKernelMatrix(vector<MatrixXi>& E, vector<vector<int> >& V_label, vector<int>& num_v, vector<int>& num_e, vector<int>& degree_max, int h_max);

double computeKernelValue(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, vector<double>& par, int kernel_type);

MatrixXd CalculateKernelPy(vector<MatrixXi>& E, vector<vector<int> >& V_label, vector<int>& V_count, vector<int>& E_count, vector<int>& D_max, vector<double>& par, int kernel_type);

// added for graphlets
int find_min(int a, int b, int c);

void card_ThreeInter(vector<int>& L1, vector<int>& L2, vector<int>& L3, vector<int>& card);
void getIndices(vector<int>& o_set1, vector<int>& o_set2, vector<int>& inter, vector<int>& diff1, vector<int>& diff2);


VectorXd countGraphletsFour(vector<vector<int>>& al, VectorXd& count_gr);
void getCardinality(vector<int>& o_set1, vector<int>& o_set2, vector<double>& card);
VectorXd countGraphletsThree(vector<vector<int>>& al, VectorXd& count_gr);
void getMinValue(MatrixXi& iam, vector<int>& idx, vector<int>& sums);

VectorXd countConnectedGraphletsFive(MatrixXi& am, vector<vector<int>>& al, VectorXd& count_gr);

VectorXd countConnectedGraphletsFour(MatrixXi& am, vector<vector<int>>& al, VectorXd& count_gr);
VectorXd countConnectedGraphletsThree(MatrixXi& am, vector<vector<int>>& al, VectorXd& count_gr);


MatrixXd CalculateGraphletKernelPy(vector<MatrixXi>& graph_adj_all, vector<vector<vector<int>>>& graph_adjlist_all, int k);
MatrixXd CalculateConnectedGraphletKernelPy(vector<MatrixXi>& graph_adj_all, vector<vector<vector<int>>>& graph_adjlist_all, int k);






