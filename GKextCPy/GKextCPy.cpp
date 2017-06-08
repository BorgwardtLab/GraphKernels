
#define Int int32_t
#include "GKextCPy.h"

using namespace Eigen;
typedef Eigen::Triplet<double> T;

using namespace std;

// =================================================================== //
// ==================== Functions used in kernels ==================== //
// =================================================================== //
// select linear kernel or Gaussian kernel in histogram kernels
double selectLinearGaussian(vector<int>& h1, vector<int>& h2, double sigma) {
  double K = 0;
  if (sigma < 0) {
    // linear kernel
    for (int i = 0; i < (int)h1.size(); i++) {
      K += (double)h1[i] * (double)h2[i];
    }
  } else {
    // Gaussian kernel
    for (int i = 0; i < (int)h1.size(); i++) {
      K += ((double)h1[i] - (double)h2[i]) * ((double)h1[i] - (double)h2[i]);
    }
    K = exp(-1.0 * K / (2.0 * sigma * sigma));
  }
  return K;
}
// map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
int productMapping(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, MatrixXi& H) {
  int n_vx = 0;
  for (int i = 0; i < (int)v1_label.size(); i++) {
    for (int j = 0; j < (int)v2_label.size(); j++) {
      if (v1_label[i] == v2_label[j]) {
	H(i, j) = n_vx;
	n_vx++;
      }
    }
  }
  return n_vx;
}





//compute the adjacency matrix Ax of the direct product graph (sparse)
MatrixXd productAdjacency(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, MatrixXi& H) {
  
  int n_vx = v1_label.size()*v2_label.size();
  
  SparseMatrix<double> Ax(n_vx, n_vx);
  MatrixXd dAx;

  vector<T> v;

  
  for (int i = 0; i < e1.rows(); i++) {
    for (int j = 0; j < e2.rows(); j++) {      
      if (   v1_label[e1(i, 0)] == v2_label[e2(j, 0)]
    && v1_label[e1(i, 1)] == v2_label[e2(j, 1)]
    && e1(i, 2) == e2(j, 2)) {
  v.push_back(T(H(e1(i, 0), e2(j, 0)), H(e1(i, 1), e2(j, 1)), 1.0));
  v.push_back(T(H(e1(i, 1), e2(j, 1)), H(e1(i, 0), e2(j, 0)), 1.0));
    }
      if (   v1_label[e1(i, 0)] == v2_label[e2(j, 1)]
    && v1_label[e1(i, 1)] == v2_label[e2(j, 0)]
    && e1(i, 2) == e2(j, 2)) {
  v.push_back(T(H(e1(i, 0), e2(j, 1)), H(e1(i, 1), e2(j, 0)), 1.0));
  v.push_back(T(H(e1(i, 1), e2(j, 0)), H(e1(i, 0), e2(j, 1)), 1.0));
      }
    }
 }

  Ax.setFromTriplets(v.begin(), v.end());
  dAx = MatrixXd(Ax);
  
  return dAx;
}



// bucket sort used in Weisfeiler-Leiman graph kernel
void bucketsort(vector<int>& x, vector<int>& index, int label_max) {
  vector<vector<int> > buckets;
  buckets.resize(label_max + 1);

  for (vector<int>::iterator itr = index.begin(), end = index.end(); itr != end; ++itr) {
    buckets[ x[*itr] ].push_back(*itr);
  }

  int counter = 0;
  for (vector<vector<int> >::iterator itr = buckets.begin(), end = buckets.end(); itr != end; ++itr) {
    for (vector<int>::iterator itr2 = (*itr).begin(), end2 = (*itr).end(); itr2 != end2; ++itr2) {
      index[counter] = *itr2;
      counter++;
    }
  }
}




// ===================================================== //
// ==================== Each kernel ==================== //
// ===================================================== //
// edge histogram karnel
double edgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, double sigma) {
  int e_label_max = 0;
  for (int i = 0; i < e1.rows(); i++) {
    if (e1(i, 2) > e_label_max) e_label_max = e1(i, 2);
  }
  for (int i = 0; i < e2.rows(); i++) {
    if (e2(i, 2) > e_label_max) e_label_max = e2(i, 2);
  }

  vector<int> h1(e_label_max + 1, 0);
  vector<int> h2(e_label_max + 1, 0);

  for (int i = 0; i < e1.rows(); i++) {
    (h1[e1(i, 2)])++;
  }
  for (int i = 0; i < e2.rows(); i++) {
    (h2[e2(i, 2)])++;
  }

  return selectLinearGaussian(h1, h2, sigma);
}



// vertex histogram karnel
double vertexHistogramKernel(vector<int>& v1_label, vector<int>& v2_label, double sigma) {
  int v1_label_max = *max_element(v1_label.begin(), v1_label.end());
  int v2_label_max = *max_element(v2_label.begin(), v2_label.end());
  int v_label_max = v1_label_max > v2_label_max ? v1_label_max : v2_label_max;

  vector<int> h1(v_label_max + 1, 0);
  vector<int> h2(v_label_max + 1, 0);

  for (int i = 0; i < (int)v1_label.size(); i++) {
    (h1[v1_label[i]])++;
  }
  for (int i = 0; i < (int)v2_label.size(); i++) {
    (h2[v2_label[i]])++;
  }

  return selectLinearGaussian(h1, h2, sigma);
}


// vertex-edge histogram karnel
double vertexEdgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double sigma) {
  int e_label_max = 0;
  for (int i = 0; i < e1.rows(); i++) {
    if (e1(i, 2) > e_label_max) e_label_max = e1(i, 2);
  }
  for (int i = 0; i < e2.rows(); i++) {
    if (e2(i, 2) > e_label_max) e_label_max = e2(i, 2);
  }
  e_label_max++;

  int v1_label_max = *max_element(v1_label.begin(), v1_label.end());
  int v2_label_max = *max_element(v2_label.begin(), v2_label.end());
  int v_label_max = v1_label_max > v2_label_max ? v1_label_max : v2_label_max;
  v_label_max++;

  vector<int> h1(v_label_max * v_label_max * e_label_max, 0);
  vector<int> h2(v_label_max * v_label_max * e_label_max, 0);

  int v1, v2;
  for (int i = 0; i < e1.rows(); i++) {
    v1 = e1(i, 0);
    v2 = e1(i, 1);
    if (v2 > v1) { int v_tmp = v1; v1 = v2; v2 = v_tmp; }
    (h1[v1_label[v1] + v1_label[v2] * v_label_max + e1(i, 2) * v_label_max * v_label_max])++;
  }
  for (int i = 0; i < e2.rows(); i++) {
    v1 = e2(i, 0);
    v2 = e2(i, 1);
    if (v2 > v1) { int v_tmp = v1; v1 = v2; v2 = v_tmp; }
    (h2[v2_label[v1] + v2_label[v2] * v_label_max + e2(i, 2) * v_label_max * v_label_max])++;
  }

  return selectLinearGaussian(h1, h2, sigma);
}

// vertex-vertex-edge histogram karnel
double vertexVertexEdgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double lambda) {
  return vertexHistogramKernel(v1_label, v2_label, -1.0) + lambda * vertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, -1.0);
}






// geometric random walk karnel
double geometricRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double lambda) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());  
  int n_vx = productMapping(e1, e2, v1_label, v2_label, H);

  // prepare identity matrix
  SparseMatrix<double> I(n_vx, n_vx);
  I.setIdentity();

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  MatrixXd dAx(n_vx, n_vx);
  
  dAx = productAdjacency(e1, e2, v1_label, v2_label, H);
  Ax = dAx.sparseView();

  // inverse of I - lambda * Ax by fixed-poInt iterations
  VectorXd I_vec(n_vx);
  for (int i  = 0; i < n_vx; i++) I_vec[i] = 1;
  VectorXd x = I_vec;
  VectorXd x_pre(n_vx); x_pre.setZero();

  double eps = pow(10, -10);
  int count = 0;
  while ((x - x_pre).squaredNorm() > eps) {
    if (count > 100) {
      // cout << "does not converge until " << count - 1 << " iterations" << endl;
      break;
    }
    x_pre = x;
    x = I_vec + lambda * Ax * x_pre;
    count++;
  }
  return x.sum();
}


// exponential random walk karnel
double exponentialRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, double beta) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());  
  int n_vx = productMapping(e1, e2, v1_label, v2_label, H);

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  MatrixXd dAx(n_vx, n_vx);

  dAx = productAdjacency(e1, e2, v1_label, v2_label, H);
  Ax = dAx.sparseView();

  // compute e^{beta * Ax}
  SelfAdjointEigenSolver<MatrixXd> es(Ax);
  VectorXd x = (beta * es.eigenvalues()).array().exp();
  MatrixXd D = x.asDiagonal();  
  MatrixXd V = es.eigenvectors();

  MatrixXd I(n_vx, n_vx);
  I.setIdentity();
  FullPivLU<MatrixXd> solver(V);
  MatrixXd V_inv = solver.solve(I);
  MatrixXd Res = V * D * V_inv;

  // compute the total sum
  double K = 0;
  for (int i = 0; i < Res.rows(); i++) {
    for (int j = 0; j < Res.cols(); j++) {
      K += Res(i, j);
    }
  }

  return K;
}


// k-step product graph karnel
double kstepRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, vector<double>& lambda_list) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());
  int n_vx = productMapping(e1, e2, v1_label, v2_label, H);

  // prepare identity matrix
  SparseMatrix<double> I(n_vx, n_vx);
  I.setIdentity();

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  MatrixXd dAx(n_vx, n_vx);

  dAx = productAdjacency(e1, e2, v1_label, v2_label, H);
  Ax = dAx.sparseView();

  // compute products until k
  int k_max = (int)lambda_list.size() - 1;
  SparseMatrix<double> Ax_pow = I;
  SparseMatrix<double> Sum = lambda_list[0] * I;
  for (int k = 1; k <= k_max; k++) {
    Ax_pow = Ax * Ax_pow;
    Sum += lambda_list[k] * Ax_pow;
  }

  // compute the total sum
  double K = 0;
  for (int i = 0; i < Sum.outerSize(); ++i) {
    for (SparseMatrix<double>::InnerIterator it(Sum, i); it; ++it) {
      K += it.value();
    }
  }

  return K;
}



// Weisfeiler-Leiman graph kernel
MatrixXd WLKernelMatrix(vector<MatrixXi>& E, vector<vector<int> >& V_label, vector<int>& num_v, vector<int>& num_e, vector<int>& degree_max, int h_max) {
  // K_mat.setZero();

  MatrixXd K_mat(V_label.size(), V_label.size());


  int n = (int)E.size();
  int v_all = accumulate(num_v.begin(), num_v.end(), 0);
  int degree_max_all = *max_element(degree_max.begin(), degree_max.end());
  vector<int> label_max_vec(n);
  for (int i = 0; i < n; i++) {
    label_max_vec[i] = *max_element(V_label[i].begin(), V_label[i].end());
  }
  int label_max = *max_element(label_max_vec.begin(), label_max_vec.end());

  int raise = 0;
  vector<int> counter(*max_element(num_v.begin(), num_v.end()));
  MatrixXi nei_list(v_all, degree_max_all + 1);
  MatrixXi label_list(v_all, h_max + 1);
  vector<int> x(v_all);
  vector<int> index(v_all);
  vector<int> index_org(v_all);
  vector<int> graph_index(v_all);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < num_v[i]; j++) {
      label_list(j + raise, 0) = V_label[i][j];
      graph_index[j + raise] = i;
    }
    raise += num_v[i];
  }

  // ===== Increment kernel values using the initial vertex labels =====
  // radix sort
  for (int i = 0; i < v_all; i++) {
    index[i] = i;
    index_org[i] = i;
    x[i] = label_list(i, 0);
  }
  bucketsort(x, index, label_max);
  // add kernel values
  vector<int> count(n);
  set<int> count_index;
  int k_value;
  for (int i = 0; i < v_all; i++) {
    count_index.insert(graph_index[index_org[index[i]]]);
    count[graph_index[index_org[index[i]]]]++;
    if (i == v_all - 1 || label_list(index[i], 0) != label_list(index[i + 1], 0)) {
      for (set<Int>::iterator itr = count_index.begin(), end = count_index.end(); itr != end; ++itr) {
	for (set<int>::iterator itr2 = itr, end2 = count_index.end(); itr2 != end2; ++itr2) {
	  k_value = count[*itr] * count[*itr2];
	  K_mat(*itr, *itr2) += k_value;
	  K_mat(*itr2, *itr) += k_value;
	}
	count[*itr] = 0;
      }
      count_index.clear();
    }
  }

  int v_raised_1, v_raised_2;
  for (int h = 0; h < h_max; h++) {
    nei_list.setZero();

    // first put vertex label
    nei_list.col(0) = label_list.col(h);
    // second put neibor labels
    raise = 0;
    for (int i = 0; i < n; i++) {
      fill(counter.begin(), counter.end(), 1);
      for (int j = 0; j < num_e[i]; j++) {
	v_raised_1 = E[i](j, 0) + raise;
	v_raised_2 = E[i](j, 1) + raise;
	nei_list(v_raised_1, counter[E[i](j, 0)]) = label_list(v_raised_2, h);
	nei_list(v_raised_2, counter[E[i](j, 1)]) = label_list(v_raised_1, h);
	counter[E[i](j, 0)]++;
	counter[E[i](j, 1)]++;
      }
      raise += num_v[i];
    }

    // radix sort
    for (int i = 0; i < v_all; i++) {
      index[i] = i;
      index_org[i] = i;
    }
    for (int k = nei_list.cols() - 1; k >= 0; k--) {
      for (int i = 0; i < v_all; i++) {
	x[i] = nei_list(i, k);
      }
      bucketsort(x, index, label_max);
    }

    // re-labeling and increment kernel values
    label_max++;
    for (int i = 0; i < v_all; i++) {
      label_list(index_org[index[i]], h + 1) = label_max;
      count_index.insert(graph_index[index_org[index[i]]]);
      count[graph_index[index_org[index[i]]]]++;
      if (i == v_all - 1 ||
	  (nei_list.row(index[i]) - nei_list.row(index[i + 1])).sum() != 0) {
	for (set<int>::iterator itr = count_index.begin(), end = count_index.end(); itr != end; ++itr) {
	  for (set<int>::iterator itr2 = itr, end2 = count_index.end(); itr2 != end2; ++itr2) {
	    k_value = count[*itr] * count[*itr2];
	    K_mat(*itr, *itr2) += k_value;
	    K_mat(*itr2, *itr) += k_value;
	  }
	  count[*itr] = 0;
	}
	count_index.clear();	
	label_max++;
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    K_mat(i, i) /= 2;
  }

  return K_mat;
}



// compute a kernel value of a pair of graphs
double computeKernelValue(MatrixXi& e1, MatrixXi& e2, vector<int>& v1_label, vector<int>& v2_label, vector<double>& par, int kernel_type) {
  double Kval;
  switch (kernel_type) {
  case 1: // edge histogram kernel
    Kval = edgeHistogramKernel(e1, e2, -1.0);
    break;
  case 2: // vertex histogram kernel
    Kval = vertexHistogramKernel(v1_label, v2_label, -1.0);
    break;
  case 3: // vertex-edge histogram kernel
    Kval = vertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, -1.0);
    break;
  case 4: // vertex-vertex-edge histogram kernel
    Kval = vertexVertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, par[0]);
    break;
  case 5: // edge histogram kernel (Gaussian)
    Kval = edgeHistogramKernel(e1, e2, par[0]);
    break;
  case 6: // vertex histogram kernel (Gaussian)
    Kval = vertexHistogramKernel(v1_label, v2_label, par[0]);
    break;
  case 7: // vertex-edge histogram kernel (Gaussian)
    Kval = vertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, par[0]);
    break;
  
  case 8: // geometric random walk kernel
    Kval = geometricRandomWalkKernel(e1, e2, v1_label, v2_label, par[0]);
    break;
  case 9: // exponential random walk kernel
    Kval = exponentialRandomWalkKernel(e1, e2, v1_label, v2_label, par[0]);
    break;
  case 10: // k-step random walk kernel
    Kval = kstepRandomWalkKernel(e1, e2, v1_label, v2_label, par);
    break;
  default:
    Kval = 0;
    break;
  }
  return Kval;
}




MatrixXd CalculateKernelPy(vector<MatrixXi>& E, vector<vector<int> >& V_label, vector<int>& V_count, vector<int>& E_count, vector<int>& D_max, vector<double>& par, int kernel_type)
{  

  MatrixXd K(V_label.size(), V_label.size());
  

  vector<int> idx(V_label.size());
   iota(idx.begin(), idx.end(), 0);
    for (auto&& i : idx) {
      for (auto&& j : idx) {
          K(i, j) = computeKernelValue(E[i], E[j], V_label[i], V_label[j], par, kernel_type);
          K(j, i) = K(i, j);
      }

  }

  return K;

}






 

