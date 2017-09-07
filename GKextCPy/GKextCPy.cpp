
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


// ========================================================= //
// ==================== Graphlet kernel ==================== //
// ========================================================= //
// ===== graphlet kernel for k = 4 ===== //
int find_min(int a, int b, int c) {
  int m;
  int mini = a;
  if (b < mini) mini = b;
  if (c < mini) mini = c;
  if (mini == a) {
    if (mini == b) {
      if (mini == c) {
  m = 7;
      } else {
  m = 4;
      }
    } else {
      if (mini == c) {
  m = 5;
      } else {
  m = 1;
      }
    }
  } else {
    if (mini == b) {
      if (mini == c) {
  m = 6;
      } else {
  m = 2;
      }
    } else {
      m = 3;
    }
  }
  return m;
}


void card_ThreeInter(vector<int>& L1, vector<int>& L2, vector<int>& L3, vector<int>& card) {
  card.resize(7);
  fill(card.begin(), card.end(), 0);
  int i = 0, j = 0, k = 0;

  while (i < (int)L1.size() && j < (int)L2.size() && k < (int)L3.size()) {
    int m = find_min(L1[i], L2[j], L3[k]);
    card[m - 1] += 1;
    switch(m) {
    case 1:
      i++; break;
    case 2:
      j++; break;
    case 3:
      k++; break;
    case 4:
      i++; j++; break;
    case 5:
      i++; k++; break;
    case 6:
      j++; k++; break;
    case 7:
      i++; j++; k++; break;
    }
  }

  if (i < (int)L1.size() || j < (int)L2.size() || k < (int)L3.size()) {
    if (i >= (int)L1.size() && j >= (int)L2.size()) {
      card[2] += (int)L3.size() - k;
      k = (int)L3.size();
    } else {
      if (i >= (int)L1.size() && k >= (int)L3.size()) {
  card[1] += (int)L2.size() - j;
  j = (int)L2.size();
      } else {
  if (j >= (int)L2.size() && k >= (int)L3.size()) {
    card[0] += (int)L1.size() - i;
    i = (int)L1.size();
  } else {
    if (i >= (int)L1.size()) {
      while (j < (int)L2.size() && k < (int)L3.size()) {
        if (L2[j] < L3[k]) {
    card[1]++;
    j++;
        } else {
    if (L2[j] > L3[k]) {
      card[2]++;
      k++;
    } else {
      card[5]++;
      j++;
      k++;
    }
        }
      }
    } else {
      if (j >= (int)L2.size()) {
        while (i < (int)L1.size() && k < (int)L3.size()) {
    if (L1[i] < L3[k]) {
      card[0]++;
      i++;
    } else {
      if (L1[i] > L3[k]) {
        card[2]++;
        k++;
      } else {
        card[4]++;
        i++;
        k++;
      }
    }
        }
      } else {
        if (k >= (int)L3.size()) {
    while (i < (int)L1.size() && j < (int)L2.size()) {
      if (L1[i] < L2[j]) {
        card[0]++;
        i++;
      } else {
        if (L1[i] > L2[j]) {
          card[1]++;
          j++;
        } else {
          card[3]++;
          i++;
          j++;
        }
      }
    }
        }
      }
    }
  }
      }
    }
  }
  if (i < (int)L1.size() || j < (int)L2.size() || k < (int)L3.size()) {
    if (i >= (int)L1.size() && j >= (int)L2.size()) {
      card[2] += (int)L3.size() - k;
    } else if (i >= (int)L1.size() && k >= (int)L3.size()) {
      card[1] += (int)L2.size() - j;
    } else if (j >= (int)L2.size() && k >= (int)L3.size()) {
      card[0] += (int)L1.size() - i;
    }
  }
}


void getIndices(vector<int>& o_set1, vector<int>& o_set2, vector<int>& inter, vector<int>& diff1, vector<int>& diff2) {
  vector<int> inter_(min(o_set1.size(), o_set2.size()), -1);
  vector<int> diff1_(max(o_set1.size(), o_set2.size()), -1);
  vector<int> diff2_(max(o_set1.size(), o_set2.size()), -1);

  int i = 0, j = 0;
  while (i < (int)o_set1.size() && j < (int)o_set2.size()) {
    if (o_set1[i] < o_set2[j]) {
      diff1_[i] = o_set1[i];
      i++;
    } else if (o_set1[i] > o_set2[j]) {
      diff2_[j] = o_set2[j];
      j++;
    } else {
      inter_[i] = o_set1[i];
      i++;
      j++;
    }
  }

  if (i < (int)o_set1.size()) {
    for (int k = i; k < (int)o_set1.size(); ++k) {
      diff1_[k] = o_set1[k];
    }
  } else if (j < (int)o_set2.size()) {
    for (int k = j; k < (int)o_set2.size(); ++k) {
      diff2_[k] = o_set2[k];
    }
  }

  inter.clear();
  for (auto&& x : inter_) {
    if (x >= 0) inter.push_back(x);
  }
  diff1.clear();
  for (auto&& x : diff1_) {
    if (x >= 0) diff1.push_back(x);
  }
  diff2.clear();
  for (auto&& x : diff2_) {
    if (x >= 0) diff2.push_back(x);
  }
}


//template<typename V>
VectorXd countGraphletsFour(vector<vector<int>>& al, VectorXd& count_gr) {

  double n = (double)al.size();
  vector<double> w = {1.0/12.0, 1.0/10.0, 1.0/8.0, 1.0/6.0, 1.0/8.0, 1.0/6.0, 1.0/6.0, 1.0/4.0, 1.0/4.0, 1.0/2.0, 0};
  vector<int> inter, diff1, diff2, card;
  vector<double> inter_count(11);
  vector<int> v;
  vector<int>::iterator it;

  double m = 0.0;
  for (auto&& vec : al) {
    m += (double)vec.size();
  }
  m /= 2.0;

  vector<int> v1(al.size());
  iota(v1.begin(), v1.end(), 0);
  for (auto&& i : v1) {
    for (auto&& j : al[i]) {
      double K = 0.0;
      fill(inter_count.begin(), inter_count.end(), 0.0);
      getIndices(al[i], al[j], inter, diff1, diff2);
      for (auto&& k : inter) {
  card_ThreeInter(al[i], al[j], al[k], card);
  inter_count[0] += 0.5 * (double)card[6];
  inter_count[1] += 0.5 * (double)(card[3] - 1.0);
  inter_count[1] += 0.5 * (double)(card[4] - 1.0);
  inter_count[1] += 0.5 * (double)(card[5] - 1.0);
  inter_count[2] += 0.5 * (double)card[0];
  inter_count[2] += 0.5 * (double)card[1];
  inter_count[2] += (double)card[2];
  inter_count[6] += n - (double)accumulate(card.begin(), card.end(), 0);
  K += 0.5 * (double)card[6] + 0.5 * (double)(card[4] - 1.0) + 0.5 * (double)(card[5] - 1.0) + card[2];
      }
      v.clear();
      v.resize(diff1.size());
      sort(diff1.begin(), diff1.end());
      sort(al[i].begin(), al[i].end());
      it = set_difference(diff1.begin(), diff1.end(), al[i].begin(), al[i].end(), v.begin());
      v.resize(it - v.begin());
      for (auto&& k : v) {
  card_ThreeInter(al[i], al[j], al[k], card);
  inter_count[1] += 0.5 * (double)card[6];
  inter_count[2] += 0.5 * (double)card[3];
  inter_count[2] += 0.5 * (double)card[4];
  inter_count[4] += 0.5 * (double)(card[5] - 1.0);
  inter_count[3] += 0.5 * (double)(card[0] - 2.0);
  inter_count[5] += 0.5 * (double)card[1];
  inter_count[5] += (double)card[2];
  inter_count[7] += n - (double)accumulate(card.begin(), card.end(), 0);
  K += 0.5 * (double)card[6] + 0.5 * (double)card[4] + 0.5 * (double)(card[5] - 1.0) + card[2];
      }
      v.clear();
      v.resize(diff2.size());
      sort(diff2.begin(), diff2.end());
      it = set_difference(diff2.begin(), diff2.end(), v1.begin(), v1.end(), v.begin());
      v.resize(it - v.begin());
      for (auto&& k : v) {
  card_ThreeInter(al[i], al[j], al[k], card);
  inter_count[1] += 0.5 * (double)card[6];
  inter_count[2] += 0.5 * (double)card[3];
  inter_count[4] += 0.5 * (double)(card[4] - 1.0);
  inter_count[2] += 0.5 * (double)card[5];
  inter_count[5] += 0.5 * (double)card[0];
  inter_count[3] += 0.5 * (double)(card[1] - 2.0);
  inter_count[5] += (double)card[2];
  inter_count[7] += n - (double)accumulate(card.begin(), card.end(), 0);
  K += 0.5 * (double)card[6] + 0.5 * (double)(card[4] - 1.0) + 0.5 * (double)card[5] + card[2];
      }
      inter_count[8] += m + 1.0 - (double)v1.size() - (double)al[i].size() - K;
      inter_count[9] += (n - (double)inter.size() - (double)diff1.size() - (double)diff2.size())
  * (n - (double)inter.size() - (double)diff1.size() - (double)diff2.size() - 1.0) / 2
  - (m + 1.0 - (double)v1.size() - (double)al[i].size() - K);

      for (int k = 0; k < (int)count_gr.size(); ++k) {
  count_gr(k) += inter_count[k] * w[k];
      }
    }
  }

  count_gr(10) = n * (n - 1.0) * (n - 2.0) * (n - 3.0) / (4.0 * 3.0 * 2.0) - count_gr.head(10).sum();
  return count_gr;
}



// ===== graphlet kernel for k = 3 ===== //
void getCardinality(vector<int>& o_set1, vector<int>& o_set2, vector<double>& card) {
  card.resize(3);
  fill(card.begin(), card.end(), 0.0);
  int i = 0, j = 0;
  while (i < (int)o_set1.size() && j < (int)o_set2.size()) {
    if (o_set1[i] < o_set2[j]) {
      card[0] += 1.0;
      i++;
    } else if (o_set1[i] > o_set2[j]) {
      card[1] += 1.0;
      j++;
    } else {
      i++;
      j++;
      card[2] += 1.0;
    }
  }
  card[0] += (double)((int)o_set1.size() - i);
  card[1] += (double)((int)o_set2.size() - j);
}

//template<typename V>

VectorXd countGraphletsThree(vector<vector<int>>& al, VectorXd& count_gr) {

  double n = (double)al.size();
  vector<double> w = {1.0/6.0, 1.0/4.0, 1.0/2.0};
  vector<double> card(3);

  vector<int> L1(al.size());
  iota(L1.begin(), L1.end(), 0);
  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      getCardinality(al[i], al[j], card);
      count_gr(0) += w[0] * card[2];
      count_gr(1) += w[1] * (card[0] + card[1] - 2.0);
      count_gr(2) += w[2] * (n - accumulate(card.begin(), card.end(), 0.0));
    }
  }
  count_gr(3) = n * (n - 1.0) * (n - 2.0) / 6.0 - (count_gr(0) + count_gr(1) + count_gr(2));
  return count_gr;
}


// ===== connected graphlet kernel for k = 5 ===== //
void getMinValue(MatrixXi& iam, vector<int>& idx, vector<int>& sums) {
  
  SparseMatrix<int> am;
  am = iam.sparseView();

  sums.clear();
  sums.resize(idx.size());
  fill(sums.begin(), sums.end(), 0);
  for (int i = 0; i < (int)idx.size(); ++i) {
    Int k = idx[i];
    for (SparseMatrix<int>::InnerIterator it(am, k); it; ++it) {
      if(find(idx.begin(), idx.end(), it.row()) != idx.end()) {
  sums[i] += it.value();
      }
    }
  }
  sums.push_back(1);
}


//template<typename V>

VectorXd countConnectedGraphletsFive(MatrixXi& am, vector<vector<int>>& al, VectorXd& count_gr) {
  
  int n = (int)al.size();

  vector<double> w = {1.0/120.0, 1.0/72.0, 1.0/48.0, 1.0/36.0, 1.0/28.0, 1.0/20.0, 1.0/14.0, 1.0/10.0, 1.0/12.0, 1.0/8.0, 1.0/8.0, 1.0/4.0, 1.0/2.0, 1.0/12.0, 1.0/12.0, 1.0/4.0, 1.0/4.0, 1.0/2.0, 0.0, 0.0, 0.0};
  //int n = (int)am.rows();
  vector<int> L1(n);
  iota(L1.begin(), L1.end(), 0);
  vector<int> idx(5);
  vector<int> sums;

  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
  if (k != i) {
    for (auto&& l : al[k]) {
      if (l != i && l != j) {
        for (auto&& m : al[l]) {
    if (m != i && m != j && m != k) {
      int aux = am.coeff(i, k) + am.coeff(i, l) + am.coeff(i, m) + am.coeff(j, l) + am.coeff(j, m) + am.coeff(k, m);
      if (aux == 6) {
        count_gr[0] += w[0];
      } else if (aux == 5) {
        count_gr[1] += w[1];
      } else if (aux == 4) {
        idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
        getMinValue(am, idx, sums);
        int aux1 = *min_element(sums.begin(), sums.end());
        if (aux1 == 2) {
          count_gr[3] += w[3];
        } else {
          count_gr[2] += w[2];
        }
      } else if (aux == 3) {
        idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
        getMinValue(am, idx, sums);
        sort(sums.begin(), sums.end());
        if (sums[0] == 1) {
          count_gr[8] += w[8];
        } else if (sums[1] == 3) {
          count_gr[4] += w[4];
        } else if (sums[2]== 2) {
          count_gr[13] += w[13];
        } else {
          count_gr[5] += w[5];
        }
      } else if (aux == 2) {
        idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
        getMinValue(am, idx, sums);
        vector<int> aux1;
        copy(sums.begin(), sums.end(), back_inserter(aux1));
        sort(aux1.begin(), aux1.end());
        if (aux1[0] == 1) {
          if (aux1[2] == 2) {
      count_gr[15] += w[15];
          } else {
      count_gr[9] += w[9];
          }
        } else {
          if (aux1[3] == 2) {
      count_gr[10] += w[10];
          } else {
      vector<int> ind;
      for (int ii = 0; ii < (int)sums.size(); ++ii) {
        if (sums[ii] == 3) ind.push_back(ii);
      }
      // idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
      if (am.coeff(idx[ind[0]], idx[ind[1]]) == 1) {
        count_gr[6] += w[6];
      } else {
        count_gr[14] += w[14];
      }
          }
        }
      } else if (aux == 1) {
        idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
        getMinValue(am, idx, sums);
        vector<int> aux1;
        copy(sums.begin(), sums.end(), back_inserter(aux1));
        sort(aux1.begin(), aux1.end());
        if (aux1[0] == 2) {
          count_gr[7] += w[7];
        } else if (aux1[1] == 1) {
          count_gr[17] += w[17];
        } else {
          vector<int> ind;
          for (int ii = 0; ii < (int)sums.size(); ++ii) {
      if (sums[ii] == 3) ind.push_back(ii);
          }
          for (int ii = 0; ii < (int)sums.size(); ++ii) {
      if (sums[ii] == 1) ind.push_back(ii);
          }
          if (am.coeff(idx[ind[0]], idx[ind[1]]) == 1) {
      count_gr[16] += w[16];
          } else {
      count_gr[11] += w[11];
          }
        }
      } else {
        count_gr[12] += w[12];
      }
    }
        }
      }
    }
  }
      }
    }
    // count graphlets of type 20
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
  if (k != i && am.coeff(i, k) == 0) {
    for (auto&& l : al[k]) {
      if (l != i && l != j && am.coeff(i, l) == 0 && am.coeff(j, l) == 0) {
        for (auto&& m : al[k]) {
    if (m != i && m != j && m != l && am.coeff(i, m) == 0 && am.coeff(j, m) == 0 && am.coeff(l, m) == 0) {
      count_gr[19] += w[19];
    }
        }
      }
    }
  }
      }
    }
    // count graphlets of type 19 and 21
    for (int j = 0; j < (int)al[i].size() - 3; ++j) {
      for (int k = j + 1; k < (int)al[i].size() - 2; ++k) {
  for (int l = k + 1; l < (int)al[i].size() - 1; ++l) {
    for (int m = l + 1; m < (int)al[i].size(); ++m) {
      int aux = am.coeff(al[i][j], al[i][k]) + am.coeff(al[i][j], al[i][l])
        + am.coeff(al[i][j], al[i][m]) + am.coeff(al[i][k], al[i][l])
        + am.coeff(al[i][k], al[i][m]) + am.coeff(al[i][l], al[i][m]);
      if (aux == 1) {
        count_gr[18]++;
      } else if (aux == 0) {
        count_gr[20]++;
      }
    }
  }
      }
    }
  }
  return count_gr;
}


//template<typename V>
// ===== connected graphlet kernel for k = 4 ===== //
VectorXd countConnectedGraphletsFour(MatrixXi& am, vector<vector<int>>& al, VectorXd& count_gr) {


  vector<double> w = {1.0/24.0, 1.0/12.0, 1.0/4.0, 0.0, 1.0/8.0, 1.0/2.0};
  int n = (int)am.rows();
  vector<int> L1(n);
  iota(L1.begin(), L1.end(), 0);

  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
  if (k != i) {
    for (auto&& l : al[k]) {
      if (l != i && l != j){
        int aux = am.coeff(i, k) + am.coeff(i, l) + am.coeff(j, l);
        if (aux == 3) {
    count_gr[0] += w[0];
        } else if (aux == 2) {
    count_gr[1] += w[1];
        } else if (aux == 1) {
    if (am.coeff(i, l) == 1) {
      count_gr[4] += w[4];
    } else {
      count_gr[2] += w[2];
    }
        } else {
    count_gr[5] += w[5];
        }
      }
    }
  }
      }
    }

    // count "stars"
    for (int j = 0; j < (int)al[i].size() - 2; ++j) {
      for (int k = j + 1; k < (int)al[i].size() - 1; ++k) {
  for (int l = k + 1; l < (int)al[i].size(); ++l) {
    if (am.coeff(al[i][j], al[i][k]) == 0 && am.coeff(al[i][j], al[i][l]) == 0 && am.coeff(al[i][k], al[i][l]) == 0) {
      count_gr[3]++;
    }
  }
      }
    }
  }
  return count_gr;
}


// ===== connected graphlet kernel for k = 3 ===== //
//template<typename V>
VectorXd countConnectedGraphletsThree(MatrixXi& am, vector<vector<int>>& al, VectorXd& count_gr) {

  vector<double> w = {1.0/2.0, 1.0/6.0};
  int n = (int)am.rows();
  vector<int> L1(n);
  iota(L1.begin(), L1.end(), 0);

  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
  if (k != i) {
    if (am.coeff(i, k) == 1) {
      count_gr[1] += w[1];
    } else {
      count_gr[0] += w[0];
    }
  }
      }
    }
  }
  return count_gr;
}


// Python function
MatrixXd CalculateGraphletKernelPy(vector<MatrixXi>& graph_adj_all, vector<vector<vector<int>>>& graph_adjlist_all, int k) {
  // decrement one to start indices from zero
  //for (auto&& X : graph_adjlist_all)
    //for (auto&& vec : X)
      //for (auto&& x : vec) x--;

  
  int freq_size;
  
  switch (k) {
    case 3: freq_size = 4; break;
    case 4: freq_size = 11; break;
  }
  
  //freq = MatrixXd::Zero(graph_adjlist_all.size(), freq_size);
  MatrixXd freq(graph_adjlist_all.size(), freq_size);

  vector<int> idx_graph(graph_adjlist_all.size());
  iota(idx_graph.begin(), idx_graph.end(), 0);

  VectorXd count_g;
  VectorXd freq_row;

  for (auto&& i : idx_graph) {
    
    freq_row = VectorXd::Zero(freq_size);

    //Eigen::Map<VectorXd>(freq_row.data(), freq.row(i).size()) = freq.row(i);
    
    if (k == 3) {
        count_g = countGraphletsThree(graph_adjlist_all[i], freq_row);
    } 
    else if (k == 4) {

        count_g = countGraphletsFour(graph_adjlist_all[i], freq_row);
    } 

    freq.row(i) = count_g;

    if (freq.row(i).sum() != 0) {
      freq.row(i) /= freq.row(i).sum();
    }
  }

  MatrixXd K = freq * freq.transpose();

  return K;
}

// Python function
MatrixXd CalculateConnectedGraphletKernelPy(vector<MatrixXi>& graph_adj_all, vector<vector<vector<int>>>& graph_adjlist_all, int k) {
  // decrement one to start indices from zero
  //for (auto&& X : graph_adjlist_all)
    //for (auto&& vec : X)
      //for (auto&& x : vec) x--;

  //MatrixXd freq;
  int freq_size;
  
  switch (k) {
    case 3: freq_size = 2; break;
    case 4: freq_size = 6; break;
    case 5: freq_size = 21; break;
  }

  MatrixXd freq(graph_adjlist_all.size(), freq_size);
  //freq = MatrixXd::Zero(graph_adjlist_all.size(), freq_size);

  vector<int> idx_graph(graph_adjlist_all.size());
  iota(idx_graph.begin(), idx_graph.end(), 0);

  VectorXd count_g;
  VectorXd freq_row;

  for (auto&& i : idx_graph) {

    //VectorXd count_g;
    //VectorXd freq_row;

    //Eigen::Map<VectorXd>(freq_row.data(), freq.row(i).size()) = freq.row(i);
    freq_row = VectorXd::Zero(freq_size);

    if (k == 3) {        
        count_g = countConnectedGraphletsThree(graph_adj_all[i], graph_adjlist_all[i], freq_row);

    } 
    else if (k == 4) {
        count_g = countConnectedGraphletsFour(graph_adj_all[i], graph_adjlist_all[i], freq_row);

    } 
    else if (k == 5) {
        count_g = countConnectedGraphletsFive(graph_adj_all[i], graph_adjlist_all[i], freq_row);
    }

    //freq.row(i) = Eigen::Map<VectorXd>(count_g.data(), count_g.size());
    freq.row(i) = count_g;
    
    if (freq.row(i).sum() != 0) {
      freq.row(i) /= freq.row(i).sum();
    }
  }
  MatrixXd K = freq * freq.transpose();
  
  return K;
}

 

