/*******************************************************************************
* krylov.cpp
*
* 2020 Gabriel A. Moreira
*
* g.antunes.moreira at gmail dot com
* https://github.com/gabmoreira
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/* Comment line below if system does not have MKL installed */
#define EIGEN_USE_MKL_ALL

#include <vector>
#include <numeric>
#include <iostream>

#include <Eigen/Eigenvalues>

#include "internal.hpp"
#include "krylov.hpp"

using std::min;
using std::max;
using std::iota;
using std::fill;
using std::sort;
using std::vector;
using std::accumulate;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Success;

typedef Eigen::SparseMatrix<double> Sparse;
typedef Eigen::PardisoLDLT<Sparse> MKL_LDLT;

/* Exceptions thrown by this module */
#define ERROR_STARTING_VECTOR 0
#define ERROR_EIGEN_CHOLESKY 1
#define ERROR_REORTHOGONALIZATION 2

#define MAXITER 300
#define TOL 1e-11
#define EPS 2.2204e-12
#define RES_TOL 3.6668e-10
#define KRYLOV_SUBSPACE_SIZE 20

namespace kry {

/*******************************************************************************
 * Solves sparse system via LDL^t decomposition
 ******************************************************************************/
void sparseSolve(MKL_LDLT* solver, const Sparse& A, const VectorXd& b, VectorXd& x) {
  solver->factorize(A);
  x = solver->solve(b);
};


/*******************************************************************************
 * Computes k eigenvectors and eigenvalues of a real and symmetric matrix
 * using Krylov-Schur method by Stewart 2001.
 ******************************************************************************/
void eigs(MKL_LDLT* solver, const Sparse& A, VectorXd& eigenvalues,
          MatrixXd& eigenvectors, unsigned int k, double sigma) {
  unsigned int n = A.rows();

  /* Define our sparse linear solver and compute the LU decomposition.
  The eigen equation is (A-sigma*I)^-1 v = lambda v. We need to compute the
  iterative powers of ((A-sigma*I)^-1)^k v. In order to avoid explicitly
  computing the inverse we compute at each stage v_{k+1} = solve(A-sigma*I, v_{k}).
  This solver is equivalent to the popular UMFPACK used by MATLAB and
  computes results to the same degree of machine precision. */
  Sparse speye(n,n);
  speye.setIdentity();
  Sparse AMinusSigmaI = A;
  AMinusSigmaI -= sigma * speye;
  solver->factorize(AMinusSigmaI);

  /* Check that the Cholesky decomposition was successful */
  if(solver->info()!=Success)
    throw ERROR_EIGEN_CHOLESKY;

  MatrixXd V = MatrixXd::Zero(n, KRYLOV_SUBSPACE_SIZE);                         // Matrix to hold the orthogonalized Krylov vectors
  unsigned int sizeV = 0;                                                       // Size of the V matrix (not actual size but size in use)

  int stop = 0;                                                                 // Stop the algorithm flag (reorthogonalization gone bad)
  unsigned int k0 = k;                                                          // Number of selected eigenvalues (may vary within the algorithm)

  vector<double> Alpha;                                                         // Store projection data
  vector<double> Beta;                                                          // Store projection data
  Alpha.reserve(KRYLOV_SUBSPACE_SIZE);
  Beta.reserve(KRYLOV_SUBSPACE_SIZE);

  MatrixXd c;

  double normRes = 0;                                                           // Stores residual norm
  bool justRestarted = false;

  vector<int> converged(k0);                                                    // Flags indicating converged eigenvalues
  unsigned int nconv = 0;                                                       // Number of converged eigenvalues

  MatrixXd ritz_vectors;
  MatrixXd ritz_values;

  vector<int> idx_vec(KRYLOV_SUBSPACE_SIZE);                                    // Used to compute argsorts

  MatrixXd H = MatrixXd::Zero(KRYLOV_SUBSPACE_SIZE, KRYLOV_SUBSPACE_SIZE);      // Hessenberg matrix

  /* Starting vector for the Krylov iterations */
  VectorXd v0 = VectorXd::Random(n,1);
  v0.normalize();

  VectorXd v = solver->solve(v0);                                               // One step of the power iteration
  v.normalize();

  /* Main loop of the Krylov-Schur method */
  for (unsigned int i = 0; i < MAXITER; ++i) {
    /* Loop to build the invariant subspace in V, H (Hessenberg) */
    for (unsigned int j = sizeV; j < KRYLOV_SUBSPACE_SIZE; ++j) {
      V.col(j) = v;                                                             // Store the Krylov vector v from the previous iteration in the Krylov matrix V
      VectorXd r = solver->solve(v);                                            // Compute Krylov power using Cholesky solver: r = (A-sigma I)^-1 * v
      double alpha = v.transpose() * r;                                         // Component of the new vector r in the direction of the previous vector v

      /* Gram-Schmidt orthogonalization */
      if (j == 0) {
        r.noalias() -= alpha * v;                                               // 1st iteration: just subtract the projection on the previous vector
      } else if (justRestarted) {
        r.noalias() -= internal::projectOnto(V, r, j+1);                        // Subtract projection on all the previous Krylov vectors of V
        justRestarted = false;
      } else {
        r.noalias() -= alpha * v;                                               // Subtract projection on the previous vector and residuals on one before that
        r.noalias() -= normRes * V.col(j-1);
      };

      /* Robust reorthogonalization of r against the columns of V. This helps
       * to prevent numerical errors (implemented in  MATLAB as  well)  */
      internal::reorthogonalize(V, r, normRes, j+1, stop);

      /* Check if it was possible to reorthogonalize. If not, we cannot find the
       * invariant subspace and must exit the program */
      if (stop)
        throw ERROR_REORTHOGONALIZATION;

      Alpha.push_back(alpha);                                                   // Save projection of r on the the previous vector v (alpha)
      Beta.push_back(normRes);                                                  // Save projection of r on the previous previous vector (beta)
      v = r;                                                                    // Set the current vector v as the orthogonalized vector r
    };
    /* End of subspace loop */

    /* Build Hessenberg matrix H */
    H = MatrixXd::Zero(KRYLOV_SUBSPACE_SIZE, KRYLOV_SUBSPACE_SIZE);             // Hessenberg matrix

    MatrixXd H1 = MatrixXd::Zero(Alpha.size(),Alpha.size());

    Eigen::Map<VectorXd> Alpha_vec(Alpha.data(), Alpha.size());
    Eigen::Map<VectorXd> Beta_vec(Beta.data(), Beta.size()-1);

    H1.diagonal(-1).array() += Beta_vec.array();
    H1.diagonal(0).array() += Alpha_vec.array();
    H1.diagonal(1).array() += Beta_vec.array();

    H.bottomRightCorner(H1.rows(), H1.cols()) += H1;

    if (ritz_values.rows() > 0)
      H.topLeftCorner(ritz_values.rows(),ritz_values.rows()) += ritz_values.asDiagonal();

    if (ritz_values.rows() > 0) {
      for (unsigned int l = 0; l < k; ++l) {
        H(l,k+1) = c(0,l);
        H(k+1,l) = c(0,l);
      };
    };
    /* Hessenberg matrix complete  */

    Alpha.resize(0);
    Beta.resize(0);

    // Compute Ritz pairs (Eigenpairs of the Hessenberg matrix)
    Eigen::EigenSolver<MatrixXd> eig(H);
    ritz_vectors = eig.eigenvectors().real();
    ritz_values = eig.eigenvalues().real();

    /* Compute residuals to find out what has converged */
    vector<double> res(ritz_vectors.cols());
    for (unsigned int l = 0; l < ritz_vectors.cols(); ++l) { res[l] = abs(normRes * ritz_vectors(ritz_vectors.rows()-1, l)); };

    /* Argsort eigenvalues and residuals in descending order (based on eigenvalues) */
    iota(idx_vec.begin(), idx_vec.end(), 0);
    sort(idx_vec.begin(), idx_vec.end(), [&](size_t a, size_t b) { return abs(ritz_values(a)) > abs(ritz_values(b)); });

    /* Check which ones have converged */
    for (unsigned int l = 0; l < k0; ++l) {
       converged[l] = (int) ( res[idx_vec[l]] < TOL * max(RES_TOL, abs(ritz_values(idx_vec[l]))) );
    };

    nconv = accumulate(converged.begin(), converged.end(), 0);                  // Number of eigenvalues which have converged

    if (nconv >= k0) {
      break;                                                                    // More than k0 eigenvalues converged: We're good to go!
    } else {
      k = k0 + min( (float) nconv, (float) floor(((float) (KRYLOV_SUBSPACE_SIZE-k0)) / 2.0f) );     // Adjust k to prevent stagnation
      if (k == 1)
        k = floor(KRYLOV_SUBSPACE_SIZE / 2.0f);
    };

    MatrixXd _ritz_vectors = MatrixXd::Zero(KRYLOV_SUBSPACE_SIZE, k);
    VectorXd _ritz_values = VectorXd::Zero(KRYLOV_SUBSPACE_SIZE, 1);

    for (unsigned int i = 0; i < k; ++i) {
      _ritz_vectors.col(i) += ritz_vectors.col(idx_vec[i]);
      _ritz_values(i) += ritz_values(idx_vec[i]);
    };

    ritz_vectors = _ritz_vectors;
    ritz_values = _ritz_values;

    V.leftCols(k) = V * ritz_vectors;                                           // Store variables for the next iteration
    c = normRes * ritz_vectors.bottomRows<1>();

    justRestarted = true;
    sizeV = k;
  };   /* End of the iteration loop */

  /* More reasonable people would check for convergence here and return only
   * eigenpairs which converged. We, on the other hand, do not have time for
   * such nonsense */
  eigenvectors = MatrixXd::Zero(ritz_vectors.rows(), k0);
  eigenvalues = VectorXd::Zero(k0);

  for (unsigned int l = 0; l < k0; ++l) {
    eigenvalues(l) += ritz_values(idx_vec[l]);
    eigenvectors.col(l) += ritz_vectors.col(idx_vec[l]);
  };

  c = normRes * eigenvectors.bottomRows<1>();

  eigenvectors = V * eigenvectors;
  eigenvectors = eigenvectors.array().rowwise() * eigenvalues.transpose().array();
  eigenvectors.noalias() += v * c;

  /* Eigenvector normalization */
  VectorXd norms = eigenvectors.colwise().norm();
  eigenvectors.array().rowwise() /= norms.array().transpose();

  // Eigenvalue shift and invert
  for (unsigned int l = 0; l < eigenvalues.rows(); ++l) { eigenvalues(l) = 1.0f / eigenvalues(l) + sigma; };
}; /* End of Krylov-Schur for symmetric matrices*/
};
