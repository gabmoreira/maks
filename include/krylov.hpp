/*******************************************************************************
* krylov.hpp
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
#ifndef KRYLOV_H
#define KRYLOV_H

#define EIGEN_USE_MKL_ALL

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

namespace kry {

/*******************************************************************************
 * Solves Sparse linear system A*x = b based on LDL^t decomposition
 *
 * @param A Input sparse matrix
 * @param b Input vector (right-hand side of the equation)
 * @param x Output vector
 * @return Nothing
 ******************************************************************************/
void sparseSolve(Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>* solver,
                const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
                Eigen::VectorXd& x);


/*******************************************************************************
 * For a real sparse symmetric matrix A, this function computes k real eigenvalues
 * near a real target sigma using the Krylov-Schur method (G. W. Stewart 2001).
 * Behaves very similarly to MATLAB's eigs in terms of precision. The matrix must
 * have a large degree of sparsity for this method to be performant. Make sure
 * it is built via triplets and not by casting a dense matrix to sparse.
 *
 * @param A Input sparse matrix of doubles used to compute eigenvalues and
 *          eigenvectors.
 * @param eigenvalues Output eigen vector of size k used to store the computed
 *                    eigenvalues.
 * @param eigenvectors Output eigen matrix size n x k. The normalized
 *                     eigenvectors are stored in its columns.
 * @param k Input number of eigenvalues and eigenvectors to be computed. Must be in
 *          accordance with the dimensions of the matrices eigenvalues and
 *          eigenvectors.
 * @param sigma Input eigenvalue real target
 * @return Nothing
 ******************************************************************************/
void eigs(Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>* solver,
          const Eigen::SparseMatrix<double>& A, Eigen::VectorXd& eigenvalues,
          Eigen::MatrixXd& eigenvectors, unsigned int k, double sigma);
};
#endif
