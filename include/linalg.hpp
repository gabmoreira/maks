/**
 * linalg.hpp
 *
 * 2021 Gabriel Moreira
 *
 * https://github.com/gabmoreira/maks
 *
 * This software and the related documents  are provided as  is,  with no express
 * or implied  warranties,  other  than those  that are  expressly stated  in the
 * License.
 *
 * Copyright © 2021 Gabriel Moreira. All rights reserved.
 */

#ifndef LINALG_HPP
#define LINALG_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

/* Linear algebra namespace */
namespace alg {


/**
 * Converts two vectors i and j to matrix indices (in-place).
 *
 * @param ei (input) Eigen::VectorXi - row indices.
 * @param ej (input) Eigen::VectorXi - column indices.
 * @param num_nodes (output) int - number of nodes.
 */
void convertToIdx(Eigen::Ref<Eigen::VectorXi> ei,
                  Eigen::Ref<Eigen::VectorXi> ej,
                  int& num_nodes);


/**
 * Builds sparse graph adjacency matrix.
 *
 * @param ei (input) Eigen::VectorXi - row indices.
 * @param ej (input) Eigen::VectorXi - column indices.
 * @param num_edges (input) int - number of edges.
 * @param sym (input) bool - whether to make the matrix symmetric.
 * @param dst (output) Eigen::SparseMatrix<double> - output adjacency matrix.
 */
void adjacency(const Eigen::Ref<Eigen::VectorXi> ei,
               const Eigen::Ref<Eigen::VectorXi> ej,
               int   num_edges,
               bool  sym,
               Eigen::SparseMatrix<double>& dst);


/**
 * Builds sparse block matrix.
 *
 * @param ei (input) Eigen::VectorXi - row indices.
 * @param ej (input) Eigen::VectorXi - column indices.
 * @param blocks (input) Eigen::MatrixXd - row block matrix containing contiguous blocks.
 * @param block_height (input) int - block height.
 * @param block_width (input) int - block width.
 * @param num_blocks (input) int - number of contiguous blocks.
 * @param sym (input) bool - whether to make the matrix symmetric.
 * @param dst (output) Eigen::SparseMatrix<double> - output adjacency matrix.
 */
void sparseBlocks(const Eigen::Ref<Eigen::VectorXi> ei,
                  const Eigen::Ref<Eigen::VectorXi> ej,
                  const Eigen::Ref<Eigen::MatrixXd> blocks,
                  int   block_height,
                  int   block_width,
                  int   num_blocks,
                  bool  sym,
                  Eigen::SparseMatrix<double>& dst);


/**
 * Solves the orthogonal Procrustes problem in SO(3) (in-place).
 *
 * @param mat (input/output) Eigen::MatrixXd - 3x3 matrix
 */
void orthoProcrustesSO3(Eigen::Ref<Eigen::Matrix3d> mat);


/**
 * Solves the orthogonal Procrustes problem in O(3) (in-place).
 *
 * @param mat (input/output) Eigen::MatrixXd - 3x3 matrix.
 */
void orthoProcrustesO3(Eigen::Ref<Eigen::Matrix3d> mat);


/**
 * Finds the projection of a vector on a subpace spanned by the columns of a matrix.
 *
 * @param subspace_basis (input) Eigen::MatrixXd - Dense matrix whose columns span a subspace
 *        where we want to project the vector v.
 * @param idx (input) unsigned int - Sub-selects of the first idx columns of A only.
 * @param v (input) Eigen::VectorXd - Vector to project.
 * @return The component of the vector in the subspace spanned by the first
 *         idx columns of the input matrix A.
 */
Eigen::VectorXd projectionOnto(const Eigen::Ref<Eigen::MatrixXd> subspace_basis,
                               unsigned int idx,
                               const Eigen::Ref<Eigen::VectorXd> v);


/**
 * Projects a vector onto a subpace spanned by the columns of a matrix (in-place)
 *
 * @param subspace_basis (input) Eigen::MatrixXd - Dense matrix whose columns span a subspace where we
 *        want to project the vector v.
 * @param idx (input) unsigned int - Sub-selects of the first idx columns of A only.
 * @param v (input/output) Eigen::VectorXd - Vector to project.
 */
void projectOnto(const Eigen::Ref<Eigen::MatrixXd> subspace_basis,
                 unsigned int idx,
                 Eigen::Ref<Eigen::VectorXd> v);


/**
 * Symmetric sparse solver via LDLT decomposition.
 *
 * Solves sparse system Ax = b via LDLt decomposition.
 *
 * @param solver (input) Eigen::PardisoLDLT<Sparse> - LDLt Intel PARDISO solver
 * @param mat (input) Eigen::MatrixXd - matrix to factorize.
 * @param rhs (input) Eigen::VectorXd - right hand side.
 * @param x (output) Eigen::VectorXd - solution of the system.
 */
void sparseSolverLDL(Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>* solver,
                     const Eigen::SparseMatrix<double>& mat,
                     const Eigen::Ref<Eigen::MatrixXd> rhs,
                     Eigen::Ref<Eigen::MatrixXd> x);


/**
 * Robust reorthogonalization. For an input matrix V and an input
 * vector r, tries to iteratively orthogonalize r against the j first columns
 * of V. This prevents numerical errors and ensures orthogonality. If, by any
 * chance, and after multiple iterations of subtracting projections, the vector
 * has a norm smaller than 1/sqrt(2) of its starting norm, it cannot be
 * reorthogonalized. This method then attempts to find another orthogonal vector
 * through random restarts. The flag stop is activated if everything fails.
 * (See G. W. Stewart 2001)
 *
 * @param V (input) Dense matrix of doubles whose columns contain orthogonal basis
 *        vectors
 * @param r (input/output) Vector to orthogonalize
 * @param residual_norm (output) Residual norm or r after it has been orthogonalized against the first j columns of V.
 * @param idx (input) Index specifying the number of columns of V used in the orthogonalization.
 * @param stop (output) Stop the algorithm flag. Indicates that it is impossible to find, to machine precision, a vector r orthogonal to the first j cols of V.
 */
void reorthogonalize(const Eigen::Ref<Eigen::MatrixXd> V,
                     Eigen::Ref<Eigen::VectorXd> r,
                     double& residual_norm,
                     unsigned int idx,
                     int& stop);


/**
 * Symmetric Krylov-Schur LDLt eigensolver.
 *
 * For a real sparse symmetric matrix A, this function computes k real eigenvalues
 * near a real target sigma using the Krylov-Schur method (G. W. Stewart 2001).
 *
 * @param solver (input) Eigen::PardisoLDLT<Sparse> - LDLt Intel PARDISO solver.
 * @param mat (input) Eigen::SparseMatrix<double> - sparse matrix to compute eigenvalues and eigenvectors.
 * @param eigenvalues (output) Eigen::VectorXd - k x 1 matrix to store the computed eigenvalues.
 * @param eigenvectors (output) Eigen::MatrixXd - n x k matrix to store normalized eigenvectors.
 * @param k (input) int - number of eigenvalues and eigenvectors to be computed.
 * @param sigma (input) double - spectral shift / eigenvalue target.
*/
void symKrylovSchurLDL(Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>* solver,
                       const Eigen::SparseMatrix<double>& mat,
                       Eigen::Ref<Eigen::VectorXd> eigenvalues,
                       Eigen::Ref<Eigen::MatrixXd> eigenvectors,
                       unsigned int k,
                       double sigma);

};

#endif /* LINALG_HPP */
