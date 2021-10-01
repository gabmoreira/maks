/**
 * ravg.hpp
 *
 * 2021 Gabriel A. Moreira
 *
 * gmoreira at isr.tecnico.ulisboa.pt
 * https://github.com/gabmoreira/maks
 *
 * This software and the related documents  are provided as  is,  with no express
 * or implied  warranties,  other  than those  that are  expressly stated  in the
 * License.
 */

#ifndef RAVG_HPP
#define RAVG_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>


/**
 * Cycle Solver in SO(3).
 *
 * Closed-form cycle graph rotation averaging solver in SO(3).
 *
 * @param pairwise (input) Eigen::MatrixXd - 3 x 3n Eigen matrix with SO(3) blocks.
 * @param num_nodes (input) int - number of variables.
 * @param R (output) Eigen::MatrixXd - 3n x 3 Eigen matrix containing optimal solution.
 */
void cycleSolverSO3(const Eigen::Ref<Eigen::MatrixXd> pairwise,
                    int num_nodes,
                    Eigen::Ref<Eigen::MatrixXd> R);


/**
 * Rotation averaging Primal-Dual method in SO(3).
 *
 * Primal-dual update method for averaging rotations in SO(3).
 *
 * @param Rtilde (input) Eigen::SparseMatrix<double> - 3n x 3n sparse rotation adjacency matrix.
 * @param R (output) Eigen::MatrixXd - 3n x 3 solution.
 * @param A (input/output) Eigen::SparseMatrix<double> - 3n x 3n symmetric and sparse graph adjacency matrix.
 * @param num_nodes (input) Int - number of variables.
 * @param maxiter (input) int - maxiter.
 * @param dual (output) Double - dual problem.
 * @param eta (input) Double - minimum eigenvalue stopping criterion.
 * @param sigma (input) Double - spectral shift.
 */
void primalDualSO3(const Eigen::SparseMatrix<double>& Rtilde,
                   const Eigen::SparseMatrix<double>& A,
                   Eigen::Ref<Eigen::MatrixXd> R,
                   int num_nodes,
                   int maxiter,
                   double& dual,
                   double eta,
                   double sigma);


#endif /* RAVG_HPP */
