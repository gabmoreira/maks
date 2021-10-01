/**
 * ravg.cpp
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

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>

#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "err.hpp"
#include "ravg.hpp"

using std::vector;
using std::string;

using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::Ref;

using namespace std::chrono;

typedef Eigen::Triplet<double>      Triplet;
typedef Eigen::SparseMatrix<double> Sparse;
typedef Eigen::PardisoLDLT<Sparse>  PardisoLDLT;


/**
 * Cycle graph solver in SO(3).
 *
 * Closed-form cycle graph rotation averaging solver in SO(3).
 *
 * @param edge_r (input) Eigen::MatrixXd - 3 x 3n Eigen matrix with SO(3) blocks. Assumes blocks are stacked as [ R_{1,2} R_{2,3} ... R_{n,1} ]
 * @param num_nodes (input) int - number of variables.
 * @param R (output) Eigen::MatrixXd - 3n x 3 Eigen matrix containing optimal solution.
 */
void cycleSolverSO3(const Ref<MatrixXd> edge_r, int num_nodes, Ref<MatrixXd> R) {
    auto start = high_resolution_clock::now();
    
    // Cycle error E = R_12 * R_23 * R_n-1,n * R_n,1
    MatrixXd E = Matrix3d::Identity();
    
    for (int i = 0; i < num_nodes; ++i)
        E = E * edge_r.block<3,3>(0, i*3);
    
    // Diagonal blocks of change-of-basis matrix U
    vector<Matrix3d> U(num_nodes);
    
    // Build first block
    Eigen::SelfAdjointEigenSolver<MatrixXd> eig(E + E.transpose());
    U[0] = eig.eigenvectors();
    
    // Build other blocks
    for (int i = 1; i < num_nodes; ++i)
        U[i] = edge_r.block<3,3>(0, (i-1)*4).transpose() * U[i-1];
    
    /* Change basis according to U^\top * Rtilde * U
     * We need only compute the corner block U_n * Rtilde_n1 * U_1 */
    MatrixXd one_param_edge_r = MatrixXd::Zero(3,3);
    one_param_edge_r += U[num_nodes-1].transpose() * edge_r.block<3,3>(0, (num_nodes-1)*3) * U[0];
    
    // Gamma is the angular cycle error
    double gamma_cos = 0.5 * (one_param_edge_r(0,0) + one_param_edge_r(1,1));
    double sine_sign = (double) ((one_param_edge_r(1,0) > 0) - (one_param_edge_r(1,0) < 0));
    double gamma     = sine_sign * acos(gamma_cos);
    
    double gamma_mean = gamma / num_nodes;
    double theoretical_min = - ((double) num_nodes) * (5.0f + 4.0f * cos(gamma_mean));
    
    // Build solution in SO(2) and change basis to produce solution in SO(3)
    R.block<3,3>(0,0) += U[0];
    VectorXd theta = VectorXd::Zero(num_nodes,1);
    
    for (int i = 1; i < num_nodes; ++i) {
        theta(i) += theta(i-1) + gamma_mean;
        double c = cos(theta(i));
        double s = sin(theta(i));
        
        R(i*3,0)   += c;
        R(i*3,1)   -= s;
        R(i*3+1,0) += s;
        R(i*3+1,1) += c;
        R(i*3+2,2) += 1;
        
        // Revert the basis-change
        R.block<3,3>(i*3,0) = U[i] * R.block<3,3>(i*3,0);
    };
    
    // Fix gauge freedom by setting R_1 = I_3
    R = R * R.block<3,3>(0,0).transpose();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    
    // Verbose
    printf("Cycle error (mean): %1.8f\n", gamma_mean);
    printf("Theoretical minimum: %1.8f\n", theoretical_min);
    printf("CPU time: %.6f\n", static_cast<long long int>(duration.count())*1e-6);
};


/**
 * Primal-Dual in SO(3)
 *
 * Primal-dual update method for averaging rotations in SO(3).
 * Finds R in SO(3)^n that minimizes the trace of R^\top Rtilde R
 *
 * @param Rtilde (input) Eigen::SparseMatrix<double> - 3n x 3n sparse rotation adjacency matrix.
 * @param A (input/output) Eigen::SparseMatrix<double> - n x n symmetric graph adjacency matrix.
 * @param R (output) Eigen::MatrixXd - 3n x 3 solution.
 * @param num_nodes (input) int - number of variables.
 * @param maxiter (input) int - maxiter.
 * @param dual (output) double - value of the dual problem.
 * @param eta (input) double - minimum eigenvalue stopping criterion.
 * @param sigma (input) double - spectral shift.
 */
void primalDualSO3(const Sparse& Rtilde,
                   const Sparse& A,
                   Ref<MatrixXd> R,
                   int num_nodes,
                   int maxiter,
                   double& dual,
                   double eta,
                   double sigma) {
    
    printf("\n%*s SO(3) Primal-dual\n", 20, " ");
    
    auto start = high_resolution_clock::now();
    
    // Graph degree vector
    VectorXd D = A * VectorXd::Ones(num_nodes);
            
    /* Initialize symmetric block diagonal Lagrange multiplier.
     * Should initialize all 3x3 blocks instead of just the main
     * diagonal to ensure the correct NNZ pattern. */
    vector<Triplet> lambda_buffer(num_nodes * 9);
    
    unsigned int k = 0;
    unsigned int j = 0;
    for (unsigned int i = 0; i < num_nodes; ++i) {
        lambda_buffer[j]   = Triplet( k,   k,   D(i) );
        lambda_buffer[j+1] = Triplet( k+1, k+1, D(i) );
        lambda_buffer[j+2] = Triplet( k+2, k+2, D(i) );
        lambda_buffer[j+3] = Triplet( k,   k+1, 1e-15 );
        lambda_buffer[j+4] = Triplet( k,   k+2, 1e-15 );
        lambda_buffer[j+5] = Triplet( k+1, k,   1e-15 );
        lambda_buffer[j+6] = Triplet( k+1, k+2, 1e-15 );
        lambda_buffer[j+7] = Triplet( k+2, k,   1e-15 );
        lambda_buffer[j+8] = Triplet( k+2, k+1, 1e-15 );
        j += 9;
        k += 3;
    };
    
    Sparse Lambda(num_nodes * 3, num_nodes * 3);
    Lambda.setFromTriplets(lambda_buffer.begin(), lambda_buffer.end());
        
    Sparse L = Lambda - Rtilde;
    L.makeCompressed();
    
    // LDLt solver used in the shift-and-invert spectral transform
    PardisoLDLT ldlt_solver;
    
    // Analyze the NNZ pattern
    ldlt_solver.analyzePattern(L);

    // Store eigenvalues
    VectorXd evals = VectorXd::Zero(3, 1);
        
    // Main primal-dual loop
    for (unsigned int i = 0; i < maxiter; ++i) {
        /* Call symmetric Krylov-Schur LDLt eigensolver on the matrix
         * Lambda - Rtilde. The three eigenvectors will be stored in R. */
        alg::symKrylovSchurLDL(&ldlt_solver, L, evals, R, 3, sigma);
        
        // Correct gauge freedom
        Matrix3d gauge = R.block<3,3>(0,0);
        R = R * gauge.inverse();
        
        // SO(3) projection
        for (unsigned int j = 0; j < 3*num_nodes; j += 3)
            alg::orthoProcrustesSO3(R.block<3,3>(j, 0));
        
        // Build Lambda from triplets for the next iteration
        MatrixXd RtildeR = Rtilde * R;
        dual = 3 * num_nodes;
        for (unsigned int k = 0; k < num_nodes*3; k+=3) {
            Matrix3d Rk = R.block<3,3>(k,0);
            Rk.transposeInPlace();
            MatrixXd block = RtildeR.block<3,3>(k,0) * Rk;
            
            // Dual <- + trace(Lambda_i)
            dual += block(0,0) + block(1,1) + block(2,2);
            
            L.coeffRef(k,k)     = block(0,0);
            L.coeffRef(k+1,k+1) = block(1,1);
            L.coeffRef(k+2,k+2) = block(2,2);
            L.coeffRef(k,k+1)   = 0.5 * (block(0,1) + block(1,0));
            L.coeffRef(k+1,k)   = L.coeffRef(k,k+1);
            L.coeffRef(k,k+2)   = 0.5 * (block(0,2) + block(2,0));
            L.coeffRef(k+2,k)   = L.coeffRef(k,k+2);
            L.coeffRef(k+1,k+2) = 0.5 * (block(1,2) + block(2,1));
            L.coeffRef(k+2,k+1) = L.coeffRef(k+1,k+2);
        };
        
        // Verbose
        printf("Iteration %d   Dual %.9f   Min eigenvalue %1.3e\n", i, dual, abs(evals.minCoeff()));

        // Stopping criterion
        if (abs(evals.minCoeff()) < eta)
            break;
    };
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    printf("\nCPU time : %.6f\n", static_cast<long long int>(duration.count())*1e-6);
};
