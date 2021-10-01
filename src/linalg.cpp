/**
 * linalg.cpp
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

#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <Eigen/Eigenvalues>

#include "err.hpp"
#include "linalg.hpp"

using std::min;
using std::max;
using std::iota;
using std::fill;
using std::sort;
using std::vector;
using std::string;
using std::accumulate;

using Eigen::Ref;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::VectorXi;

typedef Eigen::Triplet<double>      Triplet;
typedef Eigen::SparseMatrix<double> Sparse;
typedef Eigen::PardisoLDLT<Sparse>  PardisoLDLT;

/* Parameters for the Krylov-Schur eigensolver */
#define KRYLOV_MAXITER           300
#define KRYLOV_TOL               1e-15
#define KRYLOV_RES_TOL           2.22e-16
#define KRYLOV_SUBSPACE_MAX_SIZE 20

namespace alg {


/**
 * Converts two vectors i and j to matrix indices (in-place).
 *
 * @param ei (input) Eigen::VectorXi - row indices.
 * @param ej (input) Eigen::VectorXi - column indices.
 * @param num_nodes (output) int - number of nodes.
 */
void convertToIdx(Ref<VectorXi> ei, Ref<VectorXi> ej, int& num_nodes) {
    
    int num_edges = (int) ei.rows();
    
    // Obtain node list from edges
    vector<int> nodes;
    
    nodes.insert(nodes.end(), ei.data(), ei.data() + num_edges);
    nodes.insert(nodes.end(), ej.data(), ej.data() + num_edges);
        
    sort(nodes.begin(), nodes.end());
    auto last = unique(nodes.begin(), nodes.end());
    
    num_nodes = (int) (last - nodes.begin());
    
    // Build edge map
    vector<int> idx_map(num_nodes);
    
    for (int i = 0; i < num_nodes ; ++i)
        idx_map[nodes[i]] = i;
    
    // Store matrix indices
    for (int i = 0; i < num_edges; ++i) {
        ei(i) = idx_map[ei(i)];
        ej(i) = idx_map[ej(i)];
    };
};


/**
 * Builds sparse graph adjacency matrix.
 *
 * @param ei (input) Eigen::VectorXi - row indices.
 * @param ej (input) Eigen::VectorXi - column indices.
 * @param num_edges (input) int - number of edges.
 * @param sym (input) bool - whether to make the matrix symmetric.
 * @param dst (output) Eigen::SparseMatrix<double> - output adjacency matrix.
 */
void adjacency(const Ref<VectorXi> ei, const Ref<VectorXi> ej, int num_edges, bool sym, Sparse& dst) {
    
    int buffer_size = num_edges;
    
    if (sym)
        buffer_size *= 2;
    
    // Allocate buffer
    vector<Triplet> buffer(buffer_size);
    
    // Fill in the data
    for (int i = 0; i < num_edges; ++i)
        buffer[i] = Triplet( ei(i), ej(i), 1.0f );

    // Duplicate entries if symmetric flag is true
    if (sym)
        for (int i = 0; i < num_edges; ++i)
            buffer[num_edges + i] = Triplet( ej(i), ei(i),  1.0f );

    dst.setFromTriplets(buffer.begin(), buffer.end());
};


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
void sparseBlocks(const Ref<VectorXi> ei,
                  const Ref<VectorXi> ej,
                  const Ref<MatrixXd> blocks,
                  int block_height,
                  int block_width,
                  int num_blocks,
                  bool sym,
                  Sparse& dst) {
        
    int buffer_size = num_blocks * block_height * block_width;
    
    if (sym)
        buffer_size *= 2;
    
    // Allocate buffer
    vector<Triplet> buffer(buffer_size);

    // Fill in the data
    int k = 0;
    for (int i = 0; i < num_blocks; ++i) {
        for (int jj = 0; jj < block_width; ++jj) {
            for (int ii = 0; ii < block_height; ++ii) {
                buffer[k] = Triplet( ei(i)*block_height + ii, ej(i)*block_width + jj, blocks(ii, i*block_width + jj) );
                k++;
            };
        };
    };

    // Duplicate entries if symmetric flag is true
    if (sym)
        for (int i = 0; i < num_blocks; ++i) {
            for (int jj = 0; jj < block_width; ++jj) {
                for (int ii = 0; ii < block_height; ++ii) {
                    buffer[k] = Triplet( ej(i)*block_width + jj, ei(i)*block_height + ii,  blocks(ii, i*block_width + jj) );
                    k++;
                };
            };
        };
    
    dst.setFromTriplets(buffer.begin(), buffer.end());
};


/**
 * Solves the orthogonal Procrustes problem in SO(3) (in-place).
 *
 * @param mat (input/output) Eigen::MatrixXd - 3x3 matrix
 */
void orthoProcrustesSO3(Ref<Matrix3d> mat) {
    Eigen::JacobiSVD<Matrix3d> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Matrix3d U = svd.matrixU();
    Matrix3d R = U * svd.matrixV().transpose();
    U.col(2) *= R.determinant();
    mat = U * svd.matrixV().transpose();
};


/**
 * Solves the orthogonal Procrustes problem in O(3) (in-place).
 *
 * @param mat (input/output) Eigen::MatrixXd - 3x3 matrix.
 */
void orthoProcrustesO3(Ref<Matrix3d> mat) {
    Eigen::JacobiSVD<Matrix3d> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat = svd.matrixU() * svd.matrixV().transpose();
};


/**
 * Finds the projection of a vector on a subpace spanned by the columns of a matrix
 *
 * @param subspace_basis (input) Eigen::MatrixXd - Dense matrix whose columns span a subspace where we want to project the vector v.
 * @param idx (input) unsigned int - Sub-selects of the first idx columns of A only.
 * @param v (input) Eigen::VectorXd - Vector to project.
 * @return Eigen::VectorXd - The component of the vector in the subspace spanned by the first
 *         idx columns of the input matrix A
 */
VectorXd projectionOnto(const Ref<MatrixXd> subspace_basis, unsigned int idx, const Ref<VectorXd> v) {
    VectorXd w = subspace_basis.leftCols(idx).transpose() * v;
    VectorXd proj = subspace_basis.leftCols(idx) * w;
    return proj;
};


/**
 * Projects a vector onto a subpace spanned by the columns of a matrix (in-place)
 *
 * @param subspace_basis (input) Eigen::MatrixXd - Dense matrix whose columns span a subspace where we want to project the vector v.
 * @param idx (input) unsigned int - Sub-selects of the first idx columns of A only.
 * @param v (input/output) Eigen::VectorXd - Vector to project.
 */
void projectOnto(const Ref<MatrixXd> subspace_basis, unsigned int idx, Ref<VectorXd> v) {
    VectorXd w = subspace_basis.leftCols(idx).transpose() * v;
    v -= subspace_basis.leftCols(idx) * w;
};


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
void sparseSolverLDL(PardisoLDLT* solver, const Sparse& mat, const Ref<MatrixXd> rhs, Ref<MatrixXd> x) {
    solver->factorize(mat);
    x = solver->solve(rhs);
};


/**
 * Robust reorthogonalization function. For an input matrix V and an input
 * vector r, tries to iteratively orthogonalize r against the j first columns
 * of V. This prevents numerical errors and ensures orthogonality. If, by any
 * chance, and after multiple iterations of subtracting projections, the vector
 * has a norm smaller than 1/sqrt(2) of its starting norm, it cannot be
 * reorthogonalized. This method then attempts to find another orthogonal vector
 * through random restarts. The flag stop is activate if everything fails.
 * (See G. W. Stewart 2001)
 *
 * @param V Dense matrix of doubles whose columns contain orthogonal basis vectors
 * @param r Vector to orthogonalize
 * @param residual_norm Residual norm or r after it has been orthogonalized against the first j columns of V.
 * @param idx Index specifying the number of columns of V used in the orthogonalization.
 * @param stop Stop the algorithm flag. Indicates that it is impossible to find, to machine precision, a vector r orthogonal to the first j cols of V.
 */
void reorthogonalize(const Ref<MatrixXd> V, Ref<VectorXd> r, double& residual_norm, unsigned int idx, int& stop) {
    stop = 0;
    double normr0 = r.norm();
    
    projectOnto(V, idx, r);
    
    residual_norm = r.norm();

    // Iteratively tries to reorthogonalize r against the columns of V
    unsigned int num_reorths = 1;
    while((residual_norm <= (1.0 / sqrt(2)) * normr0) && num_reorths < 5) {
        projectOnto(V, idx, r);
        normr0 = residual_norm;
        residual_norm = r.norm();
        num_reorths++;
    };
    
    /* Cannot reorthogonalize: Invariant subspace found. Restart with
     * a new random vector and try another 3 times */
    if (residual_norm <= (1.0 / sqrt(2)) * normr0) {
        residual_norm = 0;
        stop = 1;
        // Try another 3 times with random restarts
        for (int j = 0; j < 3; ++j) {
            r = VectorXd::Random(r.rows(),1);
            projectOnto(V, idx, r);
            r.normalize();
            
            // Reorthogonalize if necessary
            for (unsigned int k = 0; k < 5; ++k) {
                VectorXd Mr = r;
                VectorXd proj = projectionOnto(V, idx, Mr);
                double rMr = sqrt(abs(r.transpose() * Mr));
                
                if (abs(rMr - 1.0f) <= 1e-10) {
                    stop = 0;
                    break;
                };
                
                // Reorthogonalize
                r -= proj;
                r.normalize();
            };
            if (!stop)
                break;
        };
    } else {
        r /= residual_norm;
    };
};


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
void symKrylovSchurLDL(PardisoLDLT*  solver,
                       const Sparse& mat,
                       Ref<VectorXd> eigenvalues,
                       Ref<MatrixXd> eigenvectors,
                       unsigned int  k,
                       double        sigma) {
    
    long n = mat.rows();
    
    // Sparse identity matrix
    Sparse speye(n,n);
    speye.setIdentity();
    
    Sparse shiftInvert = mat;
    shiftInvert -= sigma * speye;
    
    // Get LDLt factorization of A - sigma I
    solver->factorize(shiftInvert);
    
    MKS_ASSERT(solver->info()==Eigen::Success, mksGetStatusString(mksEigenPardisoLDLErr));
    
    // Matrix to hold the orthogonal Krylov basis
    MatrixXd krylov_basis = MatrixXd::Zero(n, KRYLOV_SUBSPACE_MAX_SIZE);
    
    // Size of the V matrix (not actual size but size in use)
    unsigned int size_krylov_basis = 0;
    
    // Stop the algorithm flag (reorthogonalization gone bad)
    int stop = 0;
    
    // Number of selected eigenvalues
    unsigned int k0 = k;

    // Store projection data
    vector<double> Alpha;
    Alpha.reserve(KRYLOV_SUBSPACE_MAX_SIZE);
    
    // Store projection data
    vector<double> Beta;
    Beta.reserve(KRYLOV_SUBSPACE_MAX_SIZE);
    
    MatrixXd c;

    double residual_norm = 0;
    bool restarted = false;
    unsigned int nconv = 0;
    vector<int> converged(k0);

    MatrixXd ritz_vectors;
    MatrixXd ritz_values;

    // Used to compute argsorts
    vector<int> idx_vec(KRYLOV_SUBSPACE_MAX_SIZE);

    // H matrix
    MatrixXd H(KRYLOV_SUBSPACE_MAX_SIZE, KRYLOV_SUBSPACE_MAX_SIZE);

    // Starting vector for the Krylov iterations
    VectorXd v0 = VectorXd::Random(n,1);
    v0.normalize();
    
    // One step of the power iteration
    VectorXd v = solver->solve(v0);
    v.normalize();
    
    // Main loop of the Krylov-Schur method
    for (unsigned int i = 0; i < KRYLOV_MAXITER; ++i) {
        // Loop to build the invariant subspace in V, H
        for (unsigned int j = size_krylov_basis; j < KRYLOV_SUBSPACE_MAX_SIZE; ++j) {
            
            // Store the Krylov vector v from the previous iteration in the Krylov matrix V
            krylov_basis.col(j) = v;
            // Compute Krylov power using Cholesky solver: r = (A-sigma I)^-1 * v
            VectorXd r = solver->solve(v);
            // Component of the new vector r in the direction of the previous vector v
            double alpha = v.transpose() * r;

            // Orthogonalization
            if (j == 0) {
                // 1st iteration: just subtract the projection on the previous vector
                r.noalias() -= alpha * v;
            } else if (restarted) {
                VectorXd w = krylov_basis.leftCols(j+1).transpose() * r;
                r.noalias() -= krylov_basis.leftCols(j+1) * w;
                restarted = false;
            } else {
                // Subtract projection on the previous vector and residuals before that
                r.noalias() -= alpha * v;
                r.noalias() -= residual_norm * krylov_basis.col(j-1);
            };

            /* Robust reorthogonalization of r against the columns of V.
             * This helps prevent numerical errors */
            reorthogonalize(krylov_basis, r, residual_norm, j+1, stop);

            /* Check if it was possible to reorthogonalize. If not, we cannot find the
             * invariant subspace and must exit the program */
            MKS_ASSERT(!stop, mksGetStatusString(mksKrylovReorthogonalizationErr));
            
            // Save projection of r on the the previous vector v (alpha)
            Alpha.push_back(alpha);
            // Save projection of r on the previous previous vector (beta)
            Beta.push_back(residual_norm);
            // Set the current vector v as the orthogonalized vector r
            v = r;
        };

        // Build matrix H
        Eigen::Map<VectorXd> Alpha_vec(Alpha.data(), Alpha.size());
        Eigen::Map<VectorXd> Beta_vec(Beta.data(), Beta.size()-1);

        MatrixXd H1 = MatrixXd::Zero(Alpha.size(), Alpha.size());
        H1.diagonal(-1).array() += Beta_vec.array();
        H1.diagonal(0).array()  += Alpha_vec.array();
        H1.diagonal(1).array()  += Beta_vec.array();
        
        // H matrix
        H = MatrixXd::Zero(KRYLOV_SUBSPACE_MAX_SIZE, KRYLOV_SUBSPACE_MAX_SIZE);
        H.bottomRightCorner(H1.rows(), H1.cols()) += H1;

        if (ritz_values.rows() > 0)
            H.topLeftCorner(ritz_values.rows(),ritz_values.rows()) += ritz_values.asDiagonal();

        if (ritz_values.rows() > 0) {
            for (unsigned int l = 0; l < k; ++l) {
                H(l,k) = c(0,l);
                H(k,l) = c(0,l);
            };
        };

        Alpha.resize(0);
        Beta.resize(0);

        // Compute Ritz pairs (eigenpairs of the H matrix)
        Eigen::SelfAdjointEigenSolver<MatrixXd> eig(H);
        ritz_vectors = eig.eigenvectors();
        ritz_values  = eig.eigenvalues();

        // Compute residuals to find out what has converged
        vector<double> res(ritz_vectors.cols());
        for (unsigned int l = 0; l < ritz_vectors.cols(); ++l)
            res[l] = abs(residual_norm * ritz_vectors(ritz_vectors.rows()-1, l));

        // Argsort eigenvalues descending order
        iota(idx_vec.begin(), idx_vec.end(), 0);
        sort(idx_vec.begin(), idx_vec.end(), [&](size_t a, size_t b) { return abs(ritz_values(a)) > abs(ritz_values(b)); });

        // Check which ones have converged
        for (unsigned int l = 0; l < k0; ++l)
            converged[l] = (int) ( res[idx_vec[l]] < KRYLOV_TOL * max(KRYLOV_RES_TOL, abs(ritz_values(idx_vec[l]))) );

        // Number of eigenvalues which have converged
        nconv = accumulate(converged.begin(), converged.end(), 0);
                
        // More than k0 eigenvalues converged: We're good to go
        if (nconv >= k0) {
            break;
        } else {
            // Adjust k to prevent stagnation
            k = k0 + min( (float) nconv, (float) floor(((float) (KRYLOV_SUBSPACE_MAX_SIZE-k0)) / 2.0f) );
            if (k == 1)
                k = floor(KRYLOV_SUBSPACE_MAX_SIZE / 2.0f);
        };

        // Use previous eigenvalue argsort to sort eigenvalues and eigenvectors
        MatrixXd ritz_vectors_ = MatrixXd::Zero(KRYLOV_SUBSPACE_MAX_SIZE, k);
        VectorXd ritz_values_  = VectorXd::Zero(KRYLOV_SUBSPACE_MAX_SIZE, 1);
        
        for (unsigned int i = 0; i < k; ++i) {
            ritz_vectors_.col(i) += ritz_vectors.col(idx_vec[i]);
            ritz_values_(i) += ritz_values(idx_vec[i]);
        };
        
        ritz_vectors = ritz_vectors_;
        ritz_values = ritz_values_;

        // Store variables for the next iteration
        krylov_basis.leftCols(k) = krylov_basis * ritz_vectors;
        c = residual_norm * ritz_vectors.bottomRows<1>();

        restarted = true;
        size_krylov_basis = k;
    };
    
    MatrixXd ritz_vectors_ = MatrixXd::Zero(ritz_vectors.rows(), k0);
    eigenvalues = VectorXd::Zero(k0);

    for (unsigned int i = 0; i < k0; ++i) {
        eigenvalues(i) += ritz_values(idx_vec[i]);
        ritz_vectors_.col(i) += ritz_vectors.col(idx_vec[i]);
    };

    c = residual_norm * ritz_vectors_.bottomRows<1>();

    eigenvectors = krylov_basis * ritz_vectors_;
    eigenvectors = eigenvectors.array().rowwise() * eigenvalues.transpose().array();
    eigenvectors.noalias() += v * c;

    // Eigenvector normalization
    VectorXd norms = eigenvectors.colwise().norm();
    eigenvectors.array().rowwise() /= norms.array().transpose();

    // Eigenvalue reverse the shift and invert
    eigenvalues.array() = 1.0f / eigenvalues.array() + sigma;
};

};

