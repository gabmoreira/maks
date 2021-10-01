/**
 * ravess_example.cpp
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
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "io.hpp"
#include "ravg.hpp"
#include "linalg.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXi;

typedef Eigen::SparseMatrix<double> Sparse;

int main(int argc, char *argv[]) {
    
    int MAX_EDGES = 25000;
    
    std::string path(argv[1]);
    
    int num_nodes;
    int num_edges;
    
    VectorXi edge_i(MAX_EDGES);
    VectorXi edge_j(MAX_EDGES);
    
    // Rotation data allocation in contiguous row matrix blocks
    MatrixXd edge_r(3, 3 * MAX_EDGES);

    // Read edge rotations from .g2o file
    readG2O(path.c_str(), edge_i, edge_j, edge_r, num_edges);
    
    // Convert edges ids to matrix indices in place
    alg::convertToIdx(edge_i.topRows(num_edges), edge_j.topRows(num_edges), num_nodes);
        
    printf("Read %d nodes and %d edges from %s\n", num_nodes, num_edges, path.c_str());
    
    Sparse A(num_nodes, num_nodes);
    alg::adjacency(edge_i, edge_j, num_edges, true, A);
    
    Sparse Rtilde(num_nodes * 3, num_nodes * 3);
    alg::sparseBlocks(edge_i, edge_j, edge_r, 3, 3, num_edges, true, Rtilde);
        
    // Placeholder for the solution
    MatrixXd R(3 * num_nodes, 3);
    double dual = 0.0f;
    
    // Primal-dual method for rotation averaging
    primalDualSO3(Rtilde, A, R, num_nodes, 50, dual, 5e-14, -1e-6);
     
};
    
