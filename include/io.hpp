//
//  io.hpp
//  maks
//
//  Created by Gabriel Moreira on 09/09/2021.
//  Copyright © 2021 Gabriel Moreira. All rights reserved.
//

#ifndef IO_HPP
#define IO_HPP

#include <Eigen/Dense>

/**
 * Reads g2o file.
 *
 * @param path (input) char*  - g2o file name.
 * @param node_i (output) Eigen::VectorXi - node indices.
 * @param node_r (output) Eigen::MatrixXd - SO(3) node rotations.
 * @param node_t (output) Eigen::MatrixXd - node translations.
 * @param edge_i (output) Eigen::VectorXi - edge indices (from).
 * @param edge_j (output) Eigen::VectorXi - edge indices (to).
 * @param edge_r (output) Eigen::MatrixXd - SO(3) edge rotations.
 * @param edge_t (output) Eigen::MatrixXd - edge translations.
 * @param num_nodes (output) int - number of nodes loaded.
 * @param num_edges (output) int - number of edges loaded.
 */
void readG2O(const char*                 path,
             Eigen::Ref<Eigen::VectorXi> node_i,
             Eigen::Ref<Eigen::MatrixXd> node_r,
             Eigen::Ref<Eigen::MatrixXd> node_t,
             Eigen::Ref<Eigen::VectorXi> edge_i,
             Eigen::Ref<Eigen::VectorXi> edge_j,
             Eigen::Ref<Eigen::MatrixXd> edge_r,
             Eigen::Ref<Eigen::MatrixXd> edge_t,
             int&                        num_nodes,
             int&                        num_edges);


/**
 * Reads edge rotations from g2o file.
 *
 * @param path (input) char*  - g2o file name.
 * @param edge_i (output) Eigen::VectorXi - edge id (from).
 * @param edge_j (output) Eigen::VectorXi - edge id (to).
 * @param edge_r (output) Eigen::MatrixXd - Row block matrix containing SO(3) rotations.
 * @param num_edges (output) int - number of edges.
 */
void readG2O(const char*                 path,
             Eigen::Ref<Eigen::VectorXi> edge_i,
             Eigen::Ref<Eigen::VectorXi> edge_j,
             Eigen::Ref<Eigen::MatrixXd> edge_r,
             int&                        num_edges);

#endif /* IO_HPP */
