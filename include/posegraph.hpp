/*******************************************************************************
* posegraph.hpp
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
#ifndef POSEGRAPH_H
#define POSEGRAPH_H

#define EIGEN_USE_MKL_ALL

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

/*******************************************************************************
 * Implementation of the Graph class. It is used to support the analysis, optimization
 * and manipulation of pose graphs.
 ******************************************************************************/
class Graph {
  public:
    /*******************************************************************************
     * Default constructor
     *
     * @return Nothing
     ******************************************************************************/
    Graph();


    /*******************************************************************************
     * Adds edge (Node from = ei, Node to = ej) to the graph. Each edge has an
     * associated rigid transformation pairwise_ij (matrix in SE(3)) a RMSE and a
     * number of inliers.
     *
     * @param ei Node from.
     * @param ej Node to.
     * @param pairwise_ij Rigid transformation as a 4x4 matrix in SE(3).
     * @param rmse_ij RMSE of this observation.
     * @return Rien
     ******************************************************************************/
    void addEdge(unsigned int ei, unsigned int ej, const Eigen::Matrix4d& pairwise_ij);


    /*******************************************************************************
     * Checks if the graph has the edge (ei, ej). Also checks for (ej,ei) since they
     * are considered one and the same.
     *
     * @param ei Node from.
     * @param ej Node to.
     * @return Boolean specifying whether the edge is present.
     ******************************************************************************/
    bool containsEdge(unsigned int ei, unsigned int ej);


    /*******************************************************************************
     * Retrives a vector containing the nodes, as specified by their original ids,
     * If the graph as a single edge (10,11), the nodes will be [10, 11].
     *
     * @return Nothing
     ******************************************************************************/
    void updateNodes();


    /*******************************************************************************
     * Builds and returns the adjacency matrix of the graph.
     *
     * @return Sparse adjacency matrix.
     ******************************************************************************/
    Eigen::SparseMatrix<double> adjacency();


    /*******************************************************************************
     * Builds and returns the Laplacian matrix of the graph.
     *
     * @return Sparse Laplacian matrix.
     ******************************************************************************/
    Eigen::SparseMatrix<double> laplacian();


    /*******************************************************************************
     * Builds and returns the pairwise transformations block matrix of the graph. The
     * 4x4 block (i,j) contains the observed pairwise transformation corresponding
     * to the i-th and j-th nodes (not actually the node's ids). E.g. a graph has 2
     * edges (10,11) and (11,20). The block (1,2) corresponds to the observed
     * transformation between the nodes 10 and 11 (1st and 2nd nodes).
     *
     * @return Relative transformations sparse-block matrix.
     ******************************************************************************/
    Eigen::SparseMatrix<double> transformationMatrix();


    /*******************************************************************************
     * Pose graph optimization in SE(3)
     *
     * @return Nothing
     ******************************************************************************/
    void optimize(Eigen::Ref<Eigen::MatrixXd> evecs, double sigma, double& cf);


    /*******************************************************************************
     * Loads pairwise observations (edges) from a .g2o file.
     *
     * @return Nothing.
     ******************************************************************************/
    void readG2O(const char* path);

    bool ready;
    
    /* Number of edges. (i,j) and (j,i) are only counted once */
    unsigned int num_edges;

    /* Number of nodes */
    unsigned int num_nodes;

    /* Nodes */
    std::vector<unsigned int> nodes;
    std::map<unsigned int, unsigned int> nodes_idx_map;

    /* Nodes from */
    std::vector<unsigned int> edges_i;

    /* Nodes to */
    std::vector<unsigned int> edges_j;

    /* Matrix whose 4x4 blocks contain the transformations associated with each edge */
    Eigen::MatrixXd pairwise;
};
#endif
