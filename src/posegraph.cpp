/*******************************************************************************
* posegraph.cpp
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
#include <map>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <numeric>

#include <Eigen/Geometry>
#include <Eigen/PardisoSupport>

#include "internal.hpp"
#include "krylov.hpp"
#include "posegraph.hpp"

using namespace std::chrono;

using std::iota;
using std::cout;
using std::endl;
using std::sort;
using std::unique;
using std::vector;
using std::string;

using Eigen::Ref;
using Eigen::Matrix4d;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::Vector3d;

#define MAX_NODES 11000
#define MAX_PAIRWISE 100000
#define SE3_SIZE 4
#define SO3_SIZE 3

typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> Sparse;
typedef Eigen::PardisoLDLT<Sparse> MKL_LDLT;


/*******************************************************************************
 * Default constructor of the Graph class. Every graph is initialized with 0
 * nodes and 0 edges. To build the graph the method addEdge should be called.
 ******************************************************************************/
Graph::Graph() : pairwise(SE3_SIZE, SE3_SIZE*MAX_PAIRWISE) {
  ready = false;
  num_nodes = 0;
  num_edges = 0;

  nodes.reserve(MAX_NODES);
  edges_i.reserve(MAX_PAIRWISE);
  edges_j.reserve(MAX_PAIRWISE);

  pairwise = MatrixXd::Zero(SE3_SIZE, SE3_SIZE*MAX_PAIRWISE);
};


/*******************************************************************************
 * Checks if the graph contains the edge (ei, ej). This method checks both
 * directions, i.e. (i,j) is assumed equal to (j,i). Even though the associated
 * transformations are different (inverse of one another), they contain the
 * same information.
 ******************************************************************************/
bool Graph::containsEdge(unsigned int ei, unsigned int ej) {
  for (unsigned int i = 0; i < num_edges; ++i) {
    if (((edges_i[i] == ei) && (edges_j[i] == ej)) || ((edges_i[i] == ej) && (edges_j[i] == ei)))
      return true;
  };
  return false;
};


/*******************************************************************************
 * Adds edge to the pose graph. The first three arguments are necessary for
 * carrying out the pose graph optimization routine. However, the last two may
 * be set to -1 if no information is available.
 ******************************************************************************/
void Graph::addEdge(unsigned int ei, unsigned int ej, const Matrix4d& pairwise_ij) {
  edges_i.push_back(ei);                                                        // Node from (direction matters!)
  edges_j.push_back(ej);                                                        // Node to (direction matters!)
  pairwise.block<SE3_SIZE,SE3_SIZE>(0,num_edges*SE3_SIZE) += pairwise_ij;       // Add transformation as a 4x4 Eigen::Matrix4d
  num_edges++;
};


/*******************************************************************************
 * The Graph class is only updated by adding edges. As such, to generate data on
 * the existing nodes, this method needs to be called. It updates num_nodes,
 * the nodes vector, and the nodes_idx_map which generates, for each node, a
 * unique index in the interval [0, num_nodes].
 ******************************************************************************/
void Graph::updateNodes() {

  nodes.resize(0);
  nodes.insert(nodes.end(), edges_i.begin(), edges_i.end());
  nodes.insert(nodes.end(), edges_j.begin(), edges_j.end());

  sort(nodes.begin(), nodes.end());
  auto last = unique(nodes.begin(), nodes.end());
  nodes.erase(last, nodes.end());

  num_nodes = nodes.size();

  /* Create index map so that we can reference nodes based on their index rather
   * than on the actual node number */
  for (unsigned int i = 0; i < num_nodes; ++i) { nodes_idx_map[nodes[i]] = i; };
};


/*******************************************************************************
 * Builds the graph adjacency matrix as a Eigen::SparseMatrix<double> object
 * with dimensions num_nodes x num_nodes.
 ******************************************************************************/
Sparse Graph::adjacency() {
  Sparse A(num_nodes, num_nodes);

  vector<T> entries;
  entries.reserve(num_edges*2);

  for (unsigned int i = 0; i < num_edges; ++i) {
    entries.push_back(T(nodes_idx_map[edges_i[i]], nodes_idx_map[edges_j[i]], 1));
    entries.push_back(T(nodes_idx_map[edges_j[i]], nodes_idx_map[edges_i[i]], 1));
  };

  /* Build sparse matrix from triples (i,j,entry) -- fastest way */
  A.setFromTriplets(entries.begin(), entries.end());
  return A;
};


/*******************************************************************************
 * Builds the graph Laplacian as a Eigen::SparseMatrix<double> object
 * with dimensions num_nodes x num_nodes.
 ******************************************************************************/
Sparse Graph::laplacian() {
  Sparse L(num_nodes, num_nodes);

  vector<T> entries;
  entries.reserve(num_edges*2 + num_nodes);

  for (unsigned int i = 0; i < num_edges; ++i) {
    entries.push_back(T(nodes_idx_map[edges_i[i]], nodes_idx_map[edges_j[i]], -1));
    entries.push_back(T(nodes_idx_map[edges_j[i]], nodes_idx_map[edges_i[i]], -1));
  };

  L.setFromTriplets(entries.begin(), entries.end());

  VectorXd D = - L * VectorXd::Ones(num_nodes);

  /* Main diagonal */
  for (unsigned int i = 0; i < num_nodes; ++i) { entries.push_back(T(i, i, D(i))); };

  /* Build sparse matrix from triples (i,j,entry) -- fastest way */
  L.setFromTriplets(entries.begin(), entries.end());
  return L;
};


/*******************************************************************************
 * Builds the graph transformation block matrix.
 ******************************************************************************/
Sparse Graph::transformationMatrix() {
  Sparse M(SE3_SIZE*num_nodes, SE3_SIZE*num_nodes);

  vector<T> entries;
  entries.reserve(edges_i.size()*2*13 + num_nodes*SE3_SIZE);

  for (unsigned int i = 0; i < edges_i.size(); ++i) {
    int ii = nodes_idx_map[edges_i[i]];
    int ji = nodes_idx_map[edges_j[i]];

    // Available edges from i to j
    entries.push_back(T(ii*SE3_SIZE,   ji*SE3_SIZE,   pairwise(0,i*SE3_SIZE)));
    entries.push_back(T(ii*SE3_SIZE,   ji*SE3_SIZE+1, pairwise(0,i*SE3_SIZE+1)));
    entries.push_back(T(ii*SE3_SIZE,   ji*SE3_SIZE+2, pairwise(0,i*SE3_SIZE+2)));
    entries.push_back(T(ii*SE3_SIZE,   ji*SE3_SIZE+3, pairwise(0,i*SE3_SIZE+3)));
    entries.push_back(T(ii*SE3_SIZE+1, ji*SE3_SIZE,   pairwise(1,i*SE3_SIZE)));
    entries.push_back(T(ii*SE3_SIZE+1, ji*SE3_SIZE+1, pairwise(1,i*SE3_SIZE+1)));
    entries.push_back(T(ii*SE3_SIZE+1, ji*SE3_SIZE+2, pairwise(1,i*SE3_SIZE+2)));
    entries.push_back(T(ii*SE3_SIZE+1, ji*SE3_SIZE+3, pairwise(1,i*SE3_SIZE+3)));
    entries.push_back(T(ii*SE3_SIZE+2, ji*SE3_SIZE,   pairwise(2,i*SE3_SIZE)));
    entries.push_back(T(ii*SE3_SIZE+2, ji*SE3_SIZE+1, pairwise(2,i*SE3_SIZE+1)));
    entries.push_back(T(ii*SE3_SIZE+2, ji*SE3_SIZE+2, pairwise(2,i*SE3_SIZE+2)));
    entries.push_back(T(ii*SE3_SIZE+2, ji*SE3_SIZE+3, pairwise(2,i*SE3_SIZE+3)));
    entries.push_back(T(ii*SE3_SIZE+3, ji*SE3_SIZE+3, 1));

    // Inverse of the previous edges computed explicitly (fast)
    entries.push_back(T(ji*SE3_SIZE, ii*SE3_SIZE,   pairwise(0,i*SE3_SIZE)));
    entries.push_back(T(ji*SE3_SIZE, ii*SE3_SIZE+1, pairwise(1,i*SE3_SIZE)));
    entries.push_back(T(ji*SE3_SIZE, ii*SE3_SIZE+2, pairwise(2,i*SE3_SIZE)));
    entries.push_back(T(ji*SE3_SIZE, ii*SE3_SIZE+3,
       - pairwise(0,i*SE3_SIZE)*pairwise(0,i*SE3_SIZE+3)
       - pairwise(1,i*SE3_SIZE)*pairwise(1,i*SE3_SIZE+3)
       - pairwise(2,i*SE3_SIZE)*pairwise(3,i*SE3_SIZE+3)));

    entries.push_back(T(ji*SE3_SIZE+1, ii*SE3_SIZE,   pairwise(0,i*SE3_SIZE+1)));
    entries.push_back(T(ji*SE3_SIZE+1, ii*SE3_SIZE+1, pairwise(1,i*SE3_SIZE+1)));
    entries.push_back(T(ji*SE3_SIZE+1, ii*SE3_SIZE+2, pairwise(2,i*SE3_SIZE+1)));
    entries.push_back(T(ji*SE3_SIZE+1, ii*SE3_SIZE+3,
      - pairwise(0,i*SE3_SIZE+1)*pairwise(0,i*SE3_SIZE+3)
      - pairwise(1,i*SE3_SIZE+1)*pairwise(1,i*SE3_SIZE+3)
      - pairwise(2,i*SE3_SIZE+1)*pairwise(3,i*SE3_SIZE+3)));

    entries.push_back(T(ji*SE3_SIZE+2, ii*SE3_SIZE,   pairwise(0,i*SE3_SIZE+2)));
    entries.push_back(T(ji*SE3_SIZE+2, ii*SE3_SIZE+1, pairwise(1,i*SE3_SIZE+2)));
    entries.push_back(T(ji*SE3_SIZE+2, ii*SE3_SIZE+2, pairwise(2,i*SE3_SIZE+2)));
    entries.push_back(T(ji*SE3_SIZE+2, ii*SE3_SIZE+3,
      - pairwise(0,i*SE3_SIZE+2)*pairwise(0,i*SE3_SIZE+3)
      - pairwise(1,i*SE3_SIZE+2)*pairwise(1,i*SE3_SIZE+3)
      - pairwise(2,i*SE3_SIZE+2)*pairwise(3,i*SE3_SIZE+3)));

    entries.push_back(T(ji*SE3_SIZE+3, ii*SE3_SIZE+3, 1));
  };

  /* Main diagonal of ones */
  for (unsigned int i = 0; i < SE3_SIZE*num_nodes; ++i) { entries.push_back(T(i, i, 1)); };

  /* Build sparse matrix from triples (i,j,entry) -- fastest way */
  M.setFromTriplets(entries.begin(), entries.end());
  return M;
};



/*******************************************************************************
 * Pose graph optimization or Motion averaging in SE(3) computed in two stages:
 * 1. Rotation averaging in close-form via eigendecomposition (Krylov subspace
 * method); 2. Optimizing for translations by solving a least squares symmetric
 * indefinite problem via Cholesky's LDL^T factorization.
 ******************************************************************************/
void Graph::optimize(Ref<MatrixXd> evecs, double sigma, double& cf) {
  printf("Running MAKS solver... ");

  MKL_LDLT choleskySolver;

  auto start = high_resolution_clock::now();

  Sparse Adj = adjacency();
  VectorXd D = Adj * VectorXd::Ones(num_nodes);

  auto stop = high_resolution_clock::now();
  auto duration_graph = duration_cast<microseconds>(stop - start);

  start = high_resolution_clock::now();
  /* BUILD LAPLACIAN NORMALIZED ROTATIONS MATRIX */
  vector<T> lnrEntries;
  lnrEntries.reserve(num_edges*18 + num_nodes*SO3_SIZE);

  /* Specify entries of the lnr matrix based on triplets */
  for (unsigned int i = 0; i < num_edges; ++i) {
    int ii = nodes_idx_map[edges_i[i]];
    int ji = nodes_idx_map[edges_j[i]];

    /* Available edges from i to j */
    lnrEntries.push_back(T(ii*SO3_SIZE,   ji*SO3_SIZE,   -pairwise(0,i*SE3_SIZE)));
    lnrEntries.push_back(T(ii*SO3_SIZE,   ji*SO3_SIZE+1, -pairwise(0,i*SE3_SIZE+1)));
    lnrEntries.push_back(T(ii*SO3_SIZE,   ji*SO3_SIZE+2, -pairwise(0,i*SE3_SIZE+2)));
    lnrEntries.push_back(T(ii*SO3_SIZE+1, ji*SO3_SIZE,   -pairwise(1,i*SE3_SIZE)));
    lnrEntries.push_back(T(ii*SO3_SIZE+1, ji*SO3_SIZE+1, -pairwise(1,i*SE3_SIZE+1)));
    lnrEntries.push_back(T(ii*SO3_SIZE+1, ji*SO3_SIZE+2, -pairwise(1,i*SE3_SIZE+2)));
    lnrEntries.push_back(T(ii*SO3_SIZE+2, ji*SO3_SIZE,   -pairwise(2,i*SE3_SIZE)));
    lnrEntries.push_back(T(ii*SO3_SIZE+2, ji*SO3_SIZE+1, -pairwise(2,i*SE3_SIZE+1)));
    lnrEntries.push_back(T(ii*SO3_SIZE+2, ji*SO3_SIZE+2, -pairwise(2,i*SE3_SIZE+2)));

    /* Inverse of the previous edges computed explicitly (fast) */
    lnrEntries.push_back(T(ji*SO3_SIZE,   ii*SO3_SIZE,   -pairwise(0,i*SE3_SIZE)));
    lnrEntries.push_back(T(ji*SO3_SIZE,   ii*SO3_SIZE+1, -pairwise(1,i*SE3_SIZE)));
    lnrEntries.push_back(T(ji*SO3_SIZE,   ii*SO3_SIZE+2, -pairwise(2,i*SE3_SIZE)));
    lnrEntries.push_back(T(ji*SO3_SIZE+1, ii*SO3_SIZE,   -pairwise(0,i*SE3_SIZE+1)));
    lnrEntries.push_back(T(ji*SO3_SIZE+1, ii*SO3_SIZE+1, -pairwise(1,i*SE3_SIZE+1)));
    lnrEntries.push_back(T(ji*SO3_SIZE+1, ii*SO3_SIZE+2, -pairwise(2,i*SE3_SIZE+1)));
    lnrEntries.push_back(T(ji*SO3_SIZE+2, ii*SO3_SIZE,   -pairwise(0,i*SE3_SIZE+2)));
    lnrEntries.push_back(T(ji*SO3_SIZE+2, ii*SO3_SIZE+1, -pairwise(1,i*SE3_SIZE+2)));
    lnrEntries.push_back(T(ji*SO3_SIZE+2, ii*SO3_SIZE+2, -pairwise(2,i*SE3_SIZE+2)));
  };

  // Main diagonal of lnr
  unsigned int j = 0;
  for (unsigned int i = 0; i < num_nodes; ++i) {
    lnrEntries.push_back(T(j, j, D(i))); j++;
    lnrEntries.push_back(T(j, j, D(i))); j++;
    lnrEntries.push_back(T(j, j, D(i))); j++;
  };

  Sparse lnr(SO3_SIZE*num_nodes, SO3_SIZE*num_nodes);                           // Laplacian normalized rotations block matrix
  lnr.setFromTriplets(lnrEntries.begin(), lnrEntries.end());                    // Build sparse matrix

  stop = high_resolution_clock::now();
  auto duration_rotmat = duration_cast<microseconds>(stop - start);

  start = high_resolution_clock::now();
  choleskySolver.analyzePattern(lnr);                                           // Cholesky LDL^T analysis
  stop = high_resolution_clock::now();
  auto duration_chol = duration_cast<microseconds>(stop - start);

  /* ROTATION AVERAGING */
  start = high_resolution_clock::now();
  MatrixXd revecs(num_nodes*SO3_SIZE, SO3_SIZE);                                // Rotation averaging eigenvectors
  VectorXd revals(SO3_SIZE);                                                    // Rotation averaging eigenvalues

  kry::eigs(&choleskySolver, lnr, revals, revecs, SO3_SIZE, sigma);             // Restarted Krylov method of symmetric sparse matrix

  /* Compute rotations relative to the first one and project blocks to SO(3) */
  Matrix3d first_rotation = revecs.block<SO3_SIZE,SO3_SIZE>(0,0);
  Matrix3d inv_first_rotation = first_rotation.inverse();
  revecs = revecs * inv_first_rotation;

  /* Project to the Special Orthogonal Group SO(3) */
  for (unsigned int i = 0; i < num_nodes; ++i) { internal::projectToSO3(revecs.block<SO3_SIZE,SO3_SIZE>(i*3,0)); };

  stop = high_resolution_clock::now();
  auto duration_eigs = duration_cast<microseconds>(stop - start);

  /*---------------------------------------------------------------------------*/
  start = high_resolution_clock::now();
  vector<T> AEntries;
  vector<T> BEntries;
  AEntries.reserve(num_edges*18+num_nodes*3);
  BEntries.reserve(num_edges*6);

  /* Specify entries of the lnr matrix based on triplets */
  for (unsigned int i = 0; i < num_edges; ++i) {
    int ii = nodes_idx_map[edges_i[i]];
    int ji = nodes_idx_map[edges_j[i]];

    Matrix3d RiRj = revecs.block<3,3>(ii*3,0) * revecs.block<3,3>(ji*3,0).transpose();
    Vector3d tij = 0.5 * (pairwise.block<3,1>(0,i*4+3) + RiRj * pairwise.block<3,3>(0,i*4).transpose() * pairwise.block<3,1>(0,i*4+3));

    // Translations
    BEntries.push_back( T(ii*SO3_SIZE,   ji, -tij(0)) );
    BEntries.push_back( T(ii*SO3_SIZE+1, ji, -tij(1)) );
    BEntries.push_back( T(ii*SO3_SIZE+2, ji, -tij(2)) );
    BEntries.push_back( T(ji*SO3_SIZE,   ii,  RiRj(0,0)*tij(0) + RiRj(1,0)*tij(1) + RiRj(2,0)*tij(2)) );
    BEntries.push_back( T(ji*SO3_SIZE+1, ii,  RiRj(0,1)*tij(0) + RiRj(1,1)*tij(1) + RiRj(2,1)*tij(2)) );
    BEntries.push_back( T(ji*SO3_SIZE+2, ii,  RiRj(0,2)*tij(0) + RiRj(1,2)*tij(1) + RiRj(2,2)*tij(2)) );

    // Rotations
    AEntries.push_back( T(ii*SO3_SIZE,   ji*SO3_SIZE,   -RiRj(0,0)) );
    AEntries.push_back( T(ii*SO3_SIZE,   ji*SO3_SIZE+1, -RiRj(0,1)) );
    AEntries.push_back( T(ii*SO3_SIZE,   ji*SO3_SIZE+2, -RiRj(0,2)) );
    AEntries.push_back( T(ii*SO3_SIZE+1, ji*SO3_SIZE,   -RiRj(1,0)) );
    AEntries.push_back( T(ii*SO3_SIZE+1, ji*SO3_SIZE+1, -RiRj(1,1)) );
    AEntries.push_back( T(ii*SO3_SIZE+1, ji*SO3_SIZE+2, -RiRj(1,2)) );
    AEntries.push_back( T(ii*SO3_SIZE+2, ji*SO3_SIZE,   -RiRj(2,0)) );
    AEntries.push_back( T(ii*SO3_SIZE+2, ji*SO3_SIZE+1, -RiRj(2,1)) );
    AEntries.push_back( T(ii*SO3_SIZE+2, ji*SO3_SIZE+2, -RiRj(2,2)) );
    // Inverse of the previous edges computed explicitly (fast)
    AEntries.push_back( T(ji*SO3_SIZE,   ii*SO3_SIZE,   -RiRj(0,0)) );
    AEntries.push_back( T(ji*SO3_SIZE,   ii*SO3_SIZE+1, -RiRj(1,0)) );
    AEntries.push_back( T(ji*SO3_SIZE,   ii*SO3_SIZE+2, -RiRj(2,0)) );
    AEntries.push_back( T(ji*SO3_SIZE+1, ii*SO3_SIZE,   -RiRj(0,1)) );
    AEntries.push_back( T(ji*SO3_SIZE+1, ii*SO3_SIZE+1, -RiRj(1,1)) );
    AEntries.push_back( T(ji*SO3_SIZE+1, ii*SO3_SIZE+2, -RiRj(2,1)) );
    AEntries.push_back( T(ji*SO3_SIZE+2, ii*SO3_SIZE,   -RiRj(0,2)) );
    AEntries.push_back( T(ji*SO3_SIZE+2, ii*SO3_SIZE+1, -RiRj(1,2)) );
    AEntries.push_back( T(ji*SO3_SIZE+2, ii*SO3_SIZE+2, -RiRj(2,2)) );
  };

  /* Main constant diagonal of A */
  j = 0;
  for (unsigned int i = 0; i < nodes.size(); i++) {
    AEntries.push_back( T(j, j, D(i)) ); j++;
    AEntries.push_back( T(j, j, D(i)) ); j++;
    AEntries.push_back( T(j, j, D(i)) ); j++;
  };

  Sparse A(num_nodes*SO3_SIZE, num_nodes*SO3_SIZE);                             // Laplacian normalized estimated rotations
  A.setFromTriplets(AEntries.begin(), AEntries.end());                          // Build sparse matrix
  Sparse B(num_nodes*SO3_SIZE, num_nodes);                                      // Intermediary translation estimates
  B.setFromTriplets(BEntries.begin(), BEntries.end());                          // Build sparse matrix

  VectorXd b = -B * VectorXd::Ones(num_nodes, 1);                               // This vector is the right-hand side of Ax = b

  VectorXd translations = VectorXd::Zero(num_nodes*SO3_SIZE);
  kry::sparseSolve(&choleskySolver, A, b, translations);                        // Solve for x : Ax = b. Use the same solver

  stop = high_resolution_clock::now();
  auto duration_trans = duration_cast<microseconds>(stop - start);

  evecs = MatrixXd::Zero(num_nodes*SE3_SIZE, SE3_SIZE);                         // Set up what will be our matrix with global transformations

  /* Copy 3x3 rotations from revecs and the 3x1 translations to evecs */
  for (unsigned int i = 0; i < num_nodes; ++i) {
    evecs.block<SO3_SIZE,SO3_SIZE>(i*4,0) = revecs.block<SO3_SIZE,SO3_SIZE>(i*3,0);
    evecs.block<SO3_SIZE,1>(i*4,3) = translations.block<SO3_SIZE,1>(i*3,0);
  };

  /* Place ones to make the blocks in SE(3) */
  for (unsigned int i = 3; i < num_nodes*4; i+=4) { evecs(i,3) += 1; };

  /* Everything w.r.t to the first pose */
  Matrix4d first_transformation = evecs.block<SE3_SIZE,SE3_SIZE>(0,0);
  evecs = evecs * first_transformation.inverse();
  printf("Done!\n");
  printf("____________________________________________________\n");
  printf("                  Optimization log\n");
  printf("----------------------------------------------------\n");
  printf("Sparse solver  |   Intel MKL PARDISO Cholesky LDL\n");
  printf("----------------------------------------------------\n");
  printf("CPU Time (s)   |             Routine\n");
  printf("----------------------------------------------------\n");
  printf("   %.6f    |    Graph degree vector\n", static_cast<long long int>(duration_graph.count())*1e-6);
  printf("   %.6f    |    Created rotations block matrix\n", static_cast<long long int>(duration_rotmat.count())*1e-6);
  printf("   %.6f    |    Cholesky pattern analysis\n", static_cast<long long int>(duration_chol.count())*1e-6);
  printf("   %.6f    |    Solved rotation averaging\n", static_cast<long long int>(duration_eigs.count())*1e-6);
  printf("   %.6f    |    Optimized for translations\n", static_cast<long long int>(duration_trans.count())*1e-6);
  printf("____________________________________________________\n");
};


/*******************************************************************************
 * Reads edges (rotation expressed as a unit quaternion and translation expressed
 * as a vector in R^3) and adds them to the graph. Some .g2o files contain,
 * besides edges, vertices which may be used to initialize PGO algorithms.
 * Our method does not need to be initialized so this information is discarded.
 ******************************************************************************/
void Graph::readG2O(const char* path) {
  printf("\n____________________________________________________\n");
  printf("                   Loading data \n");
  printf("----------------------------------------------------\n");
  printf("File : %s\n", path);
  std::ifstream infile(path);

  if(!infile) { return; };

  string edge_tag = "EDGE_SE3:QUAT";
  string line;

  while (std::getline(infile, line)) {
    vector<string> tokens = internal::strSplit(line, " ");
    string tag = tokens[0];

    if (0 == tag.compare(edge_tag)) {
      int vertex_observing_id = std::stoi(tokens[1]);
      int vertex_observed_id  = std::stoi(tokens[2]);
      double edge_x = std::stod(tokens[3]);
      double edge_y = std::stod(tokens[4]);
      double edge_z = std::stod(tokens[5]);

      Eigen::Quaterniond q;
      q.x() = std::stod(tokens[6]);
      q.y() = std::stod(tokens[7]);
      q.z() = std::stod(tokens[8]);
      q.w() = std::stod(tokens[9]);

      Matrix3d edge_R = q.normalized().toRotationMatrix();
      Matrix4d edge_se3 = MatrixXd::Identity(4,4);

      edge_se3.block<3,3>(0,0) = edge_R;
      edge_se3(0,3) = edge_x;
      edge_se3(1,3) = edge_y;
      edge_se3(2,3) = edge_z;

      addEdge(vertex_observing_id, vertex_observed_id, edge_se3);
    };
  };

  printf("Relative poses :   %d\n", num_edges);
  updateNodes();
  printf("Nodes :            %d\n", num_nodes);
  printf("____________________________________________________\n\n");
  ready = true;
};


/*******************************************************************************
 * Writes optimized pose graph to disk in .g2o format. Filename convention uses
 * the original dataset's name + "_maks.g2o"
 ******************************************************************************/
void Graph::writeG2O(const char* dataset, Eigen::Ref<MatrixXd> estimate) {
  string filepath = string(dataset) + "_maks.g2o";

  std::ostringstream buffer;
  buffer.clear();

  for (int i = 0; i < num_nodes; ++i) {
    buffer << "VERTEX_SE3:QUAT";
    buffer << " " << std::to_string(i);
    buffer << " " << std::fixed << std::setprecision(6) << estimate(i*4,3);   // X
    buffer << " " << std::fixed << std::setprecision(6) << estimate(i*4+1,3); // Y
    buffer << " " << std::fixed << std::setprecision(6) << estimate(i*4+2,3); // Z

    Eigen::Quaterniond Qij(estimate.block<3,3>(i*4,0));

    buffer << " " << std::fixed << std::setprecision(6) << Qij.x();           // QX
    buffer << " " << std::fixed << std::setprecision(6) << Qij.y();           // QY
    buffer << " " << std::fixed << std::setprecision(6) << Qij.z();           // QZ
    buffer << " " << std::fixed << std::setprecision(6) << Qij.w() << "\n";   // QW
  };

  std::ofstream outfile;
  outfile.open(filepath);
  if(!outfile) {
    return;
  };

  outfile << buffer.str();
  outfile.close();
  printf("\nOptimized pose graph saved to %s\n", filepath.c_str());
};
