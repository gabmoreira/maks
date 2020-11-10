/*******************************************************************************
* 2020 Gabriel A. Moreira
*
* g.antunes.moreira at gmail dot com
* https://github.com/gabmoreira
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>

#include "posegraph.hpp"

using std::string;
using std::vector;

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(int argc, char *argv[]) {
  std::string dataset;
  if (argc == 2) {
    dataset += argv[1];
  } else {
    printf("Error: No dataset. See manual for command line arguments.\n");
    return 1;
  };

  Graph g;
  g.readG2O(dataset.c_str());

  if (g.ready) {
    double cf = -1.0;
    MatrixXd evecs = MatrixXd::Zero(g.num_nodes*4,4);
    g.optimize(evecs, 1e-10, cf);
    return 0;
  } else {
    printf("Error reading file %s\n", dataset.c_str());
    return 1;
  };
};
