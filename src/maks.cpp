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
#include "internal.hpp"

using std::string;
using std::vector;

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(int argc, char *argv[]) {
  std::string filename;

  if (argc == 2) {
    filename += argv[1];
  } else {
    printf("Error: No dataset. See manual for command line arguments.\n");
    return 1;
  };

  Graph g;
  g.readG2O(filename.c_str());

  if (g.ready) {
    double cf = -1.0;               // Cost function
    double sigma = 1e-10;           // Spectral shift

    MatrixXd estimate = MatrixXd::Zero(g.num_nodes*4,4);
    g.optimize(estimate, sigma, cf);
    
    string new_filename = filename.substr(0, filename.length()-4);
    g.writeG2O(new_filename.c_str(), estimate);
    return 0;

  } else {
    printf("Error reading file %s\n", filename.c_str());
    return 1;
  };
};
