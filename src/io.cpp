/**
 * io.cpp
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
#include <vector>
#include <fstream>

#include <Eigen/Geometry>

#include "io.hpp"

using std::vector;
using std::string;
using std::stod;
using std::stoi;

using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXi;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Ref;

#define G2O_PRECISION 15

/**
 * Splits a string based on a delimiter.
 *
 * @param str Input string to be split.
 * @param delim Input string delimiter.
 * @return vector of tokens.
 */
vector<string> strSplit(const string& str, const string& delim) {
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do {
        pos = str.find(delim, prev);
        if (pos == string::npos)
            pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty())
            tokens.push_back(token);
        prev = pos + delim.length();
    } while (pos < str.length() && prev < str.length());

    return tokens;
};


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
void readG2O(const char*   path,
             Ref<VectorXi> node_i,
             Ref<MatrixXd> node_r,
             Ref<MatrixXd> node_t,
             Ref<VectorXi> edge_i,
             Ref<VectorXi> edge_j,
             Ref<MatrixXd> edge_r,
             Ref<MatrixXd> edge_t,
             int&          num_nodes,
             int&          num_edges) {
    
    std::ifstream infile(path);
    if(!infile)
        return;

    string edge_tag = "EDGE_SE3:QUAT";
    string node_tag = "VERTEX_SE3:QUAT";

    num_nodes = 0;
    num_edges = 0;
    
    string line;
    while (std::getline(infile, line)) {
        vector<string> tokens = strSplit(line, " ");
        string tag = tokens[0];

        if (0 == tag.compare(edge_tag)) {
            int edge_i_id = stoi(tokens[1]);
            int edge_j_id = stoi(tokens[2]);
            
            // Edge translation
            Vector3d t( stod(tokens[3]), stod(tokens[4]), stod(tokens[5]) );
            
            // Edge quaternion
            Quaterniond q( stod(tokens[9]), stod(tokens[6]), stod(tokens[7]), stod(tokens[8]) );

            edge_i(num_edges) += edge_i_id;
            edge_j(num_edges) += edge_j_id;
            edge_r.block<3,3>(0,num_edges*3) += q.normalized().toRotationMatrix();
            edge_t.block<3,1>(0,num_edges) += t;
            num_edges++;
                        
        } else if (0 == tag.compare(node_tag)) {
            int node_id = stoi(tokens[1]);
            
            // Node translation
            Vector3d t( stod(tokens[2]), stod(tokens[3]), stod(tokens[4]) );

            // Node quaternion
            Quaterniond q( stod(tokens[8]), stod(tokens[5]), stod(tokens[6]), stod(tokens[7]) );

            node_i(num_nodes) += node_id;
            node_r.block<3,3>(0,num_nodes*3) += q.normalized().toRotationMatrix();
            node_t.block<3,1>(0,num_nodes) += t;
            num_nodes++;
        };
    };
};


/**
 * Reads edge rotations from g2o file.
 *
 * @param path (input) char*  - g2o file name.
 * @param edge_i (output) Eigen::VectorXi - edge id (from).
 * @param edge_j (output) Eigen::VectorXi - edge id (to).
 * @param edge_r (output) Eigen::MatrixXd - Row block matrix containing SO(3) rotations.
 * @param num_edges (output) int - number of edges.
 */
void readG2O(const char*          path,
             Eigen::Ref<VectorXi> edge_i,
             Eigen::Ref<VectorXi> edge_j,
             Eigen::Ref<MatrixXd> edge_r,
             int&                 num_edges) {
    
    std::ifstream infile(path);
    if(!infile)
        return;

    string edge_tag = "EDGE_SE3:QUAT";

    num_edges = 0;
    
    string line;
    while (std::getline(infile, line)) {
        vector<string> tokens = strSplit(line, " ");
        string tag = tokens[0];

        if (0 == tag.compare(edge_tag)) {
            int edge_i_id = stoi(tokens[1]);
            int edge_j_id = stoi(tokens[2]);
            
            /* Edge quaternion */
            Quaterniond q( stod(tokens[9]), stod(tokens[6]), stod(tokens[7]), stod(tokens[8]) );

            edge_i(num_edges) += edge_i_id;
            edge_j(num_edges) += edge_j_id;
            edge_r.block<3,3>(0,num_edges*3) += q.normalized().toRotationMatrix();
            num_edges++;
        };
    };
};




void writeG2O(const char* filename,
              const vector<int>&         node_i,
              const vector<Quaterniond>& node_q,
              const vector<Vector3d>&    node_t,
              const vector<int>&         edge_i,
              const vector<int>&         edge_j,
              const vector<Quaterniond>& edge_q,
              const vector<Vector3d>&    edge_t) {

    int num_nodes = (int) node_i.size();
    int num_edges = (int) edge_i.size();
    
    string filepath = string(filename) + ".g2o";

    std::ostringstream buffer;
    buffer.clear();

    /* Write node data */
    for (int i = 0; i < num_nodes; ++i) {
        /* Write tag */
        buffer << "VERTEX_SE3:QUAT";
        
        /* Write node id */
        buffer << " " << std::to_string(node_i[i]);
        
        /* Write node translation */
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << node_t[i](0); // x
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << node_t[i](1); // y
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << node_t[i](2); // z

        /* Write node quaternion */
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << node_q[i].x();           // qx
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << node_q[i].y();           // qy
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << node_q[i].z();           // qz
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << node_q[i].w() << "\n";   // qw
    };
    
    /* Write edge data */
    for (int i = 0; i < num_edges; ++i) {
        /* Write tag */
        buffer << "EDGE_SE3:QUAT";
        
        /* Write edge ids */
        buffer << " " << std::to_string(edge_i[i]);
        buffer << " " << std::to_string(edge_j[i]);

        /* Write node translation */
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << edge_t[i](0); // x
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << edge_t[i](1); // y
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << edge_t[i](2); // z

        /* Write node quaternion */
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << edge_q[i].x();           // qx
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << edge_q[i].y();           // qy
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << edge_q[i].z();           // qz
        buffer << " " << std::fixed << std::setprecision(G2O_PRECISION) << edge_q[i].w() << "\n";   // qw
    };

    std::ofstream outfile;
    outfile.open(filepath);
    if(!outfile) {
        return;
    };

    outfile << buffer.str();
    outfile.close();
    printf("Pose graph saved to %s\n", filepath.c_str());
};
