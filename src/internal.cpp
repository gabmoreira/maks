/*******************************************************************************
* internal.cpp
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
#include <vector>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::Ref;
using Eigen::MatrixXd;
using Eigen::Matrix3d;

namespace internal {

/*******************************************************************************
 * Splits a string based on a delimiter
 *
 * @param str Input string to be split
 * @param delim Input string delimiter
 * @return Split tokens
 *******************************************************************************/
std::vector<std::string> strSplit(const std::string& str, const std::string& delim) {
  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do {
    pos = str.find(delim, prev);
    if (pos == std::string::npos)
      pos = str.length();
    std::string token = str.substr(prev, pos-prev);
    if (!token.empty())
      tokens.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());

  return tokens;
};


/*******************************************************************************
 * Finds the projection of a vector on a subpace spanned by the columns of a matrix
 *
 * @param A Dense matrix whose columns span a subspace where we want to project the
 *          vector v.
 * @param v The vector to be projected on the column-space of A.
 * @param idx Sub-selects of the first idx columns of A only.
 * @return The component of the vector in the subspace spanned by the first
 *         idx columns of the input matrix A
 *******************************************************************************/
Eigen::VectorXd projectOnto(const Eigen::MatrixXd& A, const Eigen::VectorXd& r,
                            unsigned int idx) {
  Eigen::VectorXd w = A.leftCols(idx).transpose() * r;
  Eigen::VectorXd proj = A.leftCols(idx) * w;
  return proj;
};


/*******************************************************************************
 * Robust reorthogonalization function. For an input matrix V and an input
 * vector r, tries to iteratively orthogonalize r against the j first columns
 * of V. This prevents numerical errors and ensures orthogonality. If, by any
 * chance, and after multiple iterations of subtracting projections, the vector
 * has a norm smaller than 1/sqrt(2) of its starting norm, it cannot be
 * reorthogonalized. This method then attempts to find another orthogonal vector
 * through random restarts. The flag stop is activate if everything fails.
 * (See G. W. Stewart 2001)
 *
 * @param V Dense matrix of doubles whose columns contain orthogonal basis
 *          vectors
 * @param r Vector to orthogonalize
 * @param normRes Residual norm or r after it has been orthogonalized against
 *                the first j columns of V.
 * @param idx Index specifying the number of columns of V used in the
 *            orthogonalization.
 * @param stop Stop the algorithm flag. Indicates that it is impossible to find,
 *             to machine precision, a vector r orthogonal to the first j cols
 *             of V.
 * @return Nothing
 ******************************************************************************/
void reorthogonalize(const Eigen::MatrixXd& V, Eigen::VectorXd& r, double& normRes,
                     unsigned int idx, int& stop) {
  double normr0 = r.norm();
  stop = 0;
  r -= projectOnto(V, r, idx);
  normRes = r.norm();

  /* Iteratively tries to reorthogonalize r against the columns of V*/
  unsigned int numReorths = 1;
  while((normRes <= (1.0 / sqrt(2.0))*normr0) && numReorths < 5) {
    r -= projectOnto(V, r, idx);
    normr0 = normRes;
    normRes = r.norm();
    numReorths++;
  };

  /* Cannot reorthogonalize: Invariant subspace found. Restart with
   * a new random vector and try another 3 times */
  if (normRes <= (1.0 / sqrt(2))*normr0) {
    normRes = 0;
    stop = 1;

    // Try another 3 times with random restarts
    for (int j = 0; j < 3; ++j) {
      r = Eigen::VectorXd::Random(r.rows(),1);
      r -= projectOnto(V, r, idx);
      r.normalize();

      // Reorthogonalize if necessary
      for (unsigned int k = 0; k < 5; ++k) {
        Eigen::VectorXd Mr = r;
        Eigen::VectorXd proj = projectOnto(V, Mr, idx);
        double rMr = sqrt(abs(r.transpose() * Mr));

        if (abs(rMr - 1) <= 1e-10) {
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
    r /= normRes;
  };
};

/*******************************************************************************
 * Finds 3x3 rotation matrix closest to the provided matrix (in place)
 *
 * @param mat 3x3 matrix
 * @return Nothing
 *******************************************************************************/
void projectToSO3(Ref<MatrixXd> mat) {
  Eigen::JacobiSVD<MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Matrix3d U = svd.matrixU();
  Matrix3d V = svd.matrixV();

  Matrix3d R = U * V.transpose();
  Matrix3d reflex;
  reflex << 1, 0, 0,
            0, 1, 0,
            0, 0, R.determinant();

  R = reflex * V.transpose();
  mat = U * R;
};

};
