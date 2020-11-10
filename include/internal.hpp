/*******************************************************************************
* internal.hpp
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
#ifndef INTERNAL_H
#define INTERNAL_H

#define EIGEN_USE_MKL_ALL

#include <vector>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace internal {

/*******************************************************************************
 * Splits a string based on a delimiter
 *
 * @param str Input string to be split
 * @param delim Input string delimiter
 * @return Split tokens
 *******************************************************************************/
std::vector<std::string> strSplit(const std::string& str, const std::string& delim);

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
                            unsigned int idx);

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
                     unsigned int idx, int& stop);

/*******************************************************************************
 * Computes the 3D rotation R in SO(3) closest to mat in the least squares sense.
 * The result is stored in place (no return).
 *
 * @return Nothing.
 ******************************************************************************/
void projectToSO3(Eigen::Ref<Eigen::MatrixXd> mat);

};
#endif
