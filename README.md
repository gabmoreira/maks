## MAKS - Motion Averaging
MAKS contains routines for pose graph optimization and rotation averaging. Still under development!

### Dependencies
* [Eigen](http://eigen.tuxfamily.org) : Header-only library used for matrix operations and linear algebra solvers.
* [Intel® oneAPI MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html#gs.cq3i6h) : Free proprietary math library.

### Compiling (C++)
The first thing you need to know is that the MATLAB code is independent from C++. So no need to compile anything!

The C++ code depends on the Eigen library for matrix operations and Intel® oneApi MKL for the PARDISO sparse solvers and BLAS.

We provide a (hopefully) fast way of compiling everything. Jus run the shell script install.sh
This script will try to locate everything needed and compile the code with either icc or g++
If you run into trouble, edit the paths for Eigen and for your installation of oneAPI MKL.

### References
* Gabriel Moreira, Manuel Marques and João Paulo Costeira. Rotation Averaging in a Split Second: A Primal-Dual Method and a Closed-Form for Cycle Graphs. (To appear in) Internacional Conference on Computer Vision (ICCV), 2021. [PDF](https://arxiv.org/pdf/2109.08046.pdf)
* Gabriel Moreira, Manuel Marques and João Paulo Costeira. Fast Pose Graph Optimization via Krylov-Schur and Cholesky Factorization. IEEE/CVF Winter Conference on Applications of Computer Vision (WACV), 2021.

### Datasets

The datasets are provided as .g2o files. You will find certain datasets with a naming convention like "\_SO3.g2o". These were the ones used in our ICCV paper on Rotation Averaging. The difference between e.g., sphere.g2o and sphere_SO3.g2o is that the former was cleaned of repeated edges.

Some of the original datasets can be found [here](https://lucacarlone.mit.edu/datasets/):

* L. Carlone, R. Tron, K. Daniilidis, and F. Dellaert. Initialization Techniques for 3D SLAM: a Survey on Rotation Estimation and its Use in Pose Graph Optimization. In IEEE Intl. Conf. on Robotics and Automation (ICRA), pages 4597-4604, 2015.
* L. Carlone, D. M. Rosen, G. C. Calafiore, J. J. Leonard, and F. Dellaert. Lagrangian Duality in 3D SLAM: Verification Techniques and Optimal Solutions. In IEEE/RSJ Intl. Conf. on Intelligent Robots and Systems (IROS), 2015.

### Author
* Gabriel Moreira, October 2021

