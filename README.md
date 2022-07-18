## MAKS - Motion Averaging
MAKS contains routines for pose graph optimization and rotation averaging.

### Dependencies
* [Eigen](http://eigen.tuxfamily.org) : Header-only library used for matrix operations and linear algebra solvers.
* [Intel® oneAPI MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html#gs.cq3i6h) : Free proprietary math library.

### Compiling (C++)
The first thing you need to know is that the MATLAB code is independent from C++. So no need to compile anything!  

The C++ code depends on the Eigen library for matrix operations and Intel® oneApi MKL for the PARDISO sparse solvers and BLAS.

We provide a (hopefully) fast way of compiling everything. Just run.   
`source install.sh`  

This script will try to locate everything needed and compile the code with either `icc` or `g++`.    
If you run into trouble, edit the paths for Eigen and for your installation of oneAPI MKL in the shell script.  

The executable has the name `ravess_example`.  

To run:  
`./ravess_example ../data/sphere_SO3.g2o`  

### References
[1] Gabriel Moreira, Manuel Marques and João Paulo Costeira. Rotation Averaging in a Split Second: A Primal-Dual Method and a Closed-Form for Cycle Graphs.  _Internacional Conference on Computer Vision (ICCV)_, 2021. [[PDF](https://arxiv.org/pdf/2109.08046.pdf)] [[Video](https://youtu.be/256Mk1ywGjw)]   
[2] Gabriel Moreira, Manuel Marques and João Paulo Costeira. Fast Pose Graph Optimization via Krylov-Schur and Cholesky Factorization. _IEEE/CVF Winter Conference on Applications of Computer Vision (WACV)_, 2021.  [[Video #1](https://youtu.be/lsKUetY8wkA)] [[Video #2](https://youtu.be/HVk9iLAoeN4)] [[Video #3](https://youtu.be/_S-KZcDL5Nw)]  

### Datasets

The datasets are provided as .g2o files. You will find certain datasets with a naming convention like "\_SO3.g2o". These were the ones used in our ICCV paper [1]. The difference between e.g., sphere.g2o and sphere_SO3.g2o is that the latter was cleaned of repeated edges and has no vertex/node information, only edge measurements. 

Some of the original datasets can be found [here](https://lucacarlone.mit.edu/datasets/):  

[3] L. Carlone, R. Tron, K. Daniilidis, and F. Dellaert. Initialization Techniques for 3D SLAM: a Survey on Rotation Estimation and its Use in Pose Graph Optimization. In _IEEE Intl. Conf. on Robotics and Automation (ICRA)_, pages 4597-4604, 2015.  
[4] L. Carlone, D. M. Rosen, G. C. Calafiore, J. J. Leonard, and F. Dellaert. Lagrangian Duality in 3D SLAM: Verification Techniques and Optimal Solutions. In _IEEE/RSJ Intl. Conf. on Intelligent Robots and Systems (IROS)_, 2015.  


Gabriel Moreira
July 2022

