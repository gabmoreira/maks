# MAKS
### Motion Averaging Krylov-Schur

## Dependencies
* [Eigen](http://eigen.tuxfamily.org) : Header-only library used for matrix operations and sparse solvers.
* [Intel® MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html) : Free proprietary math library.


## Compiling
Get Intel® MKL for Linux [here](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library/choose-download/linux.html) and install it.

Clone the repository  
`git clone https://github.com/gabmoreira/maks.git`  

Create a build directory  
`mkdir build; cd build`

Run cmake and set `-DMKL_ROOT` to the root directory of your MKL installation <MKL_ROOT>  
`cmake -DMKL_ROOT=<MKL_ROOT> ..`

Your <MKL_ROOT> will look something like `<...>/intel/mkl` 

Compile everything and move back to the root directory  
`make; cd ..`

## Usage
Head over to `bin`. Use the first argument to specify a `.g2o` file. Example datasets are provided in `/data`  
`./maks ../data/grid3D.g2o`

## References
* Gabriel Moreira, Manuel Marques and João Paulo Costeira. Fast Pose Graph Optimization via Krylov-Schur and Cholesky Factorization. To appear in Winter Conference on Applications of Computer Vision (WACV), 2021.

## Datasets
* L. Carlone, R. Tron, K. Daniilidis, and F. Dellaert. Initialization Techniques for 3D SLAM: a Survey on Rotation Estimation and its Use in Pose Graph Optimization. In IEEE Intl. Conf. on Robotics and Automation (ICRA), pages 4597-4604, 2015.
* L. Carlone, D. M. Rosen, G. C. Calafiore, J. J. Leonard, and F. Dellaert. Lagrangian Duality in 3D SLAM: Verification Techniques and Optimal Solutions. In IEEE/RSJ Intl. Conf. on Intelligent Robots and Systems (IROS), 2015.


## Author
* Gabriel Moreira, 2021

