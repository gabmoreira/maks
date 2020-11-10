# MAKS
### Motion Averaging via Krylov-Schur

# Dependencies
MAKS is built on top of:
* [Eigen](http://eigen.tuxfamily.org) : Header-only library used for matrix operations and sparse solvers.
* [Intel® MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html) : Free proprietary math library.


# Compilation

Get Intel® MKL for Linux [here](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library/choose-download/linux.html) and install it.

Clone the repository:  
`git clone https://github.com/gabmoreira/pipe.git`  

Create a build directory 
`mkdir build; cd build;`

Run cmake and set -DMKL_ROOT to the root directory of your MKL installation e.g., ~/intel/mkl
`cmake -DMKL_ROOT=<MKL_ROOT> ..`

Make everything
`make; cd ..`

Head over to /bin and enjoy!  
`./maks ../data/parking-garage.g2o`

# Author
* Gabriel Moreira

