# install.sh
# This software and the related documents are provided as is, with no express or implied warranties, other than those that are expressly stated in the License.

echo "maks"
echo "This software and the related documents are provided as is, with no express or implied warranties, other than those that are expressly stated in the License."
echo "Version: October 2021"
echo "Compiling for x64 architecture with Intel MKL"
echo

# Look for a compiler: either ICC or G++
COMPILER_ARRAY=("icc" "g++")

for COMPILER in "${COMPILER_ARRAY[@]}"; do
    echo "Looking for $COMPILER..."

    if ! type "$COMPILER" &> /dev/null; then
        echo "$COMPILER not found."
    else
        echo "Found $COMPILER"
        CXX=$COMPILER
        break
    fi
done

if [ "$CXX" != "icc"  ] && [ "$CXX" != "g++"  ]; then
    echo "No compiler was found. Terminating"; return
else
    echo
    echo "Using the following compiler"
    eval "$CXX --version"
fi

# -> EDIT LINE BELOW WITH YOUR EIGEN3 PATH <-
EIGEN_INCLUDE="/usr/local/include/eigen3"

if [ -d $EIGEN_INCLUDE ]; then
	echo "Found Eigen3 in $EIGEN_INCLUDE"
    echo
else
    echo "User input required. Set the path to your Eigen3 root directory:"
    read EIGEN_INCLUDE
    if [ -d $EIGEN_INCLUDE ]; then
        echo "Found $EIGEN_INCLUDE"
        echo
    else
        echo "Could not find $EIGEN_INCLUDE. Terminating."
        #return
    fi
fi

# -> EDIT LINE BELOW WITH YOUR MKL /lib PATH <-
MKL_PATH_LIBS="/opt/intel/oneapi/mkl/latest/lib"

if [ -d $MKL_PATH_LIBS ]; then
	echo "Found Intel MKL libraries in $MKL_PATH_LIBS"
    echo
else
    echo "User input required. Set the path to your Intel MKL /lib directory:"
    read MKL_PATH_LIBS
    if [ -d $MKL_PATH_LIBS ]; then
        echo "Found Intel MKL libraries in $MKL_PATH_LIBS"
    else
        echo "Could not find $MKL_PATH_LIBS. Terminating."
        #return
    fi
fi

# -> EDIT LINE BELOW WITH YOUR MKL /include PATH <-
MKL_PATH_INCLUDE="/opt/intel/oneapi/mkl/latest/include"

if [ -d $MKL_PATH_INCLUDE ]; then
    echo "Found MKL header files directory $MKL_PATH_INCLUDE"
    echo
else
    echo "User input required. Set the path to your Intel MKL /include directory:"
    read MKL_PATH_INCLUDE
    if [ -d $MKL_PATH_INCLUDE ]; then
        echo "Found Intel MKL libraries in $MKL_PATH_INCLUDE"
    else
        echo "Could not find $MKL_PATH_INCLUDE. Terminating."
        #return
    fi
fi

# The actual MKL math stuff comes from libmkl_core
# For the MKL interface we need libmkl_intel_lp64
MKL_LIBS="-lmkl_intel_lp64 -lmkl_core"

# Runtime Threading Library (RTL) - An OpenMP implementation (libgomp in GNU or libiomp5 by Intel)
if [[ "$CXX" == "icc" ]]; then
    MKL_LIBS="$MKL_LIBS -liomp5"
elif [[ "$CXX" == "g++" ]]; then
    MKL_LIBS="$MKL_LIBS -lgomp"
fi

# The threading library from above will be used by the MKL threading library.
# It chooses according to the RTL from above
case ${CXX} in
    icc)    MKL_LIBS="$MKL_LIBS -lmkl_intel_thread";;
    g++)    MKL_LIBS="$MKL_LIBS -lmkl_gnu_thread";;
    *)      echo "${CXX} compiler not supported. Terminating"; return;;
esac

# Verbose
echo "MKL dynamic libraries:"
echo "$MKL_LIBS"
echo

# Optimization flags (CHANGE IF NEEDED ACCORDING TO YOUR COMPILER OF CHOICE)
case $CXX in
    icc)   OFLAGS="-O3 -march=native -ip -parallel -lstdc++ -std=c++11" ;;
    g++)   OFLAGS="-fopenmp -O3 -march=native -lstdc++ -std=c++11" ;;
    *)     echo "Compiler not supported. Terminating."; return;;
esac

# Verbose
echo "Optimization flags:"
echo "$OFLAGS"
echo

# Set preprocessor macros
echo "Preprocessor macros:"
MACROS="-DEIGEN_USE_MKL_ALL -DMKL_DIRECT_CALL -DNDEBUG"
echo "$MACROS"
echo

# Include path for header files
INCLUDE="./include"

# Source files
SRC="./src/ravess_example.cpp ./src/linalg.cpp ./src/ravg.cpp ./src/io.cpp"

echo "Creating executable..."

eval "$CXX $OFLAGS -o ravess_example $SRC -I$INCLUDE -I$EIGEN_INCLUDE -I$MKL_PATH_INCLUDE -L$MKL_PATH_LIBS -Wl,-rpath,$MKL_PATH_LIBS $MKL_LIBS -lpthread -lm -ldl $MACROS"

return

