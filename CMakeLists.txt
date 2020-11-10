cmake_minimum_required(VERSION 3.10)

project(MAKS VERSION 1.0
        DESCRIPTION "Motion Averaging via Krylov-Schur"
        LANGUAGES CXX)

message("-- Intel Math Kernel Library (MKL) directory: ${MKL_ROOT}")
SET(MKL_INCLUDE "${MKL_ROOT}/include")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")

SET(COMPILE_FLAGS "-march=native -fopenmp -Ofast -DNDEBUG")

add_definitions(${COMPILE_FLAGS})

include_directories(include)
include_directories(${MKL_INCLUDE})
link_directories(${MKL_ROOT}/lib/intel64)

find_library(MKL_INTEL_LP64 mkl_intel_lp64 HINTS ${MKL_ROOT}/lib/intel64 REQUIRED)
find_library(MKL_GNU_THREAD mkl_gnu_thread HINTS ${MKL_ROOT}/lib/intel64 REQUIRED)
find_library(MKL_CORE mkl_core HINTS ${MKL_ROOT}/lib/intel64 REQUIRED)

file(GLOB SOURCES "src/*.cpp")

SET(MKL_LIBS -DMKL_LP64 -m64 -Wl,--start-group ${MKL_INTEL_LP64} ${MKL_GNU_THREAD} ${MKL_CORE} -Wl,--end-group -lgomp -lpthread -lm -ldl)

find_package(OpenMP REQUIRED)

add_executable(maks ${SOURCES})
target_link_libraries(maks ${MKL_LIBS})
