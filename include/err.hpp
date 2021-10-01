/**
 * err.hpp
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

#ifndef ERR_HPP
#define ERR_HPP

#ifndef NDEBUG
#   define MKS_ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define MKS_ASSERT(condition, message) do { } while (false)
#endif

#include <string>
#include <iostream>

/* ERROR CODES */

typedef int MksStatus;

#define mksSuccess                       0          
#define mksHogSuccess                    0
#define mksMemAllocErr                  -1
#define mksEigenPardisoLDLErr           -2
#define mksOutsideWindow                -3
#define mksPyramidInvalidOctave         -4
#define mksPyramidInvalidScale          -5
#define mksKrylovReorthogonalizationErr -6
#define mksKrylovInitialVectorErr       -7
#define mksError                        -8
#define mksIOErr                        -9

std::string mksGetStatusString(int status);

#endif /* ERR_HPP */
