/**
 * err.cpp
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

#include <string>
#include "err.hpp"

std::string mksGetStatusString(int status) {
    switch(status) {
        case mksMemAllocErr:
            return "Mks - Memory allocation error.";
            break;
        case mksEigenPardisoLDLErr:
            return "Mks Eigen Intel MKL Pardiso - Sparse LDL decomposition error.";
            break;
        case mksSuccess:
            return "Mks - No errors.";
            break;
        case mksPyramidInvalidOctave:
            return "Mks Pyramid - Invalid octave.";
            break;
        case mksPyramidInvalidScale:
            return "Mks Pyramid - Invalid scale.";
            break;
        case mksKrylovReorthogonalizationErr:
            return "Mks Krylov-Schur reorthogonalization error.";
            break;
        case mksKrylovInitialVectorErr:
            return "Mks Krylov-Schur initial vector error.";
            break;
        case mksIOErr:
            return "Mks IO error.";
            break;
        default:
            return "Mks - Unknown error";
    };
};
