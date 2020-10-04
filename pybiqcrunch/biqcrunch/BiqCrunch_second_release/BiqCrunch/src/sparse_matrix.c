// *****************************************************************************
// *                                                                           *
// *  BiqCrunch is a semidefinite-based solver for binary quadratic problems.  *
// *  It uses a branch-and-bound method featuring an improved semidefinite     *
// *  bounding procedure, mixed with a polyhedral approach. BiqCrunch uses     *
// *  particular input files format to describe the combinatorial problems.    *
// *                                                                           *
// *   Copyright (C) 2010-2016 Nathan Krislock, Jérôme Malick, Frédéric Roupin *
// *                                                                           *
// *                 http://www-lipn.univ-paris13.fr/BiqCrunch/                *
// *									       *
// *****************************************************************************
//									       *
//    This program is free software: you can redistribute it and/or modify     *
//    it under the terms of the GNU General Public License as published by     *
//    the Free Software Foundation, either version 3 of the License, or        *
//    (at your option) any later version.                                      *
//                                                                             *
//    This program is distributed in the hope that it will be useful,          *
//    but WITHOUT ANY WARRANTY; without even the implied warranty of           *
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
//    GNU General Public License for more details.                             *
//                                                                             *
//    You should have received a copy of the GNU General Public License        *
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.    *
//                                                                             *
// *****************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bb.h"
#include "biqcrunch.h"

void sparseMatrixA(int s, Problem *prob, Sparse *spMatrix, double *denseMatrix, int N) 
{
    int i, j, nz;
    double tmp;
    nz = 0;

    for (j = 0; j < N; j++) {
        for (i = j; i < N; i++) {
            tmp = denseMatrix[i + j * N];
            if (fabs(tmp) > 0.0) nz++;
        }
    }

    alloc_vector(prob->As[s].i, nz + 1, int);
    alloc_vector(prob->As[s].j, nz + 1, int);
    alloc_vector(prob->As[s].val, nz + 1, double);

    sparseMatrixScaled(spMatrix, denseMatrix, N, 1.0);
}

void sparseMatrixB(int s, Problem *prob, Sparse *spMatrix, double *denseMatrix, int N) 
{
    int i, j, nz;
    double tmp;
    nz = 0;

    for (j = 0; j < N; j++) {
        for (i = j; i < N; i++) {
            tmp = denseMatrix[i + j * N];
            if (fabs(tmp) > 0.0) nz++;
        }
    }

    alloc_vector(prob->Bs[s].i, nz + 1, int);
    alloc_vector(prob->Bs[s].j, nz + 1, int);
    alloc_vector(prob->Bs[s].val, nz + 1, double);

    sparseMatrixScaled(spMatrix, denseMatrix, N, 1.0);
}

void sparseMatrixScaled(Sparse *spMatrix, double *denseMatrix, int N, double scale) 
{
    int i, j, nz;
    double tmp;

    nz = 0;
    for (j = 0; j < N; j++) {
        for (i = j; i < N; i++) {
            tmp = denseMatrix[i + j * N];
            if (fabs(tmp) > 0.0) {
                spMatrix->i[nz] = i; //invert the values so j is always <= i
                spMatrix->j[nz] = j;
                spMatrix->val[nz] = tmp * scale;
                nz++;
            }
        }
    }

    spMatrix->nnz = nz;
}

void sparseMatrix(Sparse *spMatrix, double *denseMatrix, int N) 
{
    sparseMatrixScaled(spMatrix, denseMatrix, N, 1.0);
}

