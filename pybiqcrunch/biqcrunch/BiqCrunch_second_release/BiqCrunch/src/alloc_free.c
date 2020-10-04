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

#include "bb.h"
#include "biqcrunch.h"
#include "global.h"


/*****************************************************************************/

void allocProblem(Problem * prob, int N, int mA, int mB, int max_problem) 
{
    int vec_size;
    prob->n = N;
    prob->mA = mA;
    prob->mB = mB;
    prob->max_problem = max_problem;

    vec_size = (prob->n * prob->n) / 2 + prob->n;

    // Q : Objective matrix
    alloc_matrix(prob->Q, prob->n, double);

    alloc_vector(prob->Qs.i, vec_size, int);
    alloc_vector(prob->Qs.j, vec_size, int);
    alloc_vector(prob->Qs.val, vec_size, double);

    // A : Inequality constraints
    //right-hand-side vector of inequality constraints
    alloc_vector(prob->As, prob->mA, Sparse);
    // list of inequality coefficient matrices
    alloc_vector(prob->a, prob->mA, double);

    // B : Equality constraints
    //right-hand-side vector of equality constraints
    alloc_vector(prob->Bs, prob->mB, Sparse);
    // list of inequality coefficient matrices
    alloc_vector(prob->b, prob->mB, double);
}


void allocCopyProblem(Problem *d, Problem *o) 
{
    int s;
    d->n = o->n;
    d->mA = o->mA;
    d->mB = o->mB;
    d->max_problem = o->max_problem;

    alloc_vector(d->As, d->mA, Sparse);
    alloc_vector(d->a, d->mA, double);

    alloc_vector(d->Bs, d->mB, Sparse);
    alloc_vector(d->b, d->mB, double);

    alloc_matrix(d->Q, d->n, double);

    d->Qs.nnz = o->Qs.nnz;

    alloc_vector(d->Qs.i, d->Qs.nnz, int);
    alloc_vector(d->Qs.j, d->Qs.nnz, int);
    alloc_vector(d->Qs.val, d->Qs.nnz, double);

    for (s = 0; s < d->mB; s++) {

        d->Bs[s].nnz = o->Bs[s].nnz;
        alloc_vector(d->Bs[s].i, d->Bs[s].nnz, int);
        alloc_vector(d->Bs[s].j, d->Bs[s].nnz, int);
        alloc_vector(d->Bs[s].val, d->Bs[s].nnz, double);
    }

    for (s = 0; s < d->mA; s++) {

        d->As[s].nnz = o->As[s].nnz;
        alloc_vector(d->As[s].i, d->As[s].nnz, int);
        alloc_vector(d->As[s].j, d->As[s].nnz, int);
        alloc_vector(d->As[s].val, d->As[s].nnz, double);
    }
}

void allocSDPbound(Problem *SP) 
{
    int N = SP->n;
    int wa_length;
    int iwa_length;
    int nmax;   // nmax = maximal size of problem to solve

    nmax = SP->mB + SP->mA + MaxNineqAdded; // nmax =  number of equality + inequality constraints + maximum number of cuts
    f = 0.0;

    alloc_matrix(X, N, double);
    alloc_vector(Cuts, MaxNineqAdded, Inequality);
    alloc_vector(List, params.cuts, Inequality);
    alloc_vector(g, nmax, double);
    alloc_vector(y, nmax, double);
    alloc_vector(binf, nmax, double);
    alloc_vector(bsup, nmax, double);
    alloc_vector(nbd, nmax, int);

    wa_length = (2 * mmax + 5) * nmax + 11 * mmax * mmax + 8 * mmax;
    iwa_length = 3 * nmax;

    alloc_vector(wa, wa_length, double);
    alloc_vector(iwa, iwa_length, int);

    allocProj(SP);
}


void allocProj(Problem *P0) 
{
    int N = P0->n;

    sizeWORK = 26 * N;
    sizeIWORK = 10 * N;
    alloc_vector(WORK, sizeWORK, double);
    alloc_vector(IWORK, sizeIWORK, int);
    alloc_vector(W, N, double);
    alloc_matrix(Z, N, double);
    alloc_vector(ISUPPZ, 2 * N, int);
}


void allocMemory() 
{
    int sparseVectorLength = (BobPbSize*(BobPbSize + 1))/2;

    // Allocate matrices and vector used in the evaluate function
    alloc_matrix(tempMatrix,  SP->n, double);
    alloc_matrix(SubMat,      BobPbSize, double);
    alloc_vector(SubMatLin,   BobPbSize, double);
    alloc_vector(SubMatS.i,   sparseVectorLength, int);
    alloc_vector(SubMatS.j,   sparseVectorLength, int);
    alloc_vector(SubMatS.val, sparseVectorLength, double);

    // Allocate matrices and vectors used to compute the bound
    allocSDPbound(SP);

    // Allocate the heuristic global structures
    if (!BC_allocHeuristic(SP)) {
        fprintf(stderr, 
                "Error: cannot allocate heuristic data structure\n");
        exit(1);
    }
}


/*
 * Free the global memory allocated and the memory allocated to compute 
 * the bound of the nodes of the branch-and-bound tree
 */
void freeMemory() 
{
    // free resources of the heuristic
    BC_freeHeuristic();

    // free problem instance
    freeProblem(SP);
    free(SP);
    freeProblem(PP);
    free(PP);

    // free resources used to compute the bound
    freeSDPbound();

    // free resources of BiqCrunch
    free(tempMatrix);
    free(BobSol);
    free(SubMat);
    free(SubMatLin);
}



void freeProblem(Problem *prob) 
{
    int i;
    free(prob->Q);
    free(prob->Qs.i);
    free(prob->Qs.j);
    free(prob->Qs.val);

    for (i = 0; i < prob->mA; i++) {

        free(prob->As[i].i);
        free(prob->As[i].j);
        free(prob->As[i].val);
    }

    free(prob->As);
    free(prob->a);

    for (i = 0; i < prob->mB; i++) {

        free(prob->Bs[i].i);
        free(prob->Bs[i].j);
        free(prob->Bs[i].val);
    }

    free(prob->b);
    free(prob->Bs);
}


void freeSDPbound() 
{
    freeProj();
    free(Cuts);
    free(List);
    free(X);
    free(g);
    free(y);
    free(binf);
    free(bsup);
    free(nbd);
    free(wa);
    free(iwa);
}


void freeProj() 
{
    free(W);
    free(Z);
    free(ISUPPZ);
    free(WORK);
    free(IWORK);
}
