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

#include <math.h>
#include <stdio.h>
#include <sys/resource.h>

#include "bb.h"
#include "biqcrunch.h"

extern Inequality *Cuts; // vector of triangle inequality constraints
extern double *X; // Matrix to work on : modified by projSDP
extern int M; // the number of eigenvalues
extern double *W; // contains the eigenvalues
extern double *Z; // contains the eigenvectors
extern int *ISUPPZ; // dim = 2*max(1,M), support of the eigvecs. in Z
extern double *WORK; // dim = LWORK
extern int *IWORK; // dim = LIWORK
extern int sizeWORK, sizeIWORK;

void print_vector(double *vec, int N) 
{
    int i;
    for (i = 0; i < N; i++) {
        printf("%24.16e\n", vec[i]);
    }
}

void print_symmetric_matrix(double *Mat, int N) 
{
    int i, j;
    double val;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            val = (i >= j) ? Mat[i + j*N] : Mat[j + i*N];
            printf("%24.16e", val);
        }
        printf("\n");
    }
}

// computes gradient and function values for x
/*************** Simulator (evaluates function and gradient) ******************/
void sim(Problem *P, double *yy, double *ff, double *gg, double alpha) 
{
    int N = P->n;
    int N2 = N * N;

    // Let X = P->Q
    int incx = 1;
    int incy = 1;
    dcopy_(&N2, P->Q, &incx, X, &incy);

    // X = X - scaleEq * B^*(yy)
    B(TRANSP, P, yy, ff, gg, alpha);
    A(TRANSP, P, yy, ff, gg, alpha);

    // Overwrites X with X_+
    projSDP(ff, P); // returns *ff = 0.5*|X_+|_F^2

    double alphainv = 1.0 / alpha;
    *ff = *ff * 1000.0 * (alphainv / 1000.0); // now *ff = 0.5*|X_+|_F^2/alpha

    B(NOTRANSP, P, yy, ff, gg, alpha);
    A(NOTRANSP, P, yy, ff, gg, alpha);

    *ff += alpha * (0.5 * N2);
}


/***********************  projSDP  ************************/
void projSDP(double *ff, Problem *P) 
{
    int N = P->n;
    int j;
    double bornesup;
    double borneinf = 1e-8;
    double temp;

    char UPLO = 'L';
    int LDX = N;

    // Operator norm for X to get a upper bound of the spectral radius of X
    // |X|_2 \leq \sqrt{ |X|_1 |X|_inf }  (Holder's inequality)
    //          = |X|_1 = |X|_inf  (since X is symmetric)
    //
    // Frobenius norm is also an upper bound on the spectral radius of X:
    //      |X|_2 <= |X|_F
    char NORM = 'I';
    double norminf = dlansy_(&NORM, &UPLO, &N, X, &LDX, WORK);
    NORM = 'F';
    double normfro = dlansy_(&NORM, &UPLO, &N, X, &LDX, WORK);

    // bornesup = min(norminf, normfro)
    bornesup = (norminf < normfro) ? norminf : normfro;

    // Ensure that borneinf <= bornesup.
    if (bornesup < borneinf) {
        bornesup = 2.0 * borneinf;
    }

    //printf("\nX = \n"); print_symmetric_matrix(X, N);

    /* Compute the positive eigenvalues and associated eigenvectors of X.
     *
     * The M columns of Z will contain the orthonormal eigenvectors of the 
     * matrix X corresponding to the positive eigenvalues, the i-th column 
     * of Z holding the eigenvector associated with W[i].
     */
    char JOBZ = 'V';
    char RANGE = 'V';
    double VL = borneinf;
    double VU = bornesup;
    int IL = 0;
    int IU = 0;
    double ABSTOL = 1e-8;
    int LDZ = N;
    int LWORK = sizeWORK;
    int LIWORK = sizeIWORK;
    int INFO;
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, X, &LDX, &VL, &VU, &IL, &IU, &ABSTOL, 
            &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    //printf("M = %d\n", M); printf("W = \n"); print_vector(W, M);

    // Check if the eigensolver failed (i.e., if INFO != 0)
    if (INFO) {
        VL = 0.0;
        ABSTOL = 0.0;
        dsyevr_(&JOBZ, &RANGE, &UPLO, &N, X, &LDX, &VL, &VU, &IL, &IU, &ABSTOL, 
                &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
        if (INFO) {
            fprintf(stderr, 
                    "Error: eigenvalue computation failed (INFO = %d)\n", 
                    INFO);
            exit(1);
        }
    }

    // Compute ff = 0.5*||X_+||^2 = 0.5*||W||^2
    int INCW = 1;
    double normW = dnrm2_(&M, W, &INCW);
    *ff = 0.5*(normW*normW);

    // Compute Z = Z*Diag(W)^{1/2}
    int INCZ = 1;
    for (j = 0; j < M; j++) {
        // Scale jth column of Z by sqrt(W[j])
        temp = sqrt(W[j]);
        dscal_(&N, &temp, Z + j*N, &INCZ);
    }

    char TRANS = 'N';
    double ALPHA = 1.0;
    double BETA = 0.0;

    /* X = ALPHA*Z*Z^T + BETA*X = Z*Z^T
     * Only writes to lower-triangular part of X.
     * When M = 0, we will obtain X = 0.
     */
    dsyrk_(&UPLO, &TRANS, &N, &M, &ALPHA, Z, &LDZ, &BETA, X, &LDX);

    //printf("X_+ = \n"); print_symmetric_matrix(X, N);
}



/***************** A *********************/
void A(int mode, Problem *P, double *yy, double *ff, double *gg, double alpha)
    // Constraint:  A(X) <= a
    //
    // Evaluates:
    //      if mode == TRANSP:
    //          X = X + scaleIneq * A^*(yy)
    //
    //      if mode == NOTRANSP:
    //          ff = ff + scaleIneq * <a,yy>
    //          gg = scaleIneq*( a + A(X/alpha) )
    //
{
    double temp, alphainv = 1.0 / alpha, dd;
    int type, ii, jj, kk, ineq;

    int N = P->n, entry;

    if (mode == TRANSP) {

        for (ineq = 0; ineq < P->mA; ineq++) // for each equality constraint
        {
            temp = scaleEq * yy[P->mB + ineq];

            // Using sparse representation
            for (entry = 0; entry < P->As[ineq].nnz; entry++) // for each nonzero entry
            {
                ii = P->As[ineq].i[entry];
                jj = P->As[ineq].j[entry];
                dd = P->As[ineq].val[entry];
                if (ii == jj)
                    X[ii + jj * N] -= temp * dd;
                else {
                    X[ii + jj * N] -= temp * dd;
                    X[jj + ii * N] -= temp * dd;
                }
            }
        }

        for (ineq = 0; ineq < P->NIneq; ineq++) {

            type = Cuts[ineq].type;
            ii   = Cuts[ineq].i;
            jj   = Cuts[ineq].j;
            kk   = Cuts[ineq].k;

            temp = 0.5 * scaleIneq * yy[P->mB + P->mA + ineq];

            switch (type) {
                case 1:
                    X[ii + jj * N] += temp; //X[jj + ii * N] += temp;
                    X[ii + kk * N] += temp; //X[kk + ii * N] += temp;
                    X[jj + kk * N] += temp; //X[kk + jj * N] += temp;
                    break;
                case 2:
                    X[ii + jj * N] += temp; //X[jj + ii * N] += temp;
                    X[ii + kk * N] -= temp; //X[kk + ii * N] -= temp;
                    X[jj + kk * N] -= temp; //X[kk + jj * N] -= temp;
                    break;
                case 3:
                    X[ii + jj * N] -= temp; //X[jj + ii * N] -= temp;
                    X[ii + kk * N] += temp; //X[kk + ii * N] += temp;
                    X[jj + kk * N] -= temp; //X[kk + jj * N] -= temp;
                    break;
                default: //type == 4
                    X[ii + jj * N] -= temp; //X[jj + ii * N] -= temp;
                    X[ii + kk * N] -= temp; //X[kk + ii * N] -= temp;
                    X[jj + kk * N] += temp; //X[kk + jj * N] += temp;
            }
        }

    } else { // mode==NOTRANSP

        for (ineq = 0; ineq < P->mA; ineq++) // for each equality constraint
        {
            *ff += scaleEq * P->a[ineq] * yy[P->mB + ineq];

            gg[P->mB + ineq] = 0.0;
            for (entry = 0; entry < P->As[ineq].nnz; entry++) // for each nonzero entry
            {
                ii = P->As[ineq].i[entry];
                jj = P->As[ineq].j[entry];
                dd = P->As[ineq].val[entry];
                if (ii == jj)
                    gg[P->mB + ineq] -= dd * X[ii + jj * N];
                else
                    gg[P->mB + ineq] -= 2.0 * dd * X[ii + jj * N];
            }
            gg[P->mB + ineq] *= alphainv;
            gg[P->mB + ineq] += P->a[ineq];
            gg[P->mB + ineq] *= scaleEq;
        }

        for (ineq = 0; ineq < P->NIneq; ineq++) {

            type = Cuts[ineq].type;
            ii   = Cuts[ineq].i;
            jj   = Cuts[ineq].j;
            kk   = Cuts[ineq].k;

            *ff += scaleIneq * yy[P->mB + P->mA + ineq];

            switch (type) {
                case 1:
                    temp =  X[ii + jj * N] + X[ii + kk * N] + X[jj + kk * N];
                    break;
                case 2:
                    temp =  X[ii + jj * N] - X[ii + kk * N] - X[jj + kk * N];
                    break;
                case 3:
                    temp = -X[ii + jj * N] + X[ii + kk * N] - X[jj + kk * N];
                    break;
                default: // type == 4
                    temp = -X[ii + jj * N] - X[ii + kk * N] + X[jj + kk * N];
            }
            gg[P->mB + P->mA + ineq] = (temp * alphainv + 1.0) * scaleIneq;
        }

    } // end if

}




/***************** B *********************/
void B(int mode, Problem *P, double *yy, double *ff, double *gg, double alpha)
    // Constraint:  B(X) = b
    //
    // Evaluates:
    //      if mode == TRANSP:
    //          X = X - scaleEq * B^*(yy)
    //
    //      if mode == NOTRANSP:
    //          ff = ff + scaleEq * <b,yy>
    //          gg = scaleEq*( b - B(X/alpha) )
    //
{
    double alphainv = 1.0 / alpha;
    int ii, jj, eq, entry;
    double temp, dd;

    int N = P->n;

    if (mode == TRANSP) {
        for (eq = 0; eq < P->mB; eq++){ // for each equality constraint
            temp = scaleEq * yy[eq];

            // Using sparse representation
            for (entry = 0; entry < P->Bs[eq].nnz; entry++) // for each nonzero entry
            {
                ii = P->Bs[eq].i[entry];
                jj = P->Bs[eq].j[entry];
                dd = P->Bs[eq].val[entry];
                if (ii == jj)
                    X[ii + jj * N] -= temp * dd;
                else {
                    X[ii + jj * N] -= temp * dd;
                    X[jj + ii * N] -= temp * dd;
                }
            }
        }

    } else { // mode==NOTRANSP

        for (eq = 0; eq < P->mB; eq++) {// for each equality constraint

            *ff += scaleEq * P->b[eq] * yy[eq];

            gg[eq] = 0.0;
            for (entry = 0; entry < P->Bs[eq].nnz; entry++) // for each nonzero entry
            {
                ii = P->Bs[eq].i[entry];
                jj = P->Bs[eq].j[entry];
                dd = P->Bs[eq].val[entry];
                if (ii == jj)
                    gg[eq] -= dd * X[ii + jj * N];
                else
                    gg[eq] -= 2.0 * dd * X[ii + jj * N];
            }
            gg[eq] *= alphainv;
            gg[eq] += P->b[eq];
            gg[eq] *= scaleEq;
        }

    } // end if
}

