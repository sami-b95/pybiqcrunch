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
#include <math.h>

#include "bb.h"
#include "biqcrunch.h"

extern double *SubMat; // Construction matrices (size BobSizePb)
extern double *SubMatLin; // Construction vector (size BobSizePb)
extern Sparse SubMatS; //Construction matrices, sparse representation
extern double *tempMatrix; // Construction vector (size BobSizePb + 1)

extern BiqCrunchParameters params;
extern FILE *output;

/*
 * Evaluate a specific node.
 * This function computes the bound of a specific node
 * @param n: the node of the branch-and-bound tree
 * @return the bound of the node n
 */
double Evaluate(BobNode *node, Problem *SP, Problem *PP) 
{
    double bound;
    double best_sol;
    int i;
    int intbound;
    int temp_x[BobPbSize];


    Bob_STEVL(); // increment the number of evaluated nodes
    // Write node count, problem size, and fixed variables to output file
    fprintf(output, "\n");
    fprintf(output, 
            "**************************************************************************************\n");
    fprintf(output, 
            "                                     Node %d\n", Bob_EVL());
    fprintf(output, 
            "**************************************************************************************\n");
    fprintf(output, "Problem size %d\n", BobPbSize - countFixedVariables(node));
    if (BobSTAT1) {
	    fprintf(output, "Fixed variables:");
	    for (i = 0; i < BobPbSize; i++) {
	        if (node->xfixed[i]) {
	            fprintf(output, " (x[%d],%d)", i, node->sol.X[i]);
	        }
	    }
	    fprintf(output, "\n");
    }

    // Writes subproblem to PP
    createSubproblem(node, SP, PP);

    // Compute the SDP bound
    bound = SDPbound(node, SP, PP);

    // Checking int overflow: useful to test infeasibility
    if (bound >  BIG_NUMBER) bound =  BIG_NUMBER;
    if (bound < -BIG_NUMBER) bound = -BIG_NUMBER;

    //	Round the bound
    intbound = roundBound(bound, SP->max_problem);

    // Run Heuristic 3 if we are not able to prune
    if (Bob_ULBSup(intbound)) { // intbound < best_sol = Bob_ULBGet()

        for (i = 0; i < BobPbSize; i++) {
            if (node->xfixed[i]) {
                temp_x[i] = node->sol.X[i];
            }
            else {
                temp_x[i] = 0;
            }
        }

        if (params.heur_3 && !params.root) {
            BC_runHeuristic(SP, PP, node, temp_x, ROUNDING_HEUR);
        }

        updateSolution(node, temp_x, NEWSOL_HEUR3);
    }

    // Best solution found
    best_sol = (Bob_ULBGet() == BIG_NUMBER) ? INFINITY : Bob_ULBGet();

    // Save node information to the output file
    fprintf(output, "Depth = %d, Bound = %d, Best = %.0f\n", 
            node->level, 
            (SP->max_problem) ? -intbound : intbound, 
            (SP->max_problem) ? -best_sol : best_sol);
    fprintf(output, "\n");

    return bound;
}



// Count the number of fixed variables
int countFixedVariables(BobNode *node) 
{
    int i;
    int numFixedVariables = 0;

    for (i = 0; i < BobPbSize; i++) {
        if (node->xfixed[i]) {
            numFixedVariables++;
        }
    }

    return numFixedVariables;
}



/*
 * Writes subproblem to PP.
 *
 * Computes the subproblem removing the rows and the columns of the 
 * fixed variables from each matrix (Q, B[1..mB], A[1..mA])
 * 
 * SP is the original problem
 * PP is the subproblem (some variables are fixed)
 *
 * PP is made from SP
 */
void createSubproblem(BobNode *node, Problem *SP, Problem *PP) 
{
    // Subproblem size is the number of non-fixed variables in the node
    int SubN = BobPbSize - countFixedVariables(node);

    // Matrix size is SubN + 1
    PP->n = SubN + 1;

    /////////////////////////////////////////////////////////////
    // 1. Objective function of the subproblem: Q
    /////////////////////////////////////////////////////////////

    buildObjective(node, SP, PP);

    /////////////////////////////////////////////////////////////
    // 2. Equality constraints: B[1..mB]
    /////////////////////////////////////////////////////////////

    buildConstraints(node, SP->mB, SP->b, SP->Bs, PP->b, PP->Bs);

    /////////////////////////////////////////////////////////////
    // 3. Inequality constraints: A[1..mA]
    /////////////////////////////////////////////////////////////

    buildConstraints(node, SP->mA, SP->a, SP->As, PP->a, PP->As);
}


void buildObjective(BobNode *node, Problem *SP, Problem *PP) 
{
    int i, j;
    double Constant;
    int SubN = BobPbSize - countFixedVariables(node);

    // Build the quadratic part of the matrix Q
    getSubMatrix(SubMat, &(SP->Qs), node->xfixed, BobPbSize); 

    // Build the linear part of the matrix Q
    getLinear(SubMatLin, &(SP->Qs), node->sol.X, node->xfixed, BobPbSize); 

    // Build the constant part of the matrix Q
    Constant = SP->Q[(SP->n - 1) + (SP->n - 1) * SP->n]; 

    // Add linear terms to diagonal of submatrix
    for (i = 0; i < SubN; i++) {
        SubMat[i + i*SubN] += SubMatLin[i];
    }

    /* 
     * Build objective coefficient matrix PP->Q
     */

    // Constant part
    PP->Q[(PP->n - 1) + (PP->n - 1)*PP->n] = Constant;

    // Quadratic part
    for (i = 0; i < SubN; i++) {
        for (j = 0; j < SubN; j++) {
            if (i != j) {
                PP->Q[i + j*PP->n] = SubMat[i + j*SubN];
            } 
            else {
                PP->Q[i + j*PP->n] = 0.0;
            }
        }
    }

    // Linear part
    for (i = 0; i < SubN; i++) {
        PP->Q[i + (PP->n - 1)*PP->n] = SubMat[i + i*SubN]; 
        PP->Q[(PP->n - 1) + i*PP->n] = SubMat[i + i*SubN];
    }

    // Compute the sparse matrix from the dense one
    sparseMatrix(&(PP->Qs), PP->Q, PP->n); 

}



void buildConstraints(BobNode *node, int num, 
        double *RHSsource, Sparse *Matsource, 
        double *RHSdest,   Sparse *Matdest) 
{
    int i, k, s;
    int nz;
    double coeffMatNorm;
    double scale_factor;
    double Constant;
    int SubN = BobPbSize - countFixedVariables(node);

    for (s = 0; s < num; s++) {

        RHSdest[s] = RHSsource[s];

        coeffMatNorm = getSubMatrixSparse(&SubMatS, &(Matsource[s]), 
                node->xfixed, BobPbSize);

        //Build the linear part and the constant part
        Constant = getConstant(&(Matsource[s]), node->sol.X, node->xfixed, BobPbSize);
        getLinear(SubMatLin, &(Matsource[s]), node->sol.X, node->xfixed, BobPbSize);

        if (params.scaling) { // 1 / (1 + |a| + ||A||)
            scale_factor = 1.0 / (1.0 + fabs(RHSsource[s]) + coeffMatNorm); 
        }
        else {
            scale_factor = 1.0;
        }

        // Build the final submatrix
        nz = 0;

        // Quadratic part
        for (k = 0; k < SubMatS.nnz; k++) {
            if (fabs(SubMatS.val[k]) > 0.0) {
                Matdest[s].i[nz] = SubMatS.i[k];
                Matdest[s].j[nz] = SubMatS.j[k];
                Matdest[s].val[nz] = scale_factor*SubMatS.val[k];
                nz++;
            }
        }

        // Linear part
        for (i = 0; i < SubN; i++) {
            if (fabs(SubMatLin[i]) > 0.0) {
                Matdest[s].i[nz] = SubN; 
                Matdest[s].j[nz] = i;
                Matdest[s].val[nz] = scale_factor*SubMatLin[i];
                nz++;
            }
        }

        // new constant value on left-side
        if(fabs(Constant) > 0.0) {
            Matdest[s].i[nz] = SubN;   
            Matdest[s].j[nz] = SubN;
            Matdest[s].val[nz] = scale_factor*Constant;
            nz++;
        }

        // Final number of non-zero elements of the submatrix
        Matdest[s].nnz = nz;

        // Scale the right-side of the contraint
        RHSdest[s] *= scale_factor;
    }
}



/* 
 * Return the fixed value of the node.
 * The fixed value is contribution of the fixed variables to 
 * the objective value.
 */
double getFixedValue(BobNode *node, Problem *SP) 
{
    int i, j;
    int N = SP->n;
    double xi, xj;
    double fixedvalue = 0.0;

    // linear part (last column)
    j = N - 1;
    for (i = 0; i < BobPbSize; i++) {
        if (node->xfixed[i]) {
            xi = 2*node->sol.X[i] - 1.0;
            fixedvalue -= xi*(SP->Q[i + j*N] + SP->Q[j + i*N]);
        }
    }

    // quadratic part
    for (j = 0; j < BobPbSize; j++) {
        for (i = 0; i < BobPbSize; i++) {
            if (node->xfixed[i] && node->xfixed[j]) {
                xi = 2*node->sol.X[i] - 1.0;
                xj = 2*node->sol.X[j] - 1.0;
                fixedvalue -= xi*xj*SP->Q[i + j*N];
            }
        }
    }

    return fixedvalue;
}



double getSubMatrixSparse(Sparse *Mat, Sparse *sparseMat, int *fixed, int N) 
{
    double norm = 0.0;
    int i, k, ii, jj;
    int Nfixed = 0;
    int offset[N];
    int count;
    double val;

    count = 0;
    for (i = 0; i < N; i++) {
        if (fixed[i]) {
            offset[i] = 0;
            count++;
            Nfixed++;
        } 
        else {
            offset[i] = count;
        }
    }

    k = 0; //current number of non-zero elements of the submatrix Mat

    for (i = 0; i < sparseMat->nnz; i++) {
        if (sparseMat->i[i] < N && sparseMat->j[i] < N) {
            if (!fixed[sparseMat->i[i]] && !fixed[sparseMat->j[i]]) {
                ii = sparseMat->i[i] - offset[sparseMat->i[i]];
                jj = sparseMat->j[i] - offset[sparseMat->j[i]];
                val = sparseMat->val[i];
                if (jj<=ii) {
                    Mat->i[k] = ii;  //[ii + jj * Nfree] = val;
                    Mat->j[k] = jj;  //[jj + ii * Nfree] = val;
                }
                else {
                    Mat->i[k] = jj;  //[ii + jj * Nfree] = val;
                    Mat->j[k] = ii;  //[jj + ii * Nfree] = val;
                } 
                Mat->val[k] = val;
                k++; // Increase the number of non-zero elements of Mat
                norm += val*val;
                if (ii != jj) {
                    norm += val*val; // off-diagonal: counted twice!
                }
            }
        }
    }
    norm = sqrt(norm);

    Mat->nnz = k; //final number of non-zero elements of the submatrix

    return norm; // return the norm of the matrix
}



double getSubMatrix(double *Mat, Sparse *sparseMat, int *fixed, int N) 
{
    double norm = 0.0;
    int i;
    int Nfixed = 0;
    int Nfree;
    int offset[N];
    int count;
    int ii, jj;
    double val;

    count = 0;
    for (i = 0; i < N; i++) {
        if (fixed[i]) {
            offset[i] = 0;
            count++;
            Nfixed++;
        } 
        else {
            offset[i] = count;
        }
    }

    Nfree = N - Nfixed;

    // initialize the submatrix
    for (i = 0; i < Nfree * Nfree; i++) {
        Mat[i] = 0.0;
    }

    for (i = 0; i < sparseMat->nnz; i++) {
        if (sparseMat->i[i] < N && sparseMat->j[i] < N) {
            if (!fixed[sparseMat->i[i]] && !fixed[sparseMat->j[i]]) {
                ii = sparseMat->i[i] - offset[sparseMat->i[i]];
                jj = sparseMat->j[i] - offset[sparseMat->j[i]];
                val = sparseMat->val[i];
                Mat[ii + jj * Nfree] = val;
                Mat[jj + ii * Nfree] = val;
                norm += (val * val);
                if (ii != jj) {
                    norm += (val * val);
                }
            }
        }
    }
    norm = sqrt(norm);

    return norm; // return the norm of the matrix
}



double getConstant(Sparse *sparseMat, int *sol, int *fixed, int N) 
{
    double C = 0.0;
    int i;
    int ii, jj;
    double val;
    double tmp;

    for (i = 0; i < sparseMat->nnz; i++) {
        ii = sparseMat->i[i];
        jj = sparseMat->j[i];
        val = sparseMat->val[i];

        if (ii < N && jj < N) {
            if (fixed[ii] && fixed[jj]) {
                tmp = (2*sol[ii] - 1.0)*(2*sol[jj] - 1.0)*val;
                C += tmp;
                if (ii != jj) {
                    C += tmp;
                }
            }
        } 
        else if (ii == N && jj == N) {
            C += val;
        } 
        else if (jj == N) {
            if (fixed[ii]) {
                C += (2*sol[ii] - 1.0)*2*val;
            }
        } 
        else if (ii == N) {
            if (fixed[jj]) {
                C += (2*sol[jj] - 1.0)*2*val;
            }
        }
    }

    return C;
}



double getLinear(double *subMatLin, Sparse *sparseMat, int *sol, 
        int *fixed, int N)
{
    double sum = 0;
    int i;

    int Nfixed = 0;
    int Nfree;
    int offset[N];
    int count;
    int ii, jj;
    double val;

    count = 0;
    for (i = 0; i < N; i++) {
        if (fixed[i]) {
            offset[i] = 0;
            count++;
            Nfixed++;
        } 
        else {
            offset[i] = count;
        }
    }

    Nfree = N - Nfixed;

    // initialize the submatrix
    for (i = 0; i < Nfree; i++) {
        subMatLin[i] = 0.0;
    }

    for (i = 0; i < sparseMat->nnz; i++) {
        ii = sparseMat->i[i];
        jj = sparseMat->j[i];
        val = sparseMat->val[i];

        if (ii < N && !fixed[ii] && jj == N) {
            subMatLin[ii - offset[ii]] += val;
            sum += val*val;	
        } 
        else if (jj < N && !fixed[jj] && ii == N) {
            subMatLin[jj - offset[jj]] += val;
            sum += val*val;
        } 
        else if (jj < N && ii < N) {
            if (!fixed[ii] && fixed[jj]) {
                subMatLin[ii - offset[ii]] += (2*sol[jj] - 1.0)*val;
                sum += (2*sol[jj] - 1.0)*val*(2*sol[jj] - 1.0)*val;
            } 
            else if (!fixed[jj] && fixed[ii]) {
                subMatLin[jj - offset[jj]] += (2*sol[ii] - 1.0)*val;
                sum += (2*sol[ii] - 1.0)*val*(2*sol[ii] - 1.0)*val;
            }

        }

    }

    sum *= 2; // Linear terms will be off-diagonal => *2

    return sum;
}

