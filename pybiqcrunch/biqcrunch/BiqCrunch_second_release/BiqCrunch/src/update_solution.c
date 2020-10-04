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

extern Problem *SP;
extern FILE *output;

int updateSolution(BobNode *node, int *x, int credit) 
{
    int i;
    int solutionAdded = 0;
    double sol_value;
    BobSolution solx;

    if (BC_isFeasibleSolution(x)) 
    {
        // Copy x into solx
        for (i = 0; i < BobPbSize; i++) {
	      solx.X[i] = x[i];
        }

        /* If new solution is better than the global solution, 
         * then print the new solution.
         * Also, updating the global best solution removes all nodes 
         * from the queue with lower bounds larger than heur_val
         */
        // Since Bob minimizes, we must negate the BC objective value
        sol_value = -BC_evaluateSolution(x);  
//        if (Bob_ULBUpd((int) sol_value, &(node->sol))) {
          if (Bob_ULBUpd((int) sol_value, &solx)) {
            solutionAdded = 1;

            switch (credit) {
                case NEWSOL_LEAF:
                    fprintf(output, "Leaf node: \t");
                    break;
                case NEWSOL_HEUR1:
                    fprintf(output, "Heuristic 1: \t");
                    break;
                case NEWSOL_HEUR2:
                    fprintf(output, "Heuristic 2: \t");
                    break;
                case NEWSOL_HEUR3:
                    fprintf(output, "Heuristic 3: \t");
                    break;
            }

            Bob_PrintSolution();
        }
    }

    return solutionAdded;
}



/*
 * Print the value of the solution.
 *
 * @param bs: user defined structure that contains all the information 
 *            about the solution
 */
void Bob_PrintSolution() 
{
    int Value;

    Value = (SP->max_problem) ? -Bob_ULBGet() : Bob_ULBGet();

    fprintf(output, "Beta updated => %d\n", Value);
    printf("Node %d Feasible solution %d\n", Bob_EVL(), Value);
}


double BC_evaluateSolution(int *sol) 
{
    int k, ii, jj;
    double val;
    double heur_val = 0.0;

    for (k = 0; k < SP->Qs.nnz; k++) {
        ii  = SP->Qs.i[k];
        jj  = SP->Qs.j[k];
        val = SP->Qs.val[k];

        if (ii < SP->n - 1 && jj < SP->n - 1) {
            heur_val += 2*val*(2*sol[ii] - 1.0)*(2*sol[jj] - 1.0);
        } 
        else if (ii == SP->n - 1 && jj < SP->n - 1) {
            heur_val += 2*val*(2*sol[jj] - 1.0);
        } 
        else if (ii == SP->n - 1 && jj < SP->n - 1) {
            heur_val += 2*val*(2*sol[ii] - 1.0);
        } 
        else if (ii == SP->n - 1 && jj == SP->n - 1) {
            heur_val += val;
        }
    }

    return heur_val;
}



int BC_isFeasibleSolution(int *values) 
{
    int m, ii, jj, kk;
    int N = SP->n;
    double val;
    double sum;
    double tmp;
    int isfeas;
    int v[N];

    // Convert values (0,1) to (-1,1)
    for (ii = 0; ii < N - 1; ii++) {
        v[ii] = 2*values[ii] - 1;
    }
    v[N - 1] = 1;

    isfeas = 1;

    for (m = 0; m < SP->mB; m++) {
        sum = 0.0;

        for (kk = 0; kk < SP->Bs[m].nnz; kk++) {
            ii = SP->Bs[m].i[kk];
            jj = SP->Bs[m].j[kk];
            val = SP->Bs[m].val[kk];

            tmp = val*v[ii]*v[jj];
            sum += tmp;
            if (ii != jj) {
                sum += tmp;
            }
        }

        if (sum != SP->b[m]) {
            isfeas = 0;
            break;
        }
    }

    for (m = 0; m < SP->mA; m++) {
        sum = 0.0;

        for (kk = 0; kk < SP->As[m].nnz; kk++) {
            ii = SP->As[m].i[kk];
            jj = SP->As[m].j[kk];
            val = SP->As[m].val[kk];

            tmp = val*v[ii]*v[jj];
            sum += tmp;
            if (ii != jj) {
                sum += tmp;
            }
        }

        if (sum > SP->a[m]) {
            isfeas = 0;
            break;
        }
    }

    return isfeas;
}

