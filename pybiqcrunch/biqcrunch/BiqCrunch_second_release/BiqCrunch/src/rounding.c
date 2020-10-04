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
#include <string.h>
#include <math.h>

#include "bb.h"
#include "biqcrunch.h"

extern int M;     // the number of eigenvalues
extern double *Z; // contains the eigenvectors


/* 
 * Rounds the fractional bound to an integer bound.
 *
 * Given a fractional bound that satisfies
 *
 *     subproblem value >= bound,
 *
 * we use the fact that the subproblem value is integer to provide 
 * a tighter integer bound that satisfies
 *
 *     subproblem value >= intbound.
 */
int roundBound(double bound, int max_problem) 
{
    int intbound = (max_problem) ? (long int) bound : (long int) ceil(bound);

    return intbound;
}


double primalHeuristic(Problem *P0, int *x) 
{
    int N = P0->n - 1;
    int i;
    double random;
    double gamma;
    int temp_x[N]; // temporary vector to store heuristic solutions
    double best;

    best = -BIG_NUMBER;

    for (gamma = 0.0; gamma < 1.0; gamma += 0.01) {

        for (i = 0; i < N; i++) {
            random = ((double) rand() / (double) RAND_MAX);

            if (random > gamma) {
                temp_x[i] = 1;
            } else {
                temp_x[i] = 0;
            }
        }    

        update_best(x, temp_x, &best, P0);
    }

    return best;
}


double sdpBoundHeuristic(Problem *P0, BobNode * node, int * x) 
{
    int N = P0->n - 1;
    int i;
    int count;
    int temp_x[N]; // temporary vector to store heuristic solutions
    double best;
    double cut_off;

    best = -BIG_NUMBER;

    // Using a deterministic cut-off
    for (count=0; count<N; count++) {
        cut_off = node->fracsol[count];

        for (i = 0; i < N; i++) {
            if (node->xfixed[i]) {
                temp_x[i] = node->sol.X[i];
            } else {
                if (node->fracsol[i] > cut_off) {
                    temp_x[i] = 1;
                } else {
                    temp_x[i] = 0;
                }
            }
        }

        update_best(x, temp_x, &best, P0);
    }

    // Using a different random cut-off for each xi
    for (count=0; count<100; count++) {

        for (i = 0; i < N; i++) {
            cut_off = ((double) rand() / (double) RAND_MAX); // var cut-off
            if (node->xfixed[i]) {
                temp_x[i] = node->sol.X[i];
            } else {
                if (node->fracsol[i] > cut_off) {
                    temp_x[i] = 1;
                } else {
                    temp_x[i] = 0;
                }
            }
        }

        update_best(x, temp_x, &best, P0);
    }

    return best;
}




double roundingHeuristic(Problem *P0, BobNode * node, int * x) 
{
    return sdpBoundHeuristic(P0, node, x);
}


/*
 * The Goemans/Williamson random hyperplane heuristic.
 */
double GWheuristic(P0, P, node, x, num)
    Problem *P0; // the original problem
    Problem *P;  // the current subproblem
    BobNode * node; // Node of the tree
    int *x;  // output binary vector of the solution obtained
    int num; // Number of random hyperplanes
{

    int i, j, count, index;
    int N = P0->n - 1;
    int subN = P->n;
    int temp_x1[N]; // temporary vector to store heuristic solutions
    int temp_x2[N]; // temporary vector to store heuristic solutions

    double sca;
    double best;
    double v[M];
    double vnorm, ivnorm;

    best = -BIG_NUMBER;

    // Goemans/Williamson random hyperplane rounding
    for (count = 0; count < num; count++)
    {
        // Compute the random hyperplane v
        vnorm = 0.0;
        for (i = 0; i < M; i++) 
        {
            v[i] = 1.0 + (int) 100.0 * rand() / ((double) RAND_MAX + 1);
            vnorm += v[i] * v[i];
        }
        vnorm = sqrt(vnorm);
        ivnorm = 1.0/vnorm;
        for (i = 0; i < M; i++)
            v[i] *= ivnorm;

        // Compute cuts temp_x1 and temp_x2 generated by hyperplane v
        index = 0;
        for (i = 0; i < N; i++)
        {
            if (node->xfixed[i])
            {
                temp_x1[i] = node->sol.X[i];
                temp_x2[i] = node->sol.X[i];
            }
            else
            {
                sca = 0.0;
                for (j = 0; j < M; j++)
                    sca += v[j] * Z[j * subN + index];

                if (sca < 0) 
                {
                    temp_x1[i] = 0;
                    temp_x2[i] = 1;
                }
                else
                {
                    temp_x1[i] = 1;
                    temp_x2[i] = 0;
                }

                index++;
            }
        }

        update_best(x, temp_x1, &best, P0);
        update_best(x, temp_x2, &best, P0);
    }

    return best;
}


/*
 * Given the current best solution, xbest, and a new solution, xnew, determines
 * the objective value and feasibility of xnew, then replaces xbest with xnew if
 * xnew is better. Also updates the best objective value, best.
 */
int update_best(int *xbest, int *xnew, double *best, Problem *P0) 
{
    int success = 0;
    int N = P0->n - 1;
    double heur_val;

    heur_val = BC_evaluateSolution(xnew);

    if ((*best < heur_val) && BC_isFeasibleSolution(xnew)) {
        memcpy(xbest, xnew, sizeof(int) * N);
        *best = heur_val;
        success = 1;
    }

    return success;
}


/*
 * Performs a simple local search starting from the given feasible solution x. 
 * Returns a feasible solution x that is locally optimal.
 * The objective value of the returned x will be stored in *val.
 */
int local_search(BobNode * node, int *x, double *val, Problem *P0) 
{
    int i;
    int success = 0;
    int local_opt_identified = 0;
    int N = P0->n - 1;
    int xtmp[N];

    // Make sure *val stores the objective value of x
    *val = BC_evaluateSolution(x);

    // Repeat until a local optimum has been identified
    while (!local_opt_identified) {
        // Initialize: xtmp = x
        memcpy(xtmp, x, sizeof(int) * N);

        for (i = 0; i < N; i++) {

            if (!node->xfixed[i]) {
                // Change xtmp[i] from 0 to 1, or from 1 to 0
                xtmp[i] = (xtmp[i] == 0) ? 1 : 0;

                // If xtmp is better than x, replace x with xtmp
                if (update_best(x, xtmp, val, P0)) {
                    success = 1;
                    break;
                }
                // Change xtmp[i] back
                xtmp[i] = (xtmp[i] == 0) ? 1 : 0;
            }

            // Test if a local optimum has been identified
            if (i == N-1) {
                local_opt_identified = 1;
            }
        }
    }

    return success;
}

