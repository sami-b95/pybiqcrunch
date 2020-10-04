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
// *                                                                           *
// *****************************************************************************
//                                                                             *
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
// Authors : N. Krislock, F. Roupin 2012, updated in 2016

#include <stdio.h>

#include "bb.h"
#include "biqcrunch.h"

extern BiqCrunchParameters params;

/**
 * Name of the heuristic printed in the output file
 */
char heur_name[] = "Max-Cut";


/**
 * Allocate the user defined global structure data (called at the beginning of the algorithm)
 * @param P0: original problem in domain [-1:1]
 * @param N: number of variables
 * @return 1 for success 0 for failure
 */
int BC_allocHeuristic(Problem *P0) {
    return 1;
}
/**
 * Free the user defined global structure data (called at the end of the algorithm)
 */
void BC_freeHeuristic() {
}

// User can fix other variables as soon as x[ic]=xic (i.e. 0 or 1), if contraints in problem allows it.
// Call by function Bob_GenChild (in bob_function.c)
// Recall that ic and xic are such that:  node->xfixed[ic] = 1; node->sol.X[ic] = xic;
void BC_FixVariables(BobNode *node, int ic, int xic) {
}

/**
 * Generic call of the heuristics.
 * This function is called at different times:
 * - at the beginning of the algorithm (PRIMAL_HEUR)
 * - during the evaluation of a node bound (SDP_BOUND_HEUR)
 * - after the evaluation of a node (ROUNDING_HEUR)
 * It has to return the value obtained during the execution of the chosen
 * heuristic and the binary vector of the solution
 *
 * @param P0: the original problem in domain [-1:1]
 * @param P: the current problem in domain [-1:1] (NULL if heuristic_code equals to PRIMAL_HEUR)
 * @param node: the current node (NULL if heuristic_code equals to PRIMAL_HEUR)
 * @param x: output binary vector of the solution obtained from the heuristic in domain [0:1]
 * @param heuristic_code: code that identify the timing of the call (PRIMAL_HEUR, SDP_BOUND_HEUR, ROUNDING_HEUR)
 * @return value of the solution found with the heuristic
 */
double BC_runHeuristic(Problem *P0, Problem *P, BobNode *node, int *x, int heuristic_code) {
    double heur_val;

    switch (heuristic_code) {

        case PRIMAL_HEUR:
            heur_val = primalHeuristic(P0, x);
            break;
        case SDP_BOUND_HEUR:
            //heur_val = sdpBoundHeuristic(P0, node, x);
            heur_val = GWheuristic(P0, P, node, x, params.NBGW1*P->n);
            break;
        case ROUNDING_HEUR:
            //heur_val = roundingHeuristic(P0, node, x);
            heur_val = GWheuristic(P0, P, node, x, params.NBGW2*P->n);
            break;
        default:
            printf("Choosen heuristic doesn't exist\n");
            exit(1);
    }

    // If x is feasible, perform a local search around x for a better solution
    if (params.local_search && BC_isFeasibleSolution(x)) {
        local_search(node, x, &heur_val, P0);
    }

    return heur_val;
}
