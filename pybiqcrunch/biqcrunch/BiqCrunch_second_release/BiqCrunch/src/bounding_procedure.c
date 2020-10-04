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
#include <string.h>

#include "bb.h"
#include "biqcrunch.h"

extern FILE *output;
extern BiqCrunchParameters params;

extern double TEMPI; // initial time
extern Inequality *Cuts; // vector of triangle inequality constraints
extern double *g, f; // gradient and function values
extern double *y; // dual variables
extern double *X; // Stores current X
extern double *binf; // lower bounds on the variables y
extern double *bsup; // upper bounds on the variables y
extern int *nbd; // indicates which variables are bounded
extern double *wa;  // Double workspace for L-BFGS-B
extern int *iwa; // Integer workspace for L-BFGS-B
extern double gradInorm; // norm of the gradient of f for equality constraints
extern double gradEnorm; // norm of the gradient of f for equality constraints
extern int M; // the number of eigenvalues
extern double *Z; // contains the eigenvectors


/********************************************************************/
/************************ SDPbound **********************************/
/********************************************************************/
double SDPbound(BobNode *node, Problem *SP, Problem *PP) 
{
    int i, index;
    int BFGSsuccess;
    double bound, gap;
    double oldf;
    int x[BobPbSize];
    double ynorm;
    int ylength;
    int incy = 1;
    double minAllIneq;
    int count; // major iteration counter
    int nbit; // number of iterations in BFGS

    // The fixed value is contribution of the fixed variables to 
    // the objective value
    double fixedvalue = getFixedValue(node, SP);

    if (PP->n == 1) {
        bound = fixedvalue - PP->Q[0];
        fprintf(output, 
                "N = 1: returning bound = %.1f\n", bound);
        return bound;
    }

    // Compute SDPbeta = fixedvalue - best_sol
    double beta = fixedvalue - Bob_ULBGet();

    // start with no cuts
    PP->NIneq = 0; 
    int NumAdded = 0;
    int NumSubtracted = 0;

    // Parameters
    int MinNiter = params.minNiter;

    // Set inital values for dynamically adjusted parameters
    double alpha = params.alpha0;
    double tol = params.tol0;
    int maxNAiter = params.maxNAiter;

// Recall that the scaling factors are:
//    scaleIneq = 1.0 / (1.0 + sqrt(3.0));
//    scaleEq = 1.0;

    // initial dual vector y = 0
    ylength = PP->mB + PP->mA + PP->NIneq;
    for (i = 0; i < ylength; i++) {
        y[i] = 0.0;
    }

    // Evaluate the function
    sim(PP, y, &f, g, alpha);

    count = 0;
    int done = 0;
    int stop_SDP_bound = 0;
    int giveup = 0;
    int prune = 0;
    int nbitalpha = 0; // current number of iterations for a given alpha

    // Main loop
    while (!done) {

        // Update iteration counter
        count++;
        nbitalpha++;
        oldf = f;

        // Call BFGS solver
        BFGSsuccess = 
            calllbfgsb(y, PP, alpha, tol, beta, fixedvalue, SP->max_problem, &nbit);

        bound = (SP->max_problem) ? f - fixedvalue : fixedvalue - f;

        prune = pruneTest(bound, SP->max_problem);

        // Compute the norm of the dual variables y
        ylength = PP->mB + PP->mA + PP->NIneq;
        ynorm = dnrm2_(&ylength, y, &incy);

        ////////////////////////////////////////////////////////////////////////
        //  Heuristic 2
        ////////////////////////////////////////////////////////////////////////
        if (params.heur_2 && !prune && !params.root) {

            for (i = 0; i < BobPbSize; i++) { x[i] = 0; }

            BC_runHeuristic(SP, PP, node, x, SDP_BOUND_HEUR);

            // Test if x is better than the best known solution
            updateSolution(node, x, NEWSOL_HEUR2);

            // Update beta
            beta = fixedvalue - Bob_ULBGet();

            // Update prune
            prune = pruneTest(bound, SP->max_problem);
        }
        ////////////////////////////////////////////////////////////////////////

        gap = (Bob_ULBGet() == BIG_NUMBER) ? INFINITY : f - beta;

        if (count >= MinNiter && !params.root) { 
            giveup = 
                (count >= params.maxNiter) || 
                (
                 (fabs(oldf - f) < 1.0) && 
                 (alpha < 1e-4) && 
                 (fabs(gap) > 2.0)
                );
        }

        // Add violated inequalities, remove satisfied inequalities
        if (params.withCuts && !prune && !giveup) {
            minAllIneq = updateInequalities(PP, y, &NumAdded, &NumSubtracted, count);
        }
        else {
            minAllIneq = NAN;
            NumAdded = 0;
            NumSubtracted = 0;
        }

        stop_SDP_bound = (alpha == params.minAlpha) && (NumAdded == 0);

        if (BobSTAT1) {
            if (count == 1) {
                fprintf(output,
                        "======================================================================================\n");
                fprintf(output, 
                        "%4s  %6s  %5s  %5s  %5s  %4s  %5s  %5s  %5s  %6s  %4s  %4s  %4s\n", 
                        "Iter", 
                        "Time", 
                        (params.root) ? "bound" : "gap", 
                        "alpha", 
                        "tol", 
                        "nbit", 
                        "Enorm", 
                        "Inorm", 
                        "ynorm", 
                        "minCut", 
                        "NCut", 
                        "NSub", 
                        "NAdd");
                fprintf(output, 
                        "======================================================================================\n");
            }
            fprintf(output, 
                    "%4d  %6.1f  %5.1f  %5.0e  %5.0e  %4d  %5.0e  %5.0e  %5.0e  %6.0e  %4d  -%-3d  +%-3d\n", 
                    count, 
                    temps_CPU() - TEMPI, 
                    (params.root) ? bound : gap, 
                    alpha, 
                    tol, 
                    nbit, 
                    gradEnorm, 
                    gradInorm, 
                    ynorm, 
                    minAllIneq, 
                    PP->NIneq, 
                    NumSubtracted, 
                    NumAdded);
        }
        fflush(NULL);

        // Test stopping criteria for adding inequalities
        done = 
            prune || // can prune the B&B tree 
            giveup || // see giveup tests above
            stop_SDP_bound || // see stop_SDP_bound tests above
            !BFGSsuccess || // BFGS failed
            (nbit >= params.nitermax); // reached max number of BFGS iters allowed

        // Decrease alpha and tol when # of inequals added is less than params.minCuts or if the # iterations for a given alpha is larger than maxNAiter (defined in biqcrunch.h)
        if (!done && ((NumAdded < params.minCuts) || nbitalpha>maxNAiter)) {
            nbitalpha = 0;
            alpha *= params.scaleAlpha;
            alpha = (alpha < params.minAlpha) ? params.minAlpha : alpha;
            tol *= params.scaleTol;
            tol = (tol < params.minTol) ? params.minTol : tol;
        }

        // Store the fractional solution in the node
        index = 0;
        for (i = 0; i < BobPbSize; i++) {
            if (node->xfixed[i]) {
                node->fracsol[i] = (double) node->sol.X[i];
            }
            else {
                // convert from -1..1 to 0..1
                node->fracsol[i] = 0.5*(X[(PP->n - 1) + index*PP->n] + 1.0); 
                index++;
            }
        }

    } // end main while loop

    fprintf(output, 
            "======================================================================================\n");
    if (prune) {
        fprintf(output, "Prune! ");
    }
    else if (giveup) {
        fprintf(output, "Giving up! ");
    }
    else if (stop_SDP_bound) {
        fprintf(output, "Stop bounding procedure! ");
    }
    else if (!BFGSsuccess) {
        fprintf(output, "BFGS failed! ");
    }
    else if (nbit >= params.nitermax) {
        fprintf(output, "nbit >= nitermax! ");
    }
    fprintf(output, 
            "Bound = %.2f, Iter = %d, alpha = %.0e, tol = %.0e, cuts = %d, time = %.1f\n",
            bound, count, alpha, tol, PP->NIneq, temps_CPU() - TEMPI);
    fprintf(output, 
            "======================================================================================\n");

    bound = fixedvalue - f;

    return bound;
}





/************************ calllbfgs *********************************/
int calllbfgsb(double *y, Problem *PP, double alpha, double tol, double beta, 
        double fixedvalue, int max_problem, int *nbit) 
{
    int i;
    char task[60], csave[60];
    int lsave[4], isave[44];
    double dsave[29];
    int status = 1;
    int N = PP->n;
    int iprint = -1;
    int ineq;
    int mem = mmax;
    int prune;
    int StopBFGS;
    int nvars;
    double mintemp;
    double factr; 
    double pgtol;
    double bound;

    // Compute the number of variables in BFGS
    nvars = PP->mB + PP->mA + PP->NIneq; 

    // Initialize y
    for (ineq = 0; ineq < PP->NIneq; ineq++) {
        y[PP->mB + PP->mA + ineq] = Cuts[ineq].y;
    }

    // evaluate f,g with updated y, alpha, and Inequalities (stored in vector Cuts)
    sim(PP, y, &f, g, alpha);

    /*
     * Compute the bounds for y
     */

    // equalities
    for (i = 0; i < PP->mB; i++) {
        nbd[i] = 0; // y[i] unbounded
    }

    //  inequalities
    for (i = PP->mB; i < nvars; i++) {
        nbd[i] = 1; // y[i] bounded from below
        binf[i] = 0.0; // lower bound is zero
    }

    /*
     * BFGS main loop
     */
    strcpy(task, "START"); for (i=5; i<60; i++) task[i] = ' ';
    StopBFGS = 0;
    *nbit = 0; // counts number of function/gradient evaluates
    while (!StopBFGS) {
        factr = 0.0; // CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
        pgtol = 0.0; // NORM OF PROJECTED GRADIENT <= PGTOL

        /* This calls the main L-BFGS-B function */
        setulb_(&nvars, &mem, y, binf, bsup, nbd, &f, g, 
                &factr, &pgtol, wa, iwa, task, &iprint, 
                csave, lsave, isave, dsave);

        if (strncmp(task, "FG", 2) == 0) { // L-BFGS-B requesting new f and g
            sim(PP, y, &f, g, alpha);
            (*nbit)++;
        } 
        else if (strncmp(task, "NEW_X", 5) == 0) { // L-BFGS-B found new x
            // test if we should stop
            // (e.g. if bound is worse than best current solution)

            // Compute Infinity-norm of gradE = g[0:(PP->mB)-1]
            gradEnorm = 0.0;
            for (i = 0; i < PP->mB; i++) {
                gradEnorm = (gradEnorm < fabs(g[i])) ? fabs(g[i]) : gradEnorm;
            }
            gradEnorm /= scaleEq; // remove scaling factor when computing norm

            // Compute Infinity-norm of gradI = min( g[(PP->mB):nvars-1], 0.0 )
            gradInorm = 0.0;
            for (i = PP->mB; i < nvars; i++) {
                mintemp = (g[i] > 0.0) ? 0.0 : g[i];
                gradInorm = (gradInorm < fabs(mintemp)) ? fabs(mintemp) : gradInorm;
            }
            gradInorm /= scaleIneq; // remove scaling factor when computing norm

            bound = (max_problem) ? f - fixedvalue : fixedvalue - f;

            prune = pruneTest(bound, max_problem);

            if ( prune // can prune the B&B tree
                    || (gradEnorm < tol && gradInorm < tol)
                    || ( *nbit >= params.nitermax ) ) {
                strcpy(task, "STOP"); for (i=4; i<60; i++) task[i] = ' ';
            }
        } 
        else { // L-BFGS-B has terminated
            StopBFGS = 1;
            if (strncmp(task, "STOP", 4) != 0) {
                //task[4] = '\0';
                //fprintf(output, "task = %s\n", task);
                status = 0;
            }
        }
    }


    // Scale: X = X/alpha
    double alphainv = 1.0 / alpha;
    int N2 = N * N;
    int incx = 1;
    dscal_(&N2, &alphainv, X, &incx);  // X = X/alpha


    // Update Z so that X = Z*Z'
    double sca = sqrt(alphainv);
    int NM = N * M;
    int incz = 1;
    dscal_(&NM, &sca, Z, &incz);  // Z = Z/sqrt(alpha)

    return status;

}



int pruneTest(double bound, int max_problem) 
{
    int prune;

    if (max_problem) {
        prune = ((long int) bound <= -Bob_ULBGet());
    }
    else {
        prune = (Bob_ULBGet() <= (long int) ceil(bound));
    }

    return prune;
}
