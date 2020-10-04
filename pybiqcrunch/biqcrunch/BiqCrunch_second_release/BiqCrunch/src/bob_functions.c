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
#include <sys/resource.h>

#include "bb.h"
#include "biqcrunch.h"

extern FILE *output;
extern BiqCrunchParameters params;
extern Problem *SP;
extern Problem *PP;

extern double TEMPI;  // initial time
extern int stopped;
extern double best_eval;
//extern double current_bound;


void processCommandLineArguments(int argc, char **argv) 
{
    // Control the command line arguments
    if (argc != 2) {
        fprintf(stderr, 
                "Usage: biqcrunch [-v (0|1)] file.bc file.param\n");
        exit(1);
    }

    // Create the output file
    if (!createOutputFile(argv[0])) {
        fprintf(stderr, 
                "Error: Cannot create output file in the current folder\n");
        exit(1);
    }

    // Read the input file instance
    Read_Data(argv[0]);

    // Read the parameters from a user file
    readParameter(argv[1]);
}


void initializeBobSolution() 
{
    BobSolution bs;
    int i;

    // Since Bob is solving a minimization problem, we initialize 
    // using the objective value of +INFINITY (i.e., BIG_NUMBER)
    for (i = 0; i < BobPbSize; i++) {
        bs.X[i] = 0;
    }
    Bob_ULBInit(BIG_NUMBER, &bs);
}


/*
 * Bob user function which initialize the problem, allocate the structures 
 * and evaluate the root node
 * @param argc: number of arguments on the command line
 * @param argv: array of the command line arguments
 */
void Bob_Init(int argc, char **argv) 
{
    // Start the timer
    TEMPI = temps_CPU();

    // Process the command line arguments
    processCommandLineArguments(argc, argv);

    // Seed the random number generator
    srand(params.seed);

    // Provide Bob with an initial solution
    initializeBobSolution();

    // Allocate the memory required
    allocMemory();

    // Redirect the Fortran output to /dev/null to have cleaner output
    redirect_output_();

    // Initialize the problem
    Init_UQP();
}


//=========Initialization : root node ==================//
void Init_UQP() 
{
    int i;
    double bound;
    int x[BobPbSize];
    int intbound;

    // Since Bob is minimizing, we must multiply objective values and bounds
    // by -1 when solving a maximization problem
    int sense = (SP->max_problem) ? -1 : 1;

    // Create the root node
    BobRoot = newNode(NULL);

    // Check if user value is provided
    if (params.soln_value_provided) {
        printf("Using solution value provided by user:  %d\n", 
                params.soln_value);
        Bob_ULBUpd(sense * params.soln_value, &(BobRoot->sol));
    }

    // run the heuristic to get initial solution
    for (i = 0; i < BobPbSize; i++) { x[i] = 0; } // initialize x = 0
//    if (params.heur_1 && !params.root) {
    if (params.heur_1) {
        BC_runHeuristic(SP, NULL, BobRoot, x, PRIMAL_HEUR);
    }
    updateSolution(BobRoot, x, NEWSOL_HEUR1);
    Bob_STNEXP();
    bound = Evaluate(BobRoot, SP, PP);  
    best_eval = bound;

    intbound = roundBound(bound, SP->max_problem);

    Bob_PRICREAT(BobRoot->Pri, intbound, 0, BobRoot->level);
}



/*
 * getBranchingVariable function used in the Bob_GenChild routine to determine
 * which variable x[ic] to branch on.
 *
 * @param Nodp:  the current node of the branch-and-bound search tree
 */
int getBranchingVariable(BobNode *Nodp) 
{
    int i;
    int ic = -1;  // x[ic] is the variable to branch on
    double maxValue, minValue;

    /* 
     * Choose the branching variable x[ic] based on params.branchingStrategy
     */
    if (params.branchingStrategy == LEAST_FRACTIONAL) {
        // Branch on the variable x[ic] that has the least fractional value
        maxValue = -BIG_NUMBER;
        for (i = 0; i < BobPbSize; i++) {
            if (!(Nodp->xfixed[i]) && fabs(0.5 - Nodp->fracsol[i]) > maxValue) {
                ic = i;
                maxValue = fabs(0.5 - Nodp->fracsol[ic]);
            }
        }
    }
    else if (params.branchingStrategy == MOST_FRACTIONAL) {
        // Branch on the variable x[ic] that has the most fractional value
        minValue = BIG_NUMBER;
        for (i = 0; i < BobPbSize; i++) {
            if (!(Nodp->xfixed[i]) && fabs(0.5 - Nodp->fracsol[i]) < minValue) {
                ic = i;
                minValue = fabs(0.5 - Nodp->fracsol[ic]);
            }
        }
    }
    else if (params.branchingStrategy == CLOSEST_TO_ONE) {
        // Branch on the variable x[ic] that is the closest to 1.0
        minValue = BIG_NUMBER;
        for (i = 0; i < BobPbSize; i++) {
            if (!(Nodp->xfixed[i]) && fabs(1.0 - Nodp->fracsol[i]) < minValue) {
                ic = i;
                minValue = fabs(1.0 - Nodp->fracsol[ic]);
            }
        }
    }
    else {
        fprintf(stderr, 
                "Error: Incorrect value for params.branchingStrategy\n");
        exit(1);
    }

    return ic;
}


/*
 * Determine if node is a leaf node by counting the number of fixed variables
 */
int isLeafNode(BobNode *node) 
{
    return (countFixedVariables(node) == BobPbSize);
}


/*
 * Create a new B&B node.
 *
 * If parentNode == NULL, it will create the root node.
 * Otherwise, the new node will be a child of parentNode.
 */
BobNode *newNode(BobNode *parentNode) 
{
    BobNode *node;
    int i;

    // allocate memory for the new child node
    node = Bob_NodeAlloc(0);

    if (node == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        exit(1);
    }

    // copy the solution information from the parent node
    for (i = 0; i < BobPbSize; i++) {
        if (parentNode == NULL) {
            node->xfixed[i] = 0;
            node->sol.X[i] = 0;
        }
        else {
            node->xfixed[i] = parentNode->xfixed[i];
            node->sol.X[i] = (node->xfixed[i]) ? parentNode->sol.X[i] : 0;
        }
    }

    // child is one level deeper than parent
    node->level = (parentNode == NULL) ? 0 : parentNode->level + 1;

    return node;
}


/* 
 * Bob user defined function that generate new children.
 *
 * The function generates two new children of the current node of the 
 * branch-and-bound tree.
 *
 * If params.branchingStrategy = LEASTFRACTIONAL, the branching is made 
 * on the variable that has the least fractional value.
 *
 * If params.branchingStrategy = MOSTFRACTIONAL, the branching is made 
 * on the variable that has the most fractional value.
 *
 * @param Nodp: current node of the branch-and-bound tree
 * @param Depth: depth of the node in the branch-and-bound tree
 * @param ExpCt:
 */
void Bob_GenChild(BobNode *Nodp, int Depth, BobTExpCt *ExpCt) 
{
    BobNode *node;
    double bound;
    int ic;  // x[ic] is the variable to branch on
    int xic;
    int intbound;

    // If the algorithm stops before finding the optimal solution, search in the 
    // nodes queue for the best bound (if the algorithm stops at the root node 
    // the queue is empty and the best bound is the root node bound)
    if (stopped || params.root || 
            (params.time_limit > 0 
             && (temps_CPU() - TEMPI) > params.time_limit) ) {
        int nb = Bob_PRIEVAL(Nodp->Pri);
        stopped = 1;
//        if (current_bound > nb) {
        if (best_eval > nb) {
            best_eval = nb;
//            current_bound = nb;
        }
        Bob_NodeFree(Nodp);
        return;
    }

    // If it's a leaf of the B&B tree, add the solution and don't branch
    if (isLeafNode(Nodp)) {
        updateSolution(Nodp, Nodp->sol.X, NEWSOL_LEAF);
        Bob_NodeFree(Nodp);
        return;
    }

    // Determine the variable x[ic] to branch on
    ic = getBranchingVariable(Nodp);
    fprintf(output, "Branching on x[%d] = %.2f\n", ic, Nodp->fracsol[ic]);

    // add two nodes to the search tree
    for (xic = 0; xic <= 1; xic++) { 

        // Create a new child node from the parent node, Nodp
        node = newNode(Nodp);

        // split on node ic
        node->xfixed[ic] = 1;
        node->sol.X[ic] = xic;
        fprintf(output, "Fixing x[%d] = %d\n", ic, xic);

        /***************************************************************** 
         * Fix values of other variables in node if possible (depending on 
         * problem). By default only X[ic] is fixed. 
         * Function BC_FixVariables(BobNode *node, int ic int xic) must be 
         * completed in the corresponding problem/xxx/heur.c
         *****************************************************************/
        BC_FixVariables(node, ic, xic);

	// If it's a leaf of the B&B tree, add the solution and don't branch
        if (isLeafNode(node)) {
    	    updateSolution(node, node->sol.X, NEWSOL_LEAF);
    	    Bob_NodeFree(node);
//	    return;
 	}
	else{

        /* 
         * Get the semidefinite-based bound for this node:
         *
         *     subproblem value >= bound
         */

        bound = Evaluate(node, SP, PP);
	Bob_STNEXP(); //increment the number of explored nodes
        /* 
         * Round the bound to an integer.
         *
         * Since the objective is always integer, we have:
         *
         *     subproblem value >= intbound
         */
        intbound = roundBound(bound, SP->max_problem);

        /* 
         * Define the priority of the node using the node bound, intbound,
         * and its depth, node->level.
         */
        Bob_PRICREAT(node->Pri, intbound, 0, node->level);

        /* 
         * if intbound < best_sol 
         * (i.e., global lower bound is greater than the node bound), 
         * then we must branch since there could be a better feasible 
         * solution in this subproblem
         */
        if (Bob_ULBSup(intbound)) {
            //fprintf(output, "Branching!\n");
            /*
             * if maximal depth reached, switch to depth-first search
             * otherwise insert node into the priority queue
             */
            if (Bob_ExpCtrl(Depth + 1, ExpCt)) {  
                Bob_GenChild(node, Depth + 1, ExpCt); // depth-first search
                Bob_NodeFree(node);
            } else {
                Bob_GPQIns(node); // insert node into the priority queue
            }
        } else { // otherwise, intbound >= beta, so we can prune
            //fprintf(output, "Prune!\n");
            Bob_NodeFree(node);
        }
	}
    } // end for xic

    Bob_NodeFree(Nodp);
}


/*
 * Bob function called at the end of the execution.
 * This function print the statistics of the execution, such as the number of 
 * nodes evaluated, the best solution found, the nodes best bound, the CPU time.
 * The Bob_end function is also used to free the memory allocated by the program.
 */
void Bob_End() 
{
    /*
     * Print results to the standard output and to the output file
     */

    printFinalOutput(stdout,Bob_EVL());
    printFinalOutput(output,Bob_EVL());

    // Free the memory
    freeMemory();

    // Close the output file
    closeOutputFile();
}


void printFinalOutput(FILE *file, int evalcount) 
{
    int sense = (SP->max_problem) ? -1 : 1;

    // Best solution found
    int best_sol = sense*Bob_ULBGet();

    // Best bound in the queue of the tree nodes
    double best_bound = sense*best_eval;

    // The gap between the best solution found and the bst bound
    double gap = 100.0 * 
        (best_sol != 0 ? (best_bound - best_sol)/best_sol : 1.0);

    fprintf(file, "Nodes = %d\n", evalcount);
    fprintf(file, "Root node bound = %.5lf\n", best_bound);

    if (!stopped) { // B&B not stopped early
        if (Bob_ULBGet() == BIG_NUMBER) {
            fprintf(file, "Problem is infeasible.\n");
        } else { 
            fprintf(file, "%s value = %d\n", 
                    (SP->max_problem) ? "Maximum" : "Minimum", best_sol);
            printSolution(file);
        }
    } else { // B&B stopped early
        if ((best_sol != BIG_NUMBER) && (best_sol != -BIG_NUMBER )) fprintf(file, "Best value = %d\n", best_sol);
	else fprintf(file, "Best value = inf\n");

	if (!params.root){
	fprintf(file, "Current bound = %.5lf\n",best_bound); 
//        fprintf(file, "%s bound = %.5lf\n", 
//                (params.root) ? "Root node" : "Current", best_bound);
	}
        fprintf(file, "Gap = %.2lf%%\n", fabs(gap));
    }

    fprintf(file, "CPU time = %.4f s\n", temps_CPU() - TEMPI);
}


void printSolution(FILE *file) 
{
    int i;

    fprintf(file, "Solution = { ");
    for (i = 0; i < BobPbSize; i++) {
        if (BobSol->X[i] == 1) {
            fprintf(file, "%d ", i + 1);
        }
    }
    fprintf(file, "}\n");
}


double temps_CPU() 
{
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_utime.tv_sec + (usage.ru_utime.tv_usec * 1e-6);
}
