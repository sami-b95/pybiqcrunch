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
//
// Author : Geoffrey Kozak
// 2012
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "bb.h"
#include "biqcrunch.h"

////////////////////////////// Global variables //////////////////////////////////
extern BiqCrunchParameters params;

/**
 * Name of the heuristic printed in the output file
 */
char heur_name[] = "Max-independent-set";


//Adjacency sparse matrix of the problem.
Sparse* M_adj;

//An order in which we visit the node (used in the different heuristic methods).
int* orderNode;

//The weights of each nodes.
double* weights;

//////////////////////////////////////////////////////////////////////////////////

////////////////////////// Function declarations /////////////////////////////////

// Compute the adjacency matrix of the graph of the input problem.
void getSparseAdjacencyMatrix(Problem *P0);

// Get the value of (i,j) in a sparse matrix.
double getSparseAdjacencyMatrixValue(Sparse* mAdj, int i, int j);

//Get the weights of each nodes and store them in a table.
void getNodeWeights(Problem *P0);

// Sort the nodes by increasing degree.
void orderNodeByDegree(Sparse* mAdj, int N);

// Sort the nodes by decreasing fractional solution.
void orderNodeByFracSol(Problem *P0, BobNode *node);

//The heuristic calls at the beginning of the program.
double first_stable_heurist(Problem *P0, int *x);

// Heuristic 1.
/* Description :
   1) Get the current solution from BiqCrunch.
   2) Sort the nodes by increasing degrees.
   3) Following the previous order computed, we try to
   add one by one each nodes.
   If he has no neighboor in the independent set we add
   it in the current independent set, else we don't add it.
   */
double bab_stable_heuristic(Problem * P0, BobNode *node, int *x);

// Heuristic 2.
/* Description :
   1) Get the current solution from BiqCrunch.
   2) Sort the nodes using the datas given by BiqCrunch.
   Each nodes has a fractional value. So we will sort our
   nodes by decreasing fractional value
   3) Following the previous order computed, we try to
   add one by one each nodes.
   If he has no neighboor in the independent set we add
   it in the current independent set, else we don't add it.
   */
double bab_stable_heuristic2(Problem * P0, BobNode *node, int *x);

// Heuristic 3.
/* Description :
   1) Get the current solution from BiqCrunch.
   2) Depending on a value called seuil, we will fixed
   our variables at 1 if its fractional value is higher
   and 0 if it is lower. (Note : we will probably have
   a non-admissible solution)
   3) We will correct the previous computed solution
   a) Build the adjacency matrix of the value fixed to
   one.
   b) Sort only the nodes we have fixed to one (we don't
   take the variable fixed by BiqCrunch) by decreasing degree.
   c) Following the previous order compute, we take the first
   node of the order (the highest degrees) we remove it from
   the current independent set. If it still exist edge in
   the current independent set problem back to step 3-a, else
   it's finished.
   */
double bab_stable_heuristic3(Problem * P0, BobNode *node, int *x);

//////////////////////////////////////////////////////////////////////////////////

////////////////////////// Function definitions //////////////////////////////////

/** Allocate memory for all fixed problem datas.
 * This function permit to compute previous datas which will be useful during
 * the heuristic calls. For some datsa we don't need to compute at each heuristic
 * calls so we compute the datas at the beginning of the program and we can use it
 * during the execution.
 * @param P0
 * @param N
 * @return a "boolean" value, 1 if the function succeed and 0 else
 */
int BC_allocHeuristic(Problem *P0) {

    int N = P0->n - 1;

    // Get the weight of each nodes and put them in the weights table
    getNodeWeights(P0);

    // Compute a sparse version of adjacency matrix
    getSparseAdjacencyMatrix(P0);

    // Sort node by increasing degree and put them in the table orderNode
    orderNodeByDegree(M_adj, N);

    return 1;
}

/** Free memory occupied by all fixed problem datas.
 * This function permit to free the space occupied by the datas computed
 * at the beginning of the program.
 * @see alloc_heuristic()
 */
void BC_freeHeuristic() {

    // Free the adjacency sparse matrix
    free(M_adj->i);
    free(M_adj->j);
    free(M_adj->val);
    free(M_adj);

    // Free the order node table
    free(orderNode);

    // Free the weights table
    free(weights);
}

// User can fix other variables as soon as x[ic]=xic (i.e. 0 or 1), if contraints in problem allows it.
// Call by function Bob_GenChild (in bob_function.c)
// Recall that ic and xic are such that:  node->xfixed[ic] = 1; node->sol.X[ic] = xic;
void BC_FixVariables(BobNode *node, int ic, int xic) {
    int j;

    if (xic == 1) {
        for (j = 0; j < BobPbSize; j++) {
            if (j != ic)
                if (getSparseAdjacencyMatrixValue(M_adj, ic, j) == 1.) {
                    // (ic,j) belongs to E so X[ic] = 1 => X[j] = 0
                    if(!node->xfixed[j]) {
                        node->xfixed[j] = 1;
                        node->sol.X[j] = 0;
                    }
                }
        }
    }
}


/** Run a heuristic depending on a heuristic code.
 * We can say that this function is THE MAIN function of this file.
 * Indeed, it will be this function which will be called when we want
 * call a heuristic method. Depending on the heuristic code given in
 * input, we will run the first heuristic, the sdp bound heuristic (during evaluation)
 * or the rounding heuristic (call for each node af the branch and bound tree).
 * @param P0 the initial Problem.
 * @param P the current Problem.
 * @param node the node of the branch and bound tree.
 * @param x the output integer solution.
 * @param heuristic_code the code of the heuristic method we want to call (PRIMAL_HEUR/SDP_BOUND_HEUR/ROUNDING_HEUR).
 * @return a double value, the heuristic value of the heuristic method called.
 */
double BC_runHeuristic(Problem *P0, Problem *P, BobNode *node, int *x, int heuristic_code)
{
    int N = P0->n - 1; // problem size
    int temp_x[N]; // temporary vector to store heuristic solutions
    double heur_val = 0; // value ofr a given heuristic
    double best; // best value 

    best = heur_val;
    switch (heuristic_code) {
        case PRIMAL_HEUR:
            // Specific MIS heuristic
            heur_val = first_stable_heurist(P0, x);
            update_best(x, temp_x, &best, P0);
            break;

        case SDP_BOUND_HEUR:
            // Specific MIS heuristic
            heur_val = bab_stable_heuristic2(P0, node, temp_x);
            update_best(x, temp_x, &best, P0);
            // GW heuristic
            heur_val = GWheuristic(P0, P, node, temp_x, params.NBGW1*P->n);
            update_best(x, temp_x, &best, P0);
            break;

        case ROUNDING_HEUR:
            // Specific MIS heuristic
            heur_val = bab_stable_heuristic2(P0, node, temp_x);
            update_best(x, temp_x, &best, P0);
            // GW heuristic
            heur_val = GWheuristic(P0, P, node, temp_x, params.NBGW2*P->n);
            update_best(x, temp_x, &best, P0);
            break;

    }

    // If x is feasible, perform a local search around x for a better solution
    if (params.local_search && BC_isFeasibleSolution(x)) {
        local_search(node, x, &best, P0);
    }

    return best;

}

/** Get the weights of each nodes and store them in a table.
 * This function get the weight of each node from the problem given in input
 * and store them in the global table of the weights called weights.
 * Note : This function is called once at the beginning of the program.
 * @param P0 a Problem.
 * @see alloc_heuristic()
 */
void getNodeWeights(Problem *P0) {
    int i;

    weights = (double*) malloc((P0->n - 1) * sizeof(double));

    for (i = 0; i < P0->n - 1; i++)
        weights[i] = 4 * P0->Q[i * P0->n + (P0->n - 1)];
}

/** Compute the adjacency matrix of the graph of the input problem.
 * This function compute a sparse adjacency matrix of the input problem graph
 * using the sparse equality constraints matrices and store it in the global
 * sparse variable M_adj.
 * Note : This function is called once at the beginning of the program.
 * @param P0 a Problem.
 */
void getSparseAdjacencyMatrix(Problem *P0) {
    int N = P0->n - 1, i, k, nz = 0; //j

    // The sparse version of the adjacency matrix is declared as a global variable
    M_adj = (Sparse *) malloc(sizeof(Sparse)); // Dynamic allocation of the sparse matrix
    M_adj->i = (int *) malloc(N * N * sizeof(int)); // Dynamic allocation of the table of indice i
    M_adj->j = (int *) malloc(N * N * sizeof(int)); // Dynamic allocation of the table of indice j
    M_adj->val = (double *) malloc(N * N * sizeof(double)); // Dynamic allocation of the table of matrix value (i,j) = val

    // Computing adjacency sparse matrix working on sparse equality constraints matrix
    // Complexity o(k*n) [k: number of equality constraint and n : number of node]
    for (i = 0; i < P0->n; i++) // For each node
        for (k = 0; k < P0->mB; k++) // For each equalities matrix constraint
        {
            if (P0->Bs[k].nnz == 4) // The equality matrix constraint describing the edge of the graph contains
            { // only 4 values non equal to zero (8 if we consider the entire matrix)

                if (P0->Bs[k].val[0] == 0.125) // Only the first value is useful indeed when we construct the sparse matrix
                { // information on the edge is the first added and O.125 because we transform
                    // the dense matrix (by multipliating it by UT and U) which is the reference for
                    // our sparse matrix

                    if (P0->Bs[k].j[0] == i) // If the indice j in the equality matrix constraint is equal to i
                    {
                        M_adj->i[nz] = P0->Bs[k].j[0]; // then we specify that between i in M_adj = j in Bs
                        M_adj->j[nz] = P0->Bs[k].i[0]; // and j in M_adj = i in Bs
                        M_adj->val[nz] = 1.; // there is an edge so we put 1 for (j,i)
                        nz++; // increase the number of non zero in sparse adjacency matrix
                    } else if (P0->Bs[k].i[0] == i) // If the indice i in the equality matrix constraint is equal to i
                    {
                        M_adj->i[nz] = P0->Bs[k].i[0]; // then we specify that between i in M_adj = i in Bs
                        M_adj->j[nz] = P0->Bs[k].j[0]; // and j in M_adj = j in Bs
                        M_adj->val[nz] = 1.; // there is an edge so we put 1 for M_adj(i,j)
                        nz++; // increase the number of non zero in sparse adjacency matrix
                    }
                }
            }
        }

    M_adj->nnz = nz; // we set the number of non zero in sparse adjacency matrix
}

/** Get the value of (i,j) in a sparse matrix.
 * This function permit to get the (i,j) value of a sparse matrix.
 * @param mAdj a sparse matrix.
 * @param i index i of the value we want.
 * @param j index j of the value we want.
 * @return a double, the value mAdj(i,j) or -1.0 if the value (i,j) isn't in the sparse matrix.
 */
double getSparseAdjacencyMatrixValue(Sparse* mAdj, int i, int j) {
    int k = 0;

    while ((mAdj->i[k] != i || mAdj->j[k] != j) && k < mAdj->nnz) //we looking for the (i,j) index k in the sparse matrix
        k++;

    if (k >= mAdj->nnz) //if we check all the sparse matrix value
        return -1.0; //return -1.0 to specify that (i,j) value doesn't exist in the sparse matrix

    return mAdj->val[k]; //return (i,j) value of the sparse matrix
}

/** Get the minimum index of the input table.
 *
 * @param t an integer table.
 * @param N the size of the input table.
 * @param max the maximum value of the table.
 * @return an integer, the index which has the smallest value in the input table.
 */
int getIndMin(int* t, int N, int max) { // Because we'll call this function several times until we'll deal with all the value in this table
    int i, indMin = 0, min = max; // we'll consider if a value in the table is lower than 0 we don't consider it (indeed we'll use this fonction
    // on a node degree table so we can do that because a degree is always higher than or equal to 0)

    for (i = 0; i < N; i++) // for each value in the table
    {
        if (t[i] < min && t[i] >= 0) // if the current value is lower than the minimum and this value is higher than or equal to 0
        {
            indMin = i; // we modify the minimum indice
            min = t[i]; // and also the minimum value
        }
    }
    return indMin; // we return the minimum indice computed
}

/** Sort the nodes by increasing degree.
 * This function compute an order. In fact we sort the nodes
 * by increasing degrees, we put the ouput order in the global
 * variable orderNode.
 * @param mAdj an adjacency SPARSE matrix.
 * @param N the number of variables (nodes).
 */
void orderNodeByDegree(Sparse* mAdj, int N) // Sort the node by increasing degree
{
    int* nodeDeg = (int*) malloc(N * sizeof(int)); // Table of node degree
    orderNode = (int*) malloc(N * sizeof(int)); // Table of node sort by increasing degree
    int i = 0, n = 0, nbnode = 0, degmax = 0;

    //Initialization
    for (i = 0; i < N; i++) {
        orderNode[i] = -1;
        nodeDeg[i] = 0;
    }

    for (i = 0; i < mAdj->nnz; i++) // For each node (indeed our computing sparse adacency matrix function gather all the edge of the
        // same head and this matrix is sort by lexicographic order 1 to n)
    {
        degmax = 0; // maximum degree of the graph
        if (mAdj->i[i] == n) // if it the same head (same node)
            nbnode++; // then increase degree of the node
        else // else it is a new head (node)
        {
            nodeDeg[n] = nbnode; // then we save the degree of the current node
            n = mAdj->i[i]; // we change the current node
            nbnode = 1; // we re-initialize the degree of node to 1
        }
        if (n > degmax) // If degree computed is higher than maximum
            degmax = n; // we change maximum value
    }
    nodeDeg[n] = nbnode; // For the last node, we don't forget to save its degree

    for (i = 0; i < N; i++) // Compute the table of node sort by increasing degree
    {
        n = getIndMin(nodeDeg, N, degmax); // Compute the indice of the minimum of table of node degree
        nodeDeg[n] = -1; // We specify that we have placed the node in the table
        orderNode[i] = n; // put the node with the minimum degree in the table
    }

    free(nodeDeg); // frees space allocated to node degree table

}

/** Sort the nodes by decreasing fractional solution.
 * This function compute an order. We get the datas computed by BiqCrunch.
 * We have for each nodes a fractional value so we will sort our nodes
 * by decreasing fractional value. The order computed will be store in the
 * global variable orderNode.
 * @param N the number of variable.
 * @param P0 the original problem.
 * @param node the current branch and bound tree node.
 */
void orderNodeByFracSol(Problem *P0, BobNode *node) {
    int N = P0->n - 1;
    int i = 0, n = 0, indMax;
    double max;
    double *tmpFracSol= (double*) malloc(N * sizeof(double));

    for (i = 0; i < N; i++)
        tmpFracSol[i] = node->fracsol[i];

    while (n != N) //While we don't fill the orderNode table
    {
        indMax = 0; //initialize the index of maximum value at 0
        max = tmpFracSol[0]; //initialize the maximum value at the the fractional value of node 0
        for (i = 0; i < N; i++) //for each nodes
            if (max < tmpFracSol[i]) //if the fractional value of node i is higher than the maximum value
            {
                indMax = i; //we set the index of maximum value at i
                max = tmpFracSol[i]; // we set the maximum value at the the fractional value of node i
            }
        //We put the fractional value of the node indmax at -2.0 because we want to be sure
        tmpFracSol[indMax] = -2.; //that we don't take the node indMax a second time in the order table
        orderNode[n] = indMax; //We put the node of maximum fractional value in the oerderNode table
        n++;

    }

    free(tmpFracSol);
}

/** Heuristic calls at the beginning of the program.
 * This heuristic is call once during the program.
 * It sort the nodes by increasing degrees and depending on this order
 * try to add each nodes one by one if it is possible (neoghboor not in
 * the independent set).
 * @param P0 the initial Problem.
 * @param x the ouput integer solution.
 * @param N the number of variables (nodes.)
 * @return the weight of the max independent set found.
 */
double first_stable_heurist(Problem *P0, int *x) {
    int N = P0->n - 1;
    int i = 0, k, j, current_node = 0, neighbourInStable;
    double heur = 0.;

    for (i = 0; i < N; i++)
        x[i] = 0;

    for (i = 0; i < N; i++) // for each node
    {
        k = 0;
        j = 0;
        neighbourInStable = 0; // boolean to check if a neighbour node is already in stable (1) or not (0)
        current_node = orderNode[i]; // get the node with the minimum degree
        while (k < M_adj->nnz && M_adj->i[k] != current_node) // looking for the edge of the current node in the adjacency sparse matrix
            k++;

        for (j = k; M_adj->i[j] == current_node; j++) // for each edge of current node
        {
            if (x[M_adj->j[j]]) // if the tail (neighbour node) is already in the stable
            {
                neighbourInStable = 1; // there is a neighbour in stable
                break; // we stop studiying this node
            }
        }
        if (neighbourInStable) // if we found a neighbour in the stable
            x[current_node] = 0; // we don't put the current node in the stable
        else
            // if no neighbour in the stable
            x[current_node] = 1; // we can put the current node in
    }

    for (i = 0; i < N; i++) // compute the heuristic value
        heur += weights[i] * (double) x[i];

    return heur;
}

/** Body of heuristic process.
 * This function will be called by two heuristics.
 * Before calling this function, we have to compute an order for the node
 * and then this function get the BiqCrunch solution and try to add nodes
 * following the previous computed order.
 * @param P0 the initial Problem.
 * @param x the integer output solution.
 * @param N the number of variables (nodes).
 * @return the weights of the max independent set found.
 */
double main_heuristic(Problem * P0, BobNode *node, int *x) {
    int N = P0->n - 1;
    double heur = 0.;
    int i, j, k;
    int current_node;
    int neighbourInStable;

    //G et the solution from BiqCrunch
    for (i = 0; i < N; i++) {
        if (node->xfixed[i])
            x[i] = node->sol.X[i];
        else
            x[i] = 0;
    }

    for (i = 0; i < N; i++) //For each nodes of the table orderNode
    {
        k = 0;
        j = 0;
        neighbourInStable = 0; //initialize boolean value at 0 to specify the current node has no neighboor in
        current_node = orderNode[i]; //the independent set and get the node i of the orderNode table
        if (!(node->xfixed[current_node])) //if the current node hasn't been fixed by BiqCrunch
        {
            while (k < M_adj->nnz && M_adj->i[k] != current_node) //we look for the neighboor list of the current node
                k++;

            for (j = k; M_adj->i[j] == current_node; j++) //for each neighboor of the current node
            {
                if (x[M_adj->j[j]]) //if the neighboor is in the independent set
                {
                    neighbourInStable = 1; //we specify that the current node has already a neighbootin the independent set
                    break; //we go out of the loop
                }
            }
            if (neighbourInStable) //if the current node has a neighboor in the independent set
                x[current_node] = 0; //we don't put it in the independent set
            else
                //else
                x[current_node] = 1; //we put it in
        }
    }

    for (i = 0; i < N; i++) //compute the heuristic value
        heur += weights[i] * (double) x[i];

    return heur;
}

/** Heuristic 1.
 * Description :
 *	1) Get the current solution from BiqCrunch.
 *	2) Sort the nodes by increasing degrees.
 *	3) Following the previous order computed, we try to
 *	add one by one each nodes.
 *	If he has no neighboor in the independent set we add
 *	it in the current independent set, else we don't add it.
 * @param P0 the initial Problem.
 * @param x the integer output solution.
 * @param N the number of variables (nodes).
 * @return the weights of the max independent set found.
 * @see main_heuristic()
 * @see orderNodeByDegree()
 */
double bab_stable_heuristic(Problem * P0, BobNode *node, int *x) {
    int N = P0->n - 1;
    orderNodeByDegree(M_adj, N);
    return main_heuristic(P0, node, x);
}

/** Heuristic 2.
 * Description :
 *	1) Get the current solution from BiqCrunch.
 *	2) Sort the nodes using the datas given by BiqCrunch.
 *	Each nodes has a fractional value. So we will sort our
 *	nodes by decreasing fractional value
 *	3) Following the previous order computed, we try to
 * 	add one by one each nodes.
 *	If he has no neighboor in the independent set we add
 *	it in the current independent set, else we don't add it.
 * @param P0 the original problem.
 * @param n the current branch and bound tree node.
 * @param x the integer output solution.
 * @param N the number of variables (nodes).
 * @return the weights of the max independent set found.
 * @see main_heuristic()
 * @see orderNodeByFracSol()
 */
double bab_stable_heuristic2(Problem * P0, BobNode *node, int *x) 
    //heuristic call during branch and bound
{
    orderNodeByFracSol(P0, node);
    return main_heuristic(P0, node, x);

}

/** Check if a value is in an input table.
 * @param val the value we are looking for.
 * @param tab the table where we looking for a value.
 * @param size the size of the input table.
 * @return 1 if we find it else 0.
 */
int in(int val, int *tab, int size) {
    int i;

    for (i = 0; i < size; i++)
        if (tab[i] == val)
            return 1;
    return 0;

}

/** Delete a value in an input table
 * Note : the size is an int* so the function modify
 * the size given in input after deleting the value.
 * @param val the value we want to delete.
 * @param tab the table where we want to delete the value.
 * @param size the size of the input table.
 * @return 0 if we don't found the value else 1.
 */
int deleteInTab(int val, int *tab, int *size) {
    int k = 0, i;

    while (k < *size && tab[k] != val)
        k++;

    if (k == *size)
        return 0;

    for (i = k; i < (*size) - 1; i++)
        tab[i] = tab[i + 1];

    (*size)--;

    return 1;
}

/** Sort a subset of node by decreasing degree.
 * This function permit to order a subset of node by decreasing degree.
 * The computed order is store in the table tab.
 * @param mAdj the sparse adjacency matrix of the input problem.
 * @param N the number of variables (nodes).
 * @param nodeSubset a subset of node we have to sort.
 * @param tab the ouput order table.
 */
void orderNodeByDecrDegree(Sparse* mAdj, int N, int *nodeSubset, int* tab) // Sort a subset of node by decreasing degree
{
    int* nodeDeg = (int*) malloc(N * sizeof(int)); // Table of node degree
    int i = 0, ind_n = 0, n = nodeSubset[0], nbnode = 0, degmax = 0;

    //Initialization
    for (i = 0; i < N; i++) {
        tab[i] = -1;
        nodeDeg[i] = 0;
    }

    degmax = 0; //initialize the maximum degree at 0

    for (i = 0; i < mAdj->nnz; i++) //we look for the node which are contained in the subset nodes table in the sparse adjacency matrix
    {

        if (in(mAdj->i[i], nodeSubset, N) == 0) //if the node isn't in the subset node
            i++; //we continue
        else //the current node is in the subset node
        {
            if (mAdj->i[i] == n) //if it is still the current node
                nbnode++; //we increase the degree of the current node
            else //else
            {
                nodeDeg[ind_n] = nbnode; //we add the degree of the node in the table of subset node degrees
                n = mAdj->i[i]; //we modify the current node
                ind_n++; //increase the index in the table of subset node degrees
                nbnode = 1; //reinitialize the degree of the new current node at one
            }
            if (nbnode > degmax) //if the current node degree is higher than themaximum degree
                degmax = nbnode; //we modify the maximum degree
        }
    }
    nodeDeg[N - 1] = nbnode; //we don't forget to had the degree of the latest node of subset node in the table of degrees

    for (i = N - 1; i >= 0; i--) //sort the nodes in subset node by decreasing degree
    {
        n = getIndMin(nodeDeg, N, degmax + 1); //get the index of minimum degree
        nodeDeg[n] = -1; //we specify that we have already add this node in the order node table
        tab[i] = nodeSubset[n]; //we pute the node at the last case unfilled
    }

    free(nodeDeg); // frees space allocated to node degree table
}

/** Heuristic 3.
 * Description :
 *	1) Get the current solution from BiqCrunch.
 *	2) Depending on a value called seuil, we will fixed
 *	our variables at 1 if its fractional value is higher
 *	and 0 if it is lower. (Note : we will probably have
 *	a non-admissible solution)
 *	3) We will correct the previous computed solution
 *		a) Build the adjacency matrix of the value fixed to
 *		one.
 *		b) Sort only the nodes we have fixed to one (we don't
 *		take the variable fixed by BiqCrunch) by decreasing degree.
 *		c) Following the previous order compute, we take the first
 *		node of the order (the highest degrees) we remove it from
 *		the current independent set. If it still exist edge in
 *		the current independent set problem back to step 3-a, else
 *		it's finished.
 * @param P0 the original problem.
 * @param n the current branch and bound tree node.
 * @param x the integer output solution.
 * @param N the number of variables (nodes).
 * @return the weights of the max independent set found.
 */
double bab_stable_heuristic3(Problem * P0, BobNode *node, int *x) {
    int N = P0->n - 1;
    double heur = 0., seuil = 0.2;
    int i, j, nz, nbNodeSetToOne = 0, nbNodeInStableTemp = 0;

    Sparse *SubM_adj;
    int *orderSubsetNode;
    int *indSetTemp = (int*) malloc(N * sizeof(int));
    int *nodeSetToOne = (int*) malloc(N * sizeof(int));

    for (i = 0; i < N; i++) {
        if (node->xfixed[i]) //if the variable is fixed we get the value from BiqCrunch datas
        {
            x[i] = node->sol.X[i];
            if (node->sol.X[i] == 1) //if the variable have been fixed to one
            {
                indSetTemp[nbNodeInStableTemp] = i; //we store the node set to one in indSetTemp
                nbNodeInStableTemp++;
            }
        } else {
            if (node->fracsol[i] >= seuil) //if the fractional solution of the current node is higher than seuil
            {
                x[i] = 1; //we fixed the variable at one
                nodeSetToOne[nbNodeSetToOne] = i; //we store this node as a node WE HAVE set to one
                indSetTemp[nbNodeInStableTemp] = i; //we store this node as a node set to one
                nbNodeInStableTemp++;
                nbNodeSetToOne++;
            } else
                //else
                x[i] = 0; //we fixed the variable at zero
        }
    }

    if (nbNodeInStableTemp && nbNodeSetToOne) {
        do {
            //------------------------ Build sparse adjacency matrix of the node set to one ------------------------------------------//

            SubM_adj = (Sparse *) malloc(sizeof(Sparse));
            SubM_adj->i = (int *) malloc(
                    nbNodeInStableTemp * nbNodeInStableTemp * sizeof(int));
            SubM_adj->j = (int *) malloc(
                    nbNodeInStableTemp * nbNodeInStableTemp * sizeof(int));
            SubM_adj->val = (double *) malloc(
                    nbNodeInStableTemp * nbNodeInStableTemp * sizeof(double));

            orderSubsetNode = (int*) malloc(nbNodeSetToOne * sizeof(int));

            nz = 0;
            for (i = 0; i < nbNodeInStableTemp; i++)
                for (j = 0; j < nbNodeInStableTemp; j++)
                    if (i != j)
                        if (getSparseAdjacencyMatrixValue(M_adj, indSetTemp[i],
                                    indSetTemp[j]) == 1.) {
                            SubM_adj->i[nz] = indSetTemp[i];
                            SubM_adj->j[nz] = indSetTemp[j];
                            SubM_adj->val[nz] = 1.;
                            nz++;
                        }

            SubM_adj->nnz = nz;

            //-------------------------------------------------------------------------------------------------------------------------//

            if (SubM_adj->nnz != 0) //If there is edges in the current independent set
            {
                orderNodeByDecrDegree(SubM_adj, nbNodeSetToOne, nodeSetToOne,
                        orderSubsetNode); //we sort the node set to one by decreasing degree

                //we take the node with the highest degree
                x[orderSubsetNode[0]] = 0; //we delete it from the current independent set
                deleteInTab(orderSubsetNode[0], indSetTemp,
                        &nbNodeInStableTemp); //we delete it from the subset node we have set to one
                deleteInTab(orderSubsetNode[0], nodeSetToOne, &nbNodeSetToOne); //we delete it from the table of nodes set to one

            }

            //free memory
            free(SubM_adj->i);
            free(SubM_adj->val);
            free(SubM_adj->j);
            free(SubM_adj);
            free(orderSubsetNode);

        } while (SubM_adj->nnz != 0 && nbNodeSetToOne > 0); //while we find an edge in the sparse adjacency matrix and the number of node we have set to one is positive
    }

    //free memory
    free(nodeSetToOne);
    free(indSetTemp);

    for (i = 0; i < N; i++) //compute the weight of the independent set found
        heur += weights[i] * (double) x[i];

    return heur;
}

//////////////////////////////////////////////////////////////////////////////////
