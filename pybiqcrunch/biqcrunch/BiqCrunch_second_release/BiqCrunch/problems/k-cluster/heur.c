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
// Authors : N. Krislock, F. Roupin 2012
// Updated version of heuristics for the k-cluster problem
// in weighted graphs F. Roupin 2013

#include <stdio.h>

#include "bb.h"
#include "biqcrunch.h"

extern BiqCrunchParameters params;

/*
 * Name of the heuristic printed in the output file
 */
char heur_name[] = "k-cluster (weighted)";


void heurist_kc(int * nsol, int * sol, int N, double *Matrix);
double primal_heuristic_kc(Problem *P0, int *X);
double sdpBoundHeuristic_kc(Problem *P0, BobNode *node, int *X);
double rounding_heuristic_kc(Problem *P0, BobNode *node, int *X);

int K;

int BC_allocHeuristic(Problem *P) {
    K = (P->b[0]) / 2; // K is the first value of the right-hand-side vector
    return 1;
}

void BC_freeHeuristic() {
    // nothing to free
}

// User can fix other variables as soon as x[ic]=xic (i.e. 0 or 1), if contraints in problem allows it.
// Call by function Bob_GenChild (in bob_function.c)
// Recall that ic and xic are such that:  node->xfixed[ic] = 1; node->sol.X[ic] = xic;
void BC_FixVariables(BobNode *node, int ic, int xic) {
}

double BC_runHeuristic(Problem *P0, Problem *P, BobNode *node, int * x, int heuristic_code)
{
    int N = P0->n - 1; // problem size
    int temp_x[N]; // temporary vector to store heuristic solutions
    double heur_val = 0; // value ofr a given heuristic
    double best; // best value 


    best = heur_val;
    switch (heuristic_code) {
        case PRIMAL_HEUR:
            printf("Primal heuristic\n");
            heur_val = primal_heuristic_kc(P0, temp_x);
            update_best(x, temp_x, &best, P0);
            break;
        case SDP_BOUND_HEUR:
            // Specific greedy heuristic for k-cluster
            heur_val = sdpBoundHeuristic_kc(P0, node, temp_x);
            update_best(x, temp_x, &best, P0);
            //		printf("H%d = %f ",heuristic_code,heur_val);
            // GW heuristic
            heur_val = GWheuristic(P0, P, node, temp_x, params.NBGW1*P->n);
            update_best(x, temp_x, &best, P0);
            //		printf(" %f\n",best);
            break;
        case ROUNDING_HEUR:
            // Specific greedy heuristic for k-cluster
            heur_val = rounding_heuristic_kc(P0, node, temp_x);
            update_best(x, temp_x, &best, P0);
            // GW heuristic
            heur_val = GWheuristic(P0, P, node, temp_x, params.NBGW2*P->n);
            update_best(x, temp_x, &best, P0);
            break;
        default:
            printf("Choosen heuristic doesn't exist\n");
            exit(1);
    }

    // If x is feasible, perform a local search around x for a better solution
    if (params.local_search && BC_isFeasibleSolution(x)) {
        local_search(node, x, &best, P0);
    }

    // return the value of the heuristic
    return best;

}

double sdpBoundHeuristic_kc(Problem *P0, BobNode *node, int *x) {
    int N = P0->n - 1;
    int i, k;
    double heur_val;

    double MaxFracSol;
    for (i = 0; i < N; i++) {
        // set entries that have not been fixed to 0
        if (node->xfixed[i]) {
            x[i] = node->sol.X[i];
        } else {
            x[i] = 0;
        }
    }

    // compute Kprime = the number of fixed variables with mysol[i]==1
    int Kprime = 0;
    for (i = 0; i < N; i++) {
        if (node->xfixed[i]) {
            Kprime += x[i];
        }
    }

    int index;
    for (k = 0; k < K - Kprime; k++) {
        MaxFracSol = -1.0; // node->fracsol is 0..1
        index = -1;
        for (i = 0; i < N; i++) {
            if (node->fracsol[i] > MaxFracSol && x[i] == 0 && !node->xfixed[i]) {
                MaxFracSol = node->fracsol[i];
                index = i;
            }
        }
        if (index == -1) {
        } else {
            x[index] = 1;
        }
    }

    // evaluate heur_val
    heur_val = BC_evaluateSolution(x);

    return heur_val;

}

double rounding_heuristic_kc(P0, node, X)
    Problem *P0; // Instance of the problem
    BobNode * node; // Node of the tree
    int *X; // the returned 0-1 vector
{
    //compute Kprime = the number of fixed variables with X[i]==1
    int Kprime = 0;
    int k, i;
    int N = P0->n - 1;

    for (i = 0; i < N; i++) {
        if (node->xfixed[i]) {
            Kprime += node->sol.X[i];
        }
    }

    //initialize X
    for (i = 0; i < BobPbSize; i++) {
        // set entries that have not been fixed to 0
        if (!(node->xfixed[i])) {
            X[i] = 0;
        } else
            X[i] = node->sol.X[i];
    }

    // find the indices corresponding to the largest K-Kprime entries of
    // node->fracsol and set those entries of node->sol.X to 1
    //    printf("Updating node->sol.X ...\n");
    double MaxFracSol;
    for (k = 0; k < K - Kprime; k++) {
        MaxFracSol = -1.0; // node->fracsol is 0..1
        int index = -1;
        // loop through indices looking for the largest node->fracsol[i] that is
        // not fixed and has not been set to 1
        for (i = 0; i < N; i++) {
            if (node->fracsol[i] > MaxFracSol && X[i] == 0 && !(node->xfixed[i])) {
                MaxFracSol = node->fracsol[i];
                index = i;
            }
        }
        if (index == -1) {
        } else {
            // index found, set X[index] = 1
            X[index] = 1;
        }
    }

    // evaluate the weight of the K-cluster
    double heur_val = BC_evaluateSolution(X);
    return heur_val;
}

//========= run_heuristic ==================//
double primal_heuristic_kc(Problem *P0, int *X) {
    /*
     * Run the heuristic (heurist_kc in heur.c)
     */
    int N = P0->n - 1;
    int * nsol;
    int * sol;

    double *Matrix;
    int i, j;
    double heur_val;

    // Compute the Matrix from P0->Q, omitting row n and column n
    Matrix = (double *) malloc(N * N * sizeof(double));

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            Matrix[i + j * N] = 0.0;

    for (i = 0; i < P0->Qs.nnz; i++) {
        int ii = P0->Qs.i[i];
        int jj = P0->Qs.j[i];
        double val = P0->Qs.val[i];

        if (ii < N && jj < N) {
            if (ii != jj) {
                Matrix[ii + jj * N] = 8.0 * val;
                Matrix[jj + ii * N] = Matrix[ii + jj * N];
            } else
                Matrix[jj + ii * N] = 0.0;
        }
    }

    // sol : vertices in the k-cluster
    // nsol : vertices not in the k-cluster
    nsol = (int *) malloc(N * sizeof(int));
    sol = (int *) malloc(N * sizeof(int));
    for (i = 0; i < N; i++) {
        nsol[i] = 0;
        sol[i] = 0;
    }

    // void heurist_kc(int * nsol, int * sol, int nbfin, int N, double *Matrix)
    heurist_kc(nsol, sol, N, Matrix);

    // copy sol (not 0/1 vector) into X (0/1 vector)
    for (i = 0; i < N; i++)
        X[i] = 0;
    for (i = 0; i < N; i++) {
        if (sol[i])
            X[sol[i] - 1] = 1;
    }

    // compute the value of heur_val
    heur_val = BC_evaluateSolution(X);

    free(Matrix);
    free(nsol);
    free(sol);

    return heur_val;
    // return 0;

}

//Structures used by the heuristic
typedef struct chaine chaine;
typedef struct sommet sommet;

struct chaine {
    int hd;
    chaine * tl;
};
struct sommet {
    double degre;
    int deg;
    chaine * voisins;
};

int nbinit; //total number of vertices
sommet ** tab; //array representing the graph

//Adds an element to a string
// c : string
// m : element to add
chaine * ajmaillon(chaine * c, int m) {
    if (c == NULL) {
        c = malloc(sizeof(chaine));
        c->hd = m;
        c->tl = NULL;
        return c;
    } else {
        chaine * res = malloc(sizeof(chaine));
        res->hd = m;
        res->tl = c;
        return res;
    }
}

//Removes a given element from a string, if existing
// c : string
// m : pattern to remove
chaine* rmmaillon(chaine * c, int m) {
    if (c != NULL) {
        if (c->hd == m) {
            chaine * tmp = c->tl;
            free(c);
            return tmp;
        } else {
            c->tl = rmmaillon(c->tl, m);
            return c;
        }
    } else
        return (NULL);
}

//Frees a string from memory
void free_chaine(chaine * c) {
    if (c != NULL) {
        free_chaine(c->tl);
        free(c);
    }
}

//Gets the data from the graph file and put it into the array 'tab'
//Mind that a vertex in the graph corresponds to the cell in the array
//  with the same number MINUS 1
void get_data(int N, double *Matrix) {
    int som1, som2, i, j, aretes;

    // nbinit = total number of vertices
    // aretes = total number of edges
    //fscanf(f, "%d %d", &nbinit, &aretes);
    nbinit = N;
    aretes = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < i; j++) {
            if (Matrix[i + j * N]) {
                aretes++;
            }
        }
    }

    tab = malloc(nbinit * sizeof(sommet*));
    for (i = 0; i < nbinit; i++) {
        tab[i] = malloc(sizeof(sommet));
        (tab[i])->degre = 0;
        (tab[i])->deg = 0;
        (tab[i])->voisins = NULL;
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < i; j++) {
            if (Matrix[i + j * N]) {

                som1 = i + 1;
                som2 = j + 1;

                (tab[som1 - 1])->degre += Matrix[i + j * N];
                ((tab[som1 - 1])->deg)++;
                (tab[som1 - 1])->voisins = ajmaillon((tab[som1 - 1])->voisins,
                        (som2 - 1));

                (tab[som2 - 1])->degre += Matrix[i + j * N];
                ((tab[som2 - 1])->deg)++;
                (tab[som2 - 1])->voisins = ajmaillon((tab[som2 - 1])->voisins,
                        (som1 - 1));
            }
        }
    }

}

//Sends back the number of cell representing the min degree in the array 'tab'
int degre_min() {
    double min = 32767;
    int res = 0;
    int i;
    for (i = 0; i < nbinit; i++) {
        if ((tab[i]) != NULL && (tab[i])->degre < min
                && (tab[i])->deg > (-1))
        {
            res = i;
            min = ((tab[i])->degre);
        }
    }

    //	printf("Degre min %d\n",res+1);

    return res;

}

//Removes a vertex from the array 'tab'
void retire_som(int ind_som, int N, double *Matrix) {
    if ((tab[ind_som])->voisins != NULL) {
        while ((tab[ind_som])->deg > 0) {
            (tab[((tab[ind_som])->voisins)->hd])->voisins = rmmaillon(
                    (tab[((tab[ind_som])->voisins)->hd])->voisins, ind_som);

            (tab[((tab[ind_som])->voisins)->hd])->degre -= Matrix[ind_som + (((tab[ind_som])->voisins)->hd) * N];
            (tab[ind_som])->degre -=Matrix[ind_som + (((tab[ind_som])->voisins)->hd) * N];
            ((tab[((tab[ind_som])->voisins)->hd])->deg)--;
            ((tab[ind_som])->deg)--;

            (tab[ind_som])->voisins = rmmaillon((tab[ind_som])->voisins,
                    ((tab[ind_som])->voisins)->hd);
        }
        (tab[ind_som])->degre = (-32767);
        (tab[ind_som])->deg = (-1);
    }
}

//Gets the final selection from 'tab' and sets 'sol' with it
void retour(int * nres, int * res) {
    int i, ind = 0, ind2 = 0;
    for (i = 0; i < nbinit; i++) {
        if ((tab[i])->deg != (-1)) {
            res[ind] = (i + 1);
            ind++;
        } else {
            nres[ind2] = (i + 1);
            ind2++;
        }
        free_chaine((tab[i])->voisins);
        free(tab[i]);
    }
    free(tab);
}

void TwoOpt(int * nres, int * res, int N, double *Matrix) {
    // local search : swap two vertices until no progress is made
    // res : vertices in the k-cluster
    // nres : vertices not in the k-cluster

    int i, j, k, imax, jmax, dummy;
    int progress = 1;
    double delta, deltamax;

    imax = -1;
    jmax = -1; // added to remove compiler warning about uninitialized variables
    // -- Nathan Krislock (21/04/2012)
    while (progress) {
        deltamax = 0;
        for (i = 0; i < K; i++) {
            for (j = 0; j < N - K - 1; j++) {
                // take one vertex in the k-cluster and another one not in the kcluster
                delta = 0.0;
                // compute the difference between the two solutions
                for (k = 0; k < K; k++)
                    if (k != i) {
                        delta += Matrix[(res[k] - 1) * N + nres[j] - 1]
                            - Matrix[(res[k] - 1) * N + res[i] - 1];
                    }
                if (delta > deltamax) {
                    deltamax = delta;
                    jmax = j;
                    imax = i;
                }
            }
        }
        // progress ?
        if (deltamax > 0) {
            dummy = nres[jmax];
            nres[jmax] = res[imax];
            res[imax] = dummy;
        } else
            progress = 0;
    }
}

// main function, computes the heuristic
// IN:
//	sol  : final solution
//	nbfin: K of the kcluster problem
//	graph: graph file
// EFFECT:
//	'sol' contains the selected subgraph
void heurist_kc(int * nsol, int * sol, int N, double *Matrix) {
    int i;

    //initializes nbinit and fills tab
    get_data(N, Matrix);

    /*	for (i=0;i<N;i++)
        {
        printf("%f\n",tab[i]->degre);
        }
        */
    //removes the vertices of least degrees, one by one
    for (i = 0; i < (nbinit - K); i++) {
        retire_som(degre_min(), N, Matrix);
    }

    //fills sol with the selected subgraph
    // sol : vertices in the k-cluster
    // nsol : vertices not in the k-cluster

    retour(nsol, sol);

    //	printf("Vertices in k-cluster :\n");
    //	for (i=0;i<N;i++) printf(" %d",sol[i]);
    //	printf("\n");
    //	printf("Vertices not in k-cluster :\n");
    //	for (i=0;i<N;i++) printf(" %d",nsol[i]);
    //	printf("\n");

    TwoOpt(nsol, sol, N, Matrix);

}
