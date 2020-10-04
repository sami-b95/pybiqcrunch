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

/********************************************************/
/************ List of all global variables **************/
/********************************************************/

/********************************************************/
BiqCrunchParameters params; //only modified by read_parameters.c. Used in read_parameters.c, bob_functions.c, bounding_procedure.c, evaluate_node.c, update_solution.c, inequalities.c
Problem *SP; // Original problem instance. Only modified by in bob_functions.c (when reading inputfile). Used in bounding_procedure.c, build_problem.c, evaluate_node.c, update_solution.c
FILE *output; // output file

/********************************************************/
int stopped = 0; // True if the algorithm stopped at root node or after a time limit
double best_eval; // Best evaluation value. This value keeps the value of the more promising node if the algorithm is interrupted at root node or after a time limit
double TEMPI; // initial time (modified once)

/********************************************************/
/*************** Specific to node ***********************/
/********************************************************/
// These working variables are (re-)used by each node. Allocation (before B&B) and free (after B&B) are done once in alloc_free.c.

/* Subproblem variables (some variables are fixed  */
Problem *PP; // Subproblem instance. Used in Bob_functions.c, bounding_procedure.c, build_problem.c, bounding_procedure.c, evaluate_node.c
double *SubMat; // Construction matrices (size BobSizePb). Used in evaluate_node.c
double *SubMatLin; // Construction vector (size BobSizePb). Used in evaluate_node.c
Sparse SubMatS; //Construction matrices, sparse representation. Used in evaluate_node.c
double *tempMatrix; // Construction vector (size BobSizePb + 1). Used in evaluate_node.c
double *X; // Stores current X (primal solution). Violated inequalities are computed from X.

/* Triangle Inequalities variables */
Inequality *Cuts; // vector (MaxNineqAdded) of triangle inequality constraints
Inequality *List; // list of inequalities

/* L-BFGS-B variables */
double gradInorm; // norm of the gradient of f for equality constraints
double gradEnorm; // norm of the gradient of f for equality constraints
double *g, f; // gradient and function values
double *y; // dual variables
double *binf; // lower bounds on the variables y
double *bsup; // upper bounds on the variables y
int *nbd; // indicates which variables are bounded
double *wa;  // Double workspace for L-BFGS-B
int *iwa; // Integer workspace for L-BFGS-B

/* dsyevr function variables (projection) */
int M = 0; // number of eigenvalues
double *W; // contains the eigenvalues 
double *Z; // contains the eigenvectors
int *ISUPPZ; // dim = 2*max(1,M), support of the eigvecs. in Z
double *WORK; // dim = LWORK
int *IWORK; // dim = LIWORK
int sizeWORK, sizeIWORK; // sizes of WORK and IWORK
