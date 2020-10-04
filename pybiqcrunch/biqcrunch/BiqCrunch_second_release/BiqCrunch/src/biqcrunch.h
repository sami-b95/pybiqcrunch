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

#include <stdio.h>
#include <math.h>

/*
 *  Maximum number of inequalities allowed to add
 */
#define MaxNineqAdded 10000

/*
 *  Constant numbers
 */

#define BIG_NUMBER 1e+9
#define mmax 10 // mmax is the maximum number of limited memory corrections (BFGS).
#define scaleIneq 1.0 / (1.0 + sqrt(3.0)) // scales the inequality constraints (and dual vars)
#define scaleEq 1.0 // scales the equality constraints (and dual vars)

/*
 * Different heuristics are used at different stages in the 
 * branch and bound method.
 */
#define PRIMAL_HEUR 1
#define SDP_BOUND_HEUR 2
#define ROUNDING_HEUR 3

#define TRANSP      1
#define NOTRANSP    0

/*
 * Branching strategies
 */
#define LEAST_FRACTIONAL  0
#define MOST_FRACTIONAL   1
#define CLOSEST_TO_ONE    2

/*
 * We track where we find a new solution
 */
#define NEWSOL_LEAF 1 // when we reach a leaf
#define NEWSOL_HEUR1 2 // the first heuristic (primal heuristic)
#define NEWSOL_HEUR2 3 // the second heuristic (rounding during the bound computation)
#define NEWSOL_HEUR3 4 // the third heuristic (rounding after the node evaluation)

// macro to handle the errors in the input reading
#define INSTANCE_ERROR(file,cond,line)\
    if((cond)){\
        fprintf(stderr, \
                "\nError: the instance doesn't respect the standard on line %d\n", \
                line);\
        fclose(file);\
        exit(1);\
    }

// macro used to print the banner in the console and in the output file
#define BANNER \
    "*******************************************************************************\n"\
    "*                            BIQ CRUNCH 2.0 Solver                            *\n"\
    "*******************************************************************************\n"\
    "| Copyright(C) 2010-2016 N. Krislock, J. Malick, F. Roupin                    |\n"\
    "| BIQ CRUNCH uses L-BFGS-B by C. Zhu, R. Byrd, J. Nocedal and BOB 1.0 by PNN  |\n"\
    "| Team of PRiSM Laboratory.                                                   |\n"\
    "|                                                                             |\n"\
    "| L-BFGS-B is distributed under the terms of the New BSD License. See the     |\n"\
    "| website http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html for more    |\n"\
    "| information. BOB is free software. For more information visit the website   |\n"\
    "| http://www.prism.uvsq.fr/~blec/index.php?p=./rech/so&lang=fr                |\n"\
    "*******************************************************************************\n"

// macros for allocating vectors and matrices
#define alloc_vector(var, size, type)\
    var = (type *) calloc((size) , sizeof(type));\
    if(var == NULL){\
        fprintf(stderr,\
                "\nError: Not enough memory allocating "#var" variable in %s line %d\n",\
                __FILE__,__LINE__);\
        exit(1);\
    }
#define alloc(var, type) alloc_vector(var, 1, type)
#define alloc_matrix(var, size, type) alloc_vector(var, size*size, type)

// Defining parameters structure and default values
#ifndef PARAM_FIELDS
#define PARAM_FIELDS \
    P(double,   alpha0,              "%lf",             1e-1) \
    P(double,   scaleAlpha,          "%lf",              0.5) \
    P(double,   minAlpha,            "%lf",             5e-5) \
    P(double,   tol0,                "%lf",             1e-1) \
    P(double,   scaleTol,            "%lf",             0.95) \
    P(double,   minTol,              "%lf",             1e-2) \
    P(int,      withCuts,            "%d",                 1) \
    P(double,   gapCuts,             "%lf",            -5e-2) \
    P(int,      cuts,                "%d",               500) \
    P(int,      minCuts,             "%d",                50) \
    P(int,      nitermax,            "%d",              2000) \
    P(int,      minNiter,            "%d",                12) \
    P(int,      maxNiter,            "%d",               100) \
    P(int,      maxNAiter,           "%d",                50) \
    P(int,      scaling,             "%d",                 1) \
    P(int,      root,                "%d",                 0) \
    P(int,      heur_1,              "%d",                 1) \
    P(int,      heur_2,              "%d",                 1) \
    P(int,      heur_3,              "%d",                 1) \
    P(int,      soln_value_provided, "%d",                 0) \
    P(int,      soln_value,          "%d",                 0) \
    P(int,      time_limit,          "%d",                 0) \
    P(int,      branchingStrategy,   "%d",   MOST_FRACTIONAL) \
    P(unsigned, seed,                "%u",              2016) \
    P(int,      NBGW1,               "%d",                 1) \
    P(int,      NBGW2,               "%d",                10) \
    P(int,      local_search,        "%d",                 1)
#endif

typedef struct BiqCrunchParameters { 
#define P(type, name, format, def_value) type name;
    PARAM_FIELDS
#undef P
} BiqCrunchParameters;

/*
 * Sparse matrices are stored using the following Sparse structure.
 */
typedef struct Sparse {
    int *i;
    int *j;
    double *val;
    int nnz;
} Sparse;

// Define structure for storing all inequalities
typedef struct Inequality {
    int i;
    int j;
    int k;
    int type;
    double value;
    double y;
} Inequality;

/*
 * The main problem and any subproblems are stored using the following
 * Problem structure.
 */
typedef struct Problem {
    /*
     * Objective matrix in DENSE format
     */
    double *Q;
    /*
     * Objective matrix in SPARSE format
     */
    Sparse Qs;
    int n; // size of Q
    Sparse *As; // list of sparse matrices for the inequality constraints
    double *a; // right-hand-side vector of inequality constraints
    int mA; // number of inequality constraints
    Sparse *Bs; // list of sparse matrices for the equality constraints
    double *b; // right-hand-side vector of equality constraints
    int mB; // number of equality constraints
    int max_problem; // 1 if it is a max problem, and 0 if it is a min problem
    int NIneq;
} Problem;


// Declaration of Bob functions to
int Bob_ULBGet(void); // get the bound
int Bob_ULBUpd(int, BobSolution *); // update the bound
void Bob_ULBInit(int, BobSolution *); // initialize the bound
int Bob_ULBSup(int); // compare the bound
int Bob_ExpCtrl(int, BobTExpCt *);
int Bob_STNEXP(); // increment and return the number of explored nodes
int Bob_NEXP(); // return the number of evaluated nodes
int Bob_STEVL(); // increment and return the number of evaluated nodes
int Bob_EVL(); // return the number of evaluated nodes

/* bob_functions.c */
void Init_UQP();
double temps_CPU();
int getBranchingVariable(BobNode *Nodp);
void printFinalOutput(FILE *file, int evalcount);
void printSolution(FILE *file);
int isLeafNode(BobNode *node);
BobNode *newNode(BobNode *parentNode);
void processCommandLineArguments(int argc, char **argv);
void initializeBobSolution();

/* build_problem */
void Read_Data(char *StrgGrp);
void addBooleanConstraint(Problem *p, double *matrix, int * count_sparse_eq);
void addProblemObjective(Problem *p, double *matrix);
int addProblemInequality(Problem *p, double *matrix, int a, int * count_sparse_ineq);
int addProblemEquality(Problem *p, double *matrix, int b, int * count_sparse_eq);

/* read_parameters.c */
void readParameter(char *path);

/* output_file.c */
int createOutputFile(char *instance);
void closeOutputFile();

/* evaluate_node.c */
int countFixedVariables(BobNode *node);
double Evaluate(BobNode *node, Problem *SP, Problem *PP);
double getSubMatrix(double *subM, Sparse *spMatrix, int *fixed, int N);
double getLinear(double *subMatrixLin, Sparse *spMatrix, int *sol, int *fixed, int N);
double getSubMatrixSparse(Sparse *subM, Sparse *spMatrix, int *fixed, int N);
double getConstant(Sparse *spMatrix, int *sol, int *fixed, int N);
double getFixedValue(BobNode *node, Problem *SP);
void createSubproblem(BobNode *node, Problem *SP, Problem *PP);
void buildObjective(BobNode *node, Problem *SP, Problem *PP);
void buildConstraints(BobNode *node, int num, double *RHSsource, Sparse *Matsource, double *RHSdest, Sparse *Matdest);

/* update_solution.c */
int updateSolution(BobNode *node, int *x, int credit);
void Bob_PrintSolution();
double BC_evaluateSolution(int *sol);
int BC_isFeasibleSolution(int *sol);

/* sparse_matrix.c */
void sparseMatrixA(int s, Problem *prob, Sparse *spMatrix, double *denseMatrix, int N);
void sparseMatrixB(int s, Problem *prob, Sparse *spMatrix, double *denseMatrix, int N);
void sparseMatrix(Sparse *spMatrix, double *denseMatrix, int N);
void sparseMatrixScaled(Sparse *spMatrix, double *denseMatrix, int N, double scale);

/* alloc_free.c */
void allocMemory();
void allocProblem(Problem *p, int N, int mA, int mB, int max_problem);
void allocCopyProblem(Problem *, Problem *);
void allocSDPbound(Problem *P0);
void allocProj(Problem *P0);
void freeMemory();
void freeProblem(Problem *p);
void freeSDPbound();
void freeProj();

/* heur.c */ // defined in each BiqcCunch/problems/ subdirectories. Used to call BiqCrunch heuristics (src/rounding.c) and to add new ones (user defined).
int BC_allocHeuristic(Problem *SP);
void BC_freeHeuristic();
double BC_runHeuristic(Problem *SP, Problem *P, BobNode *node, int *x, int heur_code);
void BC_FixVariables(BobNode *node, int ic, int xic); // Call by function Bob_GenChild (in bob_function.c)
// Empty body in problems/generic/heur.c. Example in problem/max-indep-set/heur.c
// May be written by user to take advantage of particular constraints to fix several variables when node->sol.X[ic] = xic (0 or 1)
// ic and xic are such that:  node->xfixed[ic] = 1; node->sol.X[ic] = xic;

/* rounding.c */
int roundBound(double bound, int max_problem);
double primalHeuristic(Problem *P0, int *x);
double sdpBoundHeuristic(Problem *P0, BobNode *node, int * x);
double roundingHeuristic(Problem *P0, BobNode *node, int * x);
double GWheuristic(Problem *P0, Problem *P, BobNode *n, int *x, int num);
int update_best(int *x, int *tempx, double *best, Problem *P0);
int local_search(BobNode *node, int *x, double *val, Problem *P0);

/* bounding_procedure.c */
int pruneTest(double bound, int max_problem);
double SDPbound(BobNode *node, Problem *SP, Problem *PP);
int calllbfgsb(double *y, Problem *P, double alpha, double tol, double beta, double fixedvalue, int max_problem,int *nbit);
// Declare the L-BFGS-B function
void setulb_( 
        int *n, // number of variables
        int *m, // size of the limited memory
        double *x, // current iterate (length n)
        double *l, // lower bounds (length n)
        double *u, // upper bounds (length n)
        int *nbd, // indicates which vars are bounded
        // nbd[i] = 0 : x[i] unbounded
        // nbd[i] = 1 : x[i] has only a lower bound
        // nbd[i] = 2 : x[i] has both lower and upper bounds
        // nbd[i] = 3 : x[i] has only an upper bound
        double *f, // value of function at the point x
        double *g, // value of the gradient at the point x (length n)
        double *factr, // termination tolerance
        // factr = 1.d+12 for low accuracy;
        //         1.d+7  for moderate accuracy;
        //         1.d+1  for extremely high accuracy.
        double *pgtol, // projected gradient tolerance
        // (suppress this termination test by pgtol = 0)
        double *wa, // workspace (length (2mmax + 4)nmax + 12mmax^2 + 12mmax)
        int *iwa, // workspace (length 3nmax)
        char *task, // character array (length 60)
        // 'START'    : when starting
        // 'FG'        : user must evaluate function f and gradient g
        // 'NEW_X'    : user can decide whether to continue or stop
        // 'CONV'    : termination test has been satisfied
        // 'ABNO'    : abnormal termination
        // 'ERROR'    : error in input parameters
        // 'STOP'    : set by user to stop L-BFGS-B
        int *iprint, // set level of output
        // iprint<0    no output is generated;
        // iprint=0    print only one line at the last iteration;
        // 0<iprint<99 print also f and |proj g| every iprint iterations;
        // iprint=99   print details of every iteration except n-vectors;
        // iprint=100  print also the changes of active set and final x;
        // iprint>100  print details of every iteration including x and g;
        char *csave, // character array (length 60)
        int *lsave, // logicial array (length 4)
        int *isave, // integer array (length 44)
        double *dsave // double array (length 29)
        );

/* inequalities.c */
double getViolatedCuts(double *X, int N, Inequality *List, int *ListSize);
double updateInequalities(Problem *PP, double *y, int *NumAdded, int *NumSubtracted, int count);
double eval_ineq(double *XX, int N, int type, int ii, int jj, int kk);

/* bound.c */
void sim(Problem *P, double *y, double *f, double *g, double alpha);
void A(int mode, Problem *P, double *y, double *f, double *g, double alpha);
void B(int mode, Problem *P, double *y, double *f, double *g, double alpha);
void projSDP(double *f, Problem *P);

extern double dnrm2_(
        int *n, double *x, int *incx);
extern void dscal_(
        int *n, double *da, double *x, int *incx);
extern void dcopy_(
        int *n, double *dx, int *incx, double *dy, int *incy);
extern void dsyevr_(
        char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA, 
        double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, 
        double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK, int *LWORK, 
        int *IWORK, int *LIWORK, int *INFO);
extern void dsyrk_(
        char *UPLO, char *TRANS, int *N, int *K, double *ALPHA, double *A, 
        int *LDA, double *BETA, double *C, int *LDC);
extern double dlansy_(
        char *NORM, char *UPLO, int *N, double *A, int *LDA, double *WORK);

void print_vector(double *vec, int N);
void print_symmetric_matrix(double *Mat, int N);

/* redirect_output.f */
void redirect_output_(); // Redirect the Fortran output to /dev/null
