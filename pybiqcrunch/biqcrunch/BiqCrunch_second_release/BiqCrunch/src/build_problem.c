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
#include <stdlib.h>

#include "bb.h"
#include "biqcrunch.h"

extern FILE *output; // Output file
extern char heur_name[]; // defined in heur.c (specific to problem)
extern Problem *SP; // Original problem instance
extern Problem *PP; // Subproblem instance

//============== Read data from BC file ==================//
void Read_Data(char *StrgGrp) 
{
    FILE *f;
    int i, j, s;
    double Tmp;
    char instance_line[80];
    int count_sparse_eq = 0;
    int count_sparse_ineq = 0;

// Print on screen BiqCrunch banner
    fprintf(output,BANNER);

// Open input BiqCrunch file
    f = fopen(StrgGrp, "r");
    if (f == NULL) {
        fflush(stdout);
        fprintf(stderr, "Error: input file %s does not exist\n", StrgGrp);
        exit(1);
    }
    printf("Input file:  %s\n", StrgGrp);
    fprintf(output,"Input file: %s\n", StrgGrp);

    int sense, constraints, blocks, size, n_ineq;
    int BC_MaxProb;

    int line = 0;

    // read comments
    int res;
    do {
        res = 0;
        instance_line[0] = '\0';
        line++;
        res = fscanf(f, "%2[;*#\\^\n]%*[^\n]\n", instance_line);
    } while (res == 1 || 
            instance_line[0] == ';' || 
            instance_line[0] == '\\' || 
            instance_line[0] == '#' || 
            instance_line[0] == '*');

    // read min/max sense (sense=1 for MAX problem, sense=-1 for MIN problem)
    line++;
    INSTANCE_ERROR(f, fscanf(f, "%79[^\n]\n", instance_line) != 1, line);
    INSTANCE_ERROR(f, sscanf(instance_line, "%d", &sense) != 1, line);
    INSTANCE_ERROR(f, sense != 1 && sense != -1, line);
    BC_MaxProb = (sense == 1);

    // number of matrices
    INSTANCE_ERROR(f, fscanf(f,"%79[^\n]\n",instance_line) != 1, line);
    INSTANCE_ERROR(f, sscanf(instance_line, "%d", &constraints)!=1, line);
    INSTANCE_ERROR(f, constraints < 0, line);

    // number of blocks
    line++;
    INSTANCE_ERROR(f, fscanf(f,"%79[^\n]\n",instance_line) != 1, line);
    INSTANCE_ERROR(f, sscanf(instance_line, "%d", &blocks)!=1, line);
    INSTANCE_ERROR(f, blocks < 1 || blocks > 2, line);

    // size of the matrices
    line++;
    INSTANCE_ERROR(f, fscanf(f, "%79[^\n]%*[^\n]\n", instance_line)!=1, line);
    if (blocks == 1) {
        INSTANCE_ERROR(f,
                sscanf(instance_line, "{%d}", &size)!=1 && sscanf(instance_line, "%d", &size)!=1,
                line);
        n_ineq = 0;
    } else if (blocks == 2) {
        INSTANCE_ERROR(f,
                sscanf(instance_line, "{%d, -%d}", &size, &n_ineq)!=2 && sscanf(instance_line, "%d, -%d", &size, &n_ineq)!=2,
                line);
    } else {
        INSTANCE_ERROR(f, 1, line);
    }
    INSTANCE_ERROR(f, size <= 0 || n_ineq < 0 || n_ineq > constraints, line);

    alloc(SP, Problem);
    allocProblem(SP, size, n_ineq, constraints - n_ineq + size, BC_MaxProb);
    BobPbSize = SP->n - 1;

    double *right_hand;
    alloc_vector(right_hand, SP->mB + SP->mA, double);

    // right-hand side values
    line++;
    for (s = 0; s < constraints; s++) {
        INSTANCE_ERROR(f, fscanf(f, "%lf", &(right_hand[s]))!=1, line);
    }

    int block; // block of the matrix read
    int last = -1; // last matrix index
    int ineq = 0; // if the constraint is an inequality
    int ge = 0; // if the inequality is a >= inequality
    int ii, jj;
    double * matrix;
    double * matrix0;
    int isInteger;
    double coeff;

    alloc_matrix(matrix, SP->n, double);
    alloc_matrix(matrix0, SP->n, double);

    // Set matrix = matrix0 = 0.0
    for (jj = 0; jj < SP->n; jj++) {			
        for (ii = 0; ii < SP->n; ii++) {
            matrix[ii + jj * SP->n] = 0.0;
            matrix0[ii + jj * SP->n] = 0.0;
        }
    }

    s = 0;
    while (!feof(f)) {

        last = s;

        // cell of the matrix
        line++;
        INSTANCE_ERROR(f,
                fscanf(f, "%d %d %d %d %lf \n", &s, &block, &i, &j, &Tmp) != 5,
                line);
        INSTANCE_ERROR(f,
                (n_ineq == 0 && block == 2) || (s < 0 || s > constraints),
                line);
        INSTANCE_ERROR(f,
                block == 1 && ((i < 1 || i > size) || (j < 1 || j > size)),
                line);
        INSTANCE_ERROR(f,
                block == 2 && ((i < 1 || i > n_ineq) || (j < 1 || j > n_ineq)),
                line);

        // Check that objective coefficients are integer
        if (s == 0) {
            if (i==j) coeff = Tmp; else coeff = 2*Tmp;
            isInteger = (coeff == (int) coeff);
            if (!isInteger) {
                printf("Error reading line %d of file %s.\n", line, StrgGrp);
                printf("The coefficient of x[%d]*x[%d] is %f in the objective function.\n", 
                        i, j, coeff);
                printf("A fractional coefficient  is not allowed in the objective function since\n");
                printf("the objective function must always be integer for any feasible solution.\n");
                printf("To fix this problem, multiply all objective coefficients by a sufficiently\n");
                printf("large integer, and try again.\n");
                exit(1);
            }
        }

        if (block == 2) {
            ineq = 1;
            if (Tmp < 0.0) { // I've found >=, change the sign of the values
                ge = 1;
            }

        } else { // block == 1
            if (last == s) {
                // for the OBJECTIVE function if the problem is not maximization
                if (s == 0 && !SP->max_problem) Tmp = -Tmp;

                if (i==j) {
                    matrix0[(SP->n) * (i - 1) + (i - 1)] += Tmp;
                } else {
                    matrix0[(SP->n) * (j - 1) + (i - 1)] += Tmp;
                    matrix0[(SP->n) * (i - 1) + (j - 1)] += Tmp;
                }

            } else { // last != s
                if (ge) {
                    for (jj = 0; jj < SP->n; jj++)
                        for (ii = 0; ii < SP->n; ii++)
                            matrix0[ii + jj * SP->n] = -matrix0[ii + jj * SP->n];

                    right_hand[last - 1] = -right_hand[last - 1];
                }

                //////////////////////////////////////////////////////////////////// F.R. dec 2013
                // Transform 0-1 problem into -1 1 problem 
                // Contribution of Tmp * xixj => 0.25 * Tmp * ( y0^2 + yiy0 + yjy0 + yiyj ) 
                for (jj = 0; jj < SP->n; jj++) {
                    for (ii = 0; ii < SP->n; ii++) {

                        matrix[(SP->n - 1) + (SP->n - 1) * SP->n] += 0.25*matrix0[ii + jj * SP->n];  // y0^2 = 1

                        matrix[(SP->n - 1) + ii * SP->n] += 0.125*matrix0[ii + jj * SP->n]; // yiy0
                        matrix[ii + (SP->n - 1) * SP->n] += 0.125*matrix0[ii + jj * SP->n];

                        matrix[(SP->n - 1) + jj * SP->n] += 0.125*matrix0[ii + jj * SP->n]; // yjy0
                        matrix[jj + (SP->n - 1) * SP->n] += 0.125*matrix0[ii + jj * SP->n];

                        // F.R. 2014 : handle diagonal terms correctly
                        if (ii!=jj) {
                            matrix[(SP->n) * jj + ii] += 0.125*matrix0[ii + jj * SP->n];  // yiyj
                            matrix[(SP->n) * ii + jj] += 0.125*matrix0[ii + jj * SP->n];
                        } else {
                            matrix[(SP->n - 1) + (SP->n - 1) * SP->n] += 0.25*matrix0[ii + jj * SP->n];  // y0^2 = 1
                        }

                    }
                }				
                ////////////////////////////////////////////////////////////////////

                if (last == 0) {
                    addProblemObjective(SP, matrix);
                } else {
                    if (ineq) {
                        INSTANCE_ERROR(f,
                                addProblemInequality(SP, matrix, right_hand[last - 1], &count_sparse_ineq) == 0,
                                line);
                    } else {
                        INSTANCE_ERROR(f,
                                addProblemEquality(SP, matrix, right_hand[last - 1], &count_sparse_eq)== 0,
                                line);
                    }
                }
                ineq = 0;
                ge = 0;

                // Set matrix = matrix0 = 0.0
                for (jj = 0; jj < SP->n; jj++) {
                    for (ii = 0; ii < SP->n; ii++) {
                        matrix[ii + jj * SP->n] = 0.0;
                        matrix0[ii + jj * SP->n] = 0.0;
                    }
                }

                last = s;

                if (i==j) {
                    matrix0[(SP->n) * (i - 1) + (i - 1)] += Tmp;
                } else {
                    matrix0[(SP->n) * (j - 1) + (i - 1)] += Tmp;
                    matrix0[(SP->n) * (i - 1) + (j - 1)] += Tmp;
                }

            } // END last != s

        } // END block == 1

    } //EOF

    if (last != -1) {
        if (ge) {
            for (jj = 0; jj < SP->n; jj++)
                for (ii = 0; ii < SP->n; ii++)
                    matrix0[ii + jj * SP->n] = -matrix0[ii + jj * SP->n];

            right_hand[last - 1] = -right_hand[last - 1];
        }

        //////////////////////////////////////////////////////////////////// F.R. dec 2013
        // Transform 0-1 problem into -1 1 problem 
        // Contribution of Tmp * xixj => 0.25 * Tmp * ( y0^2 + yiy0 + yjy0 + yiyj ) 
        for (jj = 0; jj < SP->n; jj++) {
            for (ii = 0; ii < SP->n; ii++) {

                matrix[(SP->n - 1) + (SP->n - 1) * SP->n] += 0.25*matrix0[ii + jj * SP->n];  // y0^2 = 1

                matrix[(SP->n - 1) + ii * SP->n] += 0.125*matrix0[ii + jj * SP->n]; // yiy0
                matrix[ii + (SP->n - 1) * SP->n] += 0.125*matrix0[ii + jj * SP->n];

                matrix[(SP->n - 1) + jj * SP->n] += 0.125*matrix0[ii + jj * SP->n]; // yjy0
                matrix[jj + (SP->n - 1) * SP->n] += 0.125*matrix0[ii + jj * SP->n];

                // F.R. 2014 : handle diagonal terms correctly
                if (ii!=jj) {
                    matrix[(SP->n) * jj + ii] += 0.125*matrix0[ii + jj * SP->n];  // yiyj
                    matrix[(SP->n) * ii + jj] += 0.125*matrix0[ii + jj * SP->n];

                } else {
                    matrix[(SP->n - 1) + (SP->n - 1) * SP->n] += 0.25*matrix0[ii + jj * SP->n];  // y0^2 = 1

                }

            }
        }								
        ////////////////////////////////////////////////////////////////////

        if (last == 0) {
            addProblemObjective(SP, matrix);
        } else {
            if (ineq) {
                INSTANCE_ERROR(f,
                        addProblemInequality(SP, matrix, right_hand[last - 1], &count_sparse_ineq) == 0,
                        line);
            } else {
                INSTANCE_ERROR(f,
                        addProblemEquality(SP, matrix, right_hand[last - 1], &count_sparse_eq) == 0,
                        line);
            }
        }
    }

    INSTANCE_ERROR(f, last != constraints, line);

    fclose(f);

    addBooleanConstraint(SP, matrix, &count_sparse_eq);

    alloc(PP, Problem);
    allocCopyProblem(PP, SP);

    free(matrix);
    free(matrix0);
    free(right_hand);

    // OUTPUT information on instance
    fprintf(output, "Solving as a %s problem\n",
            (SP->max_problem) ? "MAXIMIZATION" : "MINIMIZATION");
    fprintf(output, "Problem Size = %d Number of equalities = %d Number of inequalities = %d\n",
            BobPbSize, SP->mB, SP->mA);
    fprintf(output, "Using %s heuristic\n", heur_name);
}



/*
 * Add the boolean constraints to the model.
 * @param prob: the problem to add the constraints
 * @param matrix: matrix used to create the constraint
 */
void addBooleanConstraint(Problem * prob, double * matrix, int * count_sparse_eq) 
{ 
    // F.R. dec 2013 : now not using U and U^T
    int i = 0, u, t;

    for (i = 0; i < prob->n ; i++) {

        for (t = 0; t < prob->n; t++)
            for (u = 0; u < prob->n; u++)
                matrix[u + t * prob->n] = 0.0;

        matrix[i + i * prob->n] = 1.0;
        addProblemEquality(prob, matrix, 1.0, count_sparse_eq);
    }
}


/* 
 * addProblemObjective 
 */
void addProblemObjective(Problem * p, double * matrix) 
{
    int i, j;
    for (j = 0; j < p->n; j++) {
        for (i = 0; i < p->n; i++) {
            p->Q[i + j * p->n] = matrix[i + j * p->n];
        }
    }

    sparseMatrix(&(p->Qs), matrix, p->n);
}


/* 
 * addProblemInequality 
 */
int addProblemInequality(Problem * p, double * matrix, int a, int * count_sparse_ineq) 
{
    if (*count_sparse_ineq >= p->mA)
        return 0;

    sparseMatrixA(*count_sparse_ineq,p,&(p->As[*count_sparse_ineq]), matrix, p->n);

    p->a[*count_sparse_ineq] = a;
    (*count_sparse_ineq)++;

    return 1;
}


/* 
 * addProblemEquality 
 */
int addProblemEquality(Problem * p, double * matrix, int b, int * count_sparse_eq) 
{
    if (*count_sparse_eq >= p->mB)
        return 0;

    sparseMatrixB(*count_sparse_eq,p,&(p->Bs[*count_sparse_eq]), matrix, p->n);

    p->b[*count_sparse_eq] = b;
    (*count_sparse_eq)++;

    return 1;
}

