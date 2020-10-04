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

#include <stdlib.h>

/*
 * Maximum number of variables
 */
#define NMAX 1024

/*
 * Solution of the problem.
 * This structure defines the content of a solution of the problem.
 */
typedef struct BobSolution {
    /*
     * Vector X.
     * Binary vector that stores the solution of the branch-and-bound algorithm
     */
    int X[NMAX];
} BobSolution;

/*
 * Node of the branch-and-bound tree.
 * Structure that represent a node of the branch-and-bound tree and stores all the 
 * useful information.
 */
typedef struct BobNode {
    /*
     * Node information.
     *
     */
    BobTNdInfo BobNdInfo; // (int Size, int Off)
    /*
     * Number of fixed variables.
     * Integer variable that stores the number of fixed variables in the current node.
     */
    int level;
    int xfixed[NMAX];
    BobSolution sol;
    double fracsol[NMAX];
    BobTPri Pri; // (int Eval, int Depth)
} BobNode;
