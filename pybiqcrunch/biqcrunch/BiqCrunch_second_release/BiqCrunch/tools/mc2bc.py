#! /usr/bin/env python
# -*- coding: utf-8 -*-
# *****************************************************************************
# *                                                                           *
# *  BiqCrunch is a semidefinite-based solver for binary quadratic problems.  *
# *  It uses a branch-and-bound method featuring an improved semidefinite     *
# *  bounding procedure, mixed with a polyhedral approach. BiqCrunch uses     *
# *  particular input files format to describe the combinatorial problems.    *
# *                                                                           *
# *   Copyright (C) 2010-2016 Nathan Krislock, Jérôme Malick, Frédéric Roupin *
# *                                                                           *
# *                 http://www-lipn.univ-paris13.fr/BiqCrunch/                *
# *									      *
# *****************************************************************************
#									      *
#    This program is free software: you can redistribute it and/or modify     *
#    it under the terms of the GNU General Public License as published by     *
#    the Free Software Foundation, either version 3 of the License, or        *
#    (at your option) any later version.                                      *
#                                                                             *
#    This program is distributed in the hope that it will be useful,          *
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           *
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
#    GNU General Public License for more details.                             *
#                                                                             *
#    You should have received a copy of the GNU General Public License        *
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    *
#                                                                             *
# *****************************************************************************
'''
mc2bc reads a graph file in the following format:

    n m
    i j d
    i j d
    ...

where n is the number of nodes, m is the number of edges, and 'i j d' is a edge
from node i to node j with weight d.

mc2bc outputs the corresponding BiqCrunch file for the Max-Cut problem of the
input graph.

Usage: ./mc2bc.py graphfile > file.bc
'''

###############################################################################
from __future__ import print_function
import numpy as np
import sys

###############################################################################
def adjacency_matrix(graphfile):
    """
    A = adjacency_matrix(graphfile)

    Forms the adjacency matrix A from the graph specified in graphfile.

    The file must follow the format:

        n m
        i j d
        i j d
        ...

    Parameters
    ----------
    graphfile : string
        Contains the name of the file that specifies the graph.

    Returns
    -------
    A : (n, n) ndarray
        A is the adjacency matrix of the graph specified in graphfile.
    """

    try:
        f = open(graphfile)

    except IOError:
        sys.exit("Error: file '%s' does not exist." % graphfile)

    for line in f:
        data = line.split()

        if len(data) == 2:
            n = int(data[0])  # m = int(data[1])
            A = np.zeros((n, n), dtype=float)

        elif len(data) == 3:
            ii = int(data[0])
            jj = int(data[1])
            dd = float(data[2])
            A[ii-1, jj-1] += dd
            A[jj-1, ii-1] += dd

        else:
            print('Error:  data = ', data)

    f.close()

    return A


###############################################################################
def lagrangian(A):
    """
    L = lagragian(A)

    Computes the Lagrangian maxtrix L of the adjacency matrix A.

    Parameters
    ----------
    A : (n, n) ndarray
        A is the adjacency matrix of a graph.

    Returns
    -------
    L : (n, n) ndarray
        L is the Lagrangian matrix of the graph specified by the adjacency
        matrix A.
    """

    L = np.diag([sum(d) for d in A]) - A
    return L


###############################################################################
def BC_file(filename):

    L = lagrangian(adjacency_matrix(filename))
    n = len(L) - 1  # We set the last variable to zero, so there is one fewer
                    # variables than the size of L

    output = '1 = max problem'
    output += '\n0 = number of constraints'
    output += '\n1 = number of blocks'
    output += '\n%d' % (n + 1)

    for ii in range(n):  # ii = 0,...,n-1
        for jj in range(ii, n):  # jj = ii,...,n-1
            val = L[ii,jj]
            if val != 0.0:
                output += '\n0 1 %d %d %f' % (ii+1, jj+1, L[ii,jj])

    return output


###############################################################################
###############################################################################

# Call the main function
if __name__ == "__main__":

    if len(sys.argv) == 2:

        if sys.argv[1] == 'help':
            print(__doc__, file=sys.stderr)

        else:

            FILENAME = sys.argv[1]
            res = BC_file(FILENAME)

            if res != None:
                print(res)

            else:
                print("mc2bc: an error occured.", file=sys.stderr)

    else:
        print('Usage:', sys.argv[0], 'graphfile > file.bc', file=sys.stderr)
        print('For detailed help:', sys.argv[0], 'help', file=sys.stderr)
