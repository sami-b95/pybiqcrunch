#!/usr/bin/env python
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

import sys

# Control the arguments
if len( sys.argv ) != 2:
	sys.exit("biq2bc: missing argument. Usage: bq2bc.py input > output")	

##########################
#    INSTANCE READING    #
##########################

# Open the file
try:
	input_file  = open(sys.argv[1], 'r')
except IOError:
	sys.exit("biq2bc: file "+sys.argv[1]+" does not exist")
	
# Read all the instance
instance = input_file.readlines()
input_file.close() # Close the input file

##########################
#  DENSE MATRIX CREATION #
##########################

# Control first line: N and M
line = instance[0].split()
if len(line) != 2:
	sys.exit("biq2bc: the first line must contain the number of nodes N and the number of edges M")

# Size of the matrices
dim = int(line[0]) + 1

# Alloc the objective matrix
matrix = [[0]*dim for i in range(1,dim+1)]

# Compute the objective matrix
for cell in instance[1:]:
	# Check if the line is correct
	row = cell.split();
	if len(row) != 3:
		sys.exit("biq2bc: an entry of the graph must contain the nodes I and J and the value V")
		
	cell_split = [float(s) if '.' in s else int(s) for s in row] # Get an entry from the sparse graph
	
	i = cell_split[0] # Get first node
	j = cell_split[1] # Get second node
	t = float(cell_split[2]) # Get value
	
	# Check values obtained
	if i < 1 or i > dim or j < 1 or j > dim:
		sys.exit("biq2bc: the values exceed matrix size")

	if i == j: # Put the values on the last column/last row if they are on the diagonal
		matrix[i-1][dim-1] = t /2.0
		matrix[dim-1][i-1] = t /2.0
	else: # Put the values inside the matrix
		matrix[i-1][j-1] = t
		matrix[j-1][i-1] = t
		
##########################
#   OUTPUT BC INSTANCE   #
##########################

# Print the instance
print "1 =max problem"  # indicate that this is a maximization problem
print "0"  # 0 contraints for the max-cut problem
print "1"  # The matrices will have just one block
print str(dim)  # Size of the matrix
for i,r in enumerate(matrix):  # Print in sparse format the objective matrix
	for j, v in enumerate(r):
		if j >= i and v != 0:
			print "0 1 %d %d %f" % (i+1, j+1, v) 
