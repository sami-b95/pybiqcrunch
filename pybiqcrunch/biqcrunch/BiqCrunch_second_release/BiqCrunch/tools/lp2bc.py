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
lp2bc reads a binary quadratic optimization problem in LP format
and generates the corresponding BiqCrunch input file.

Usage: ./lp2bc.py file.lp > file.bc

lp2bc follows the CPLEX LP format specifications as described on the webpage:

http://pic.dhe.ibm.com/infocenter/cosinfoc/v12r4/topic/ilog.odms.cplex.help/CPLEX/File_formats_reference/topics/LP.html

Note however that lp2bc only recognizes the following keywords:

    min, minimize, max, maximize, subject to, st, s.t., bin, binary, 
    binaries, end

In particular, lp2bc does not recognize the following keywords:

    bound, bounds, general generals, gen, semi-continuous, semi, semis, sos

In addition, lp2bc can handle square brackets around quadratic terms, but it is
not required to have them.


Example lp file:
================

\ This is a comment. Any text following a backslash is ignored.

Maximize

  obj:  20 z1*z3 + 26 z1*z4 + 
        23 z2*z3 +  8 z2*z5 + 
        32 z3*z4 + 13 z4*z5 

Subject to

  c1: z1 + z2 + z3 + z4 + z5 = 3

  c2: 12 z1*z3 + 24 z1*z4 + 
      14 z2*z3 + 16 z2*z5 + 
      28 z3*z4 + 12 z4*z5 <= 30

Binary

  z1 z2 z3 
  z4 z5

End
'''


###############################################################################
from __future__ import print_function
import sys


###############################################################################
__date__ = "April 27, 2016"

__authors__ = "Nathan Krislock, Jérôme Malick, and Frédéric Roupin"

__copyright__ = '''
Copyright 2010-2016, Nathan Krislock, Jérôme Malick, and Frédéric Roupin'''

__license__ = '''
This file is part of BiqCrunch.

BiqCrunch is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BiqCrunch is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BiqCrunch.  If not, see <http://www.gnu.org/licenses/>.'''

__credits__ = '''
The finite state machines used for parsing text in this code were developed 
using the helpful advice found in the following article:

http://robertgawron.blogspot.fr/2011/04/finite-state-machine-for-parsing.html

Also, thanks to Olivier Huber at INRIA Grenoble for suggesting the use of 
finite state machines for parsing LP files.'''


###############################################################################
# Define constants
inside_brackets = 'inside brackets'
outside_brackets = 'outside brackets'


###############################################################################
def parse_equation(line, bracket_state=outside_brackets):
    '''
    Parse an equation specified in the string using a finite state machine.

    Returns a dictionary containing the coefficients, variables, comparison 
    type, right-hand-side, and the equation label.

    No variable, constant, or sense indicator may be split by a return.

    A variable cannot begin with a number, a period, or the letter e/E.

    Each constraint must begin on a new line.
    '''

    # Separate the equation label and the equation
    label_line = line.split(':')
    if len(label_line) == 1:
        label = None
    elif len(label_line) == 2:
        label, line = label_line
    else:
        print("Error: cannot have more than one ':' in an equation", 
                file=sys.stderr)
        return None

    # Remove all spaces and final return
    line = line.replace(' ', '')
    line = line.strip()

    # States of the finite state machine
    initial = 'initial'
    coefficient = 'coefficient'
    variable = 'variable'
    operator = 'operator'
    left_bracket = 'left_bracket'
    right_bracket = 'right_bracket'
    comparison = 'comparison'
    rhs = 'rhs'
    final = 'final'

    # List of valid characters that may appear in variable names
    valid_chars = ['!', '"', '#', '$', '%', '&', '(', ')', ',', '.', 
                   ';', '?', '@', '_', "'", '`', '{', '}', '~']

    # Initialize return objects
    coeffs, variables, in_bracket = [], [], []
    comptype, rhsval = None, None
    bracket_factor = 1.0

    # Initialize finite state machine
    state, buff = initial, ''

    for char in line + '\0':

        is_variable_char = char.isalpha() | (char in valid_chars)
        is_num_dot = char.isdigit() | (char == '.') 

        # Initial state
        if state == initial:
            # We can start with a coefficient, an operator, a variable,
            # a comparison, a left bracket, a right bracket, or end-of-line

            if is_num_dot:
                state = coefficient
                buff = char

            elif char in ['-', '+']:
                state = operator
                buff = char

            elif is_variable_char & (char not in ['e', 'E']):
                state = variable
                coeffs.append(1.0)  # Coefficient of variable is 1
                buff = char

            elif char in ['=', '>', '<']:
                state = comparison
                buff = char

            elif char == '[':
                state = left_bracket
                buff = char

                if bracket_state == outside_brackets:
                    bracket_state = inside_brackets

                else:
                    print("Error: cannot open a bracket without closing the " 
                        + "previous bracket.", file=sys.stderr)
                    return None

            elif char == ']':
                state = right_bracket
                buff = char

                if bracket_state == inside_brackets:
                    bracket_state = outside_brackets

                else:
                    print("Error: cannot close a bracket without opening a " 
                        + "bracket.", file=sys.stderr)
                    return None

            elif char == '\0':
                state = final

            else:
                print("error at character (%c)" % char, file=sys.stderr)
                return None

        # Reading coefficient term
        elif state == coefficient:
            # After a coefficient we must have a variable

            if is_num_dot | (char in ['e', 'E', '+', '-']):
                buff += char

            elif is_variable_char & (char not in ['e', 'E']):
                state = variable
                try:
                    coeffs.append(float(buff))

                except ValueError:
                    print("Error: cannot convert '%s' to a number." % buff, 
                            file=sys.stderr)
                    return None

                buff = char

            elif char in ['\\', '\0']:
                print("Error: cannot end a line with a coefficient.",
                        file=sys.stderr)
                return None

            else:
                print("error at character (%c)" % char, file=sys.stderr)
                return None

        # Reading variable name
        elif state == variable:   
            # After a variable we must have an operator, a comparison, a
            # right bracket, or end-of-line

            if is_variable_char | char.isdigit() | (char in ['^', '*']):
                buff += char

            elif char in ['-', '+']:
                state = operator
                variables.append(buff)
                in_bracket.append(bracket_state)
                buff = char

            elif char in ['=', '<', '>']:
                state = comparison
                variables.append(buff)
                in_bracket.append(bracket_state)
                buff = char

            elif char == ']':
                state = right_bracket
                variables.append(buff)
                in_bracket.append(bracket_state)
                buff = char

                if bracket_state == inside_brackets:
                    bracket_state = outside_brackets

                else:
                    print("Error: cannot close a bracket without opening a " 
                        + "bracket.", file=sys.stderr)
                    return None

            elif char == '\0':
                state = final
                variables.append(buff)
                in_bracket.append(bracket_state)

            else:
                print("error at character (%c)" % char, file=sys.stderr)
                return None

        # Reading operator
        elif state == operator:
            # After an operator we either have a coefficient, a variable,
            # a left bracket, or end-of-line (if operator is '+')

            if is_num_dot:
                state = coefficient
                buff += char  # Keep the operator in the buffer

            elif is_variable_char & (char not in ['e', 'E']):
                state = variable
                coeffs.append(float(buff + '1'))  # Coeff of variable is 1
                buff = char

            elif char == '[':
                state = left_bracket

                if bracket_state == outside_brackets:
                    bracket_state = inside_brackets

                else:
                    print("Error: cannot open a bracket without closing the " 
                        + "previous bracket.", file=sys.stderr)
                    return None

                if buff == '-':
                    bracket_factor *= -1.0
                    
                buff = char

            elif char == '\0':
                if buff == '+':
                    state = final

                else:
                    print("Error: cannot end a line with '-'.",
                        file=sys.stderr)
                    return None

            else:
                print("error at character (%c)" % char, file=sys.stderr)
                return None

        # Reading left bracket
        elif state == left_bracket:
            # After a left bracket we must have a coefficient, an operator, 
            # a variable, or end-of-line

            if is_num_dot:
                state = coefficient
                buff = char

            elif char in ['+', '-']:
                state = operator
                buff = char

            elif is_variable_char & (char not in ['e', 'E']):
                state = variable
                coeffs.append(1.0)  # Coefficient of variable is 1
                buff = char

            elif char == '\0':
                state = final

            else:
                print("error at character (%c)" % char, file=sys.stderr)
                return None

        # Reading right bracket
        elif state == right_bracket:
            # After a right bracket we can have '/2', an operator, or a
            # comparison

            if buff + char == ']/':
                buff += char

            elif buff + char == ']/2':
                buff += char
                bracket_factor *= 0.5

            elif char in ['-', '+']:
                state = operator
                buff = char
            
            elif char in ['=', '<', '>']:
                state = comparison
                buff = char
                
            elif char == '\0':
                state = final

            else:
                print("error at character (%c)" % char, file=sys.stderr)
                return None

        # Reading comparison
        elif state == comparison:
            # After a comparison we must have a RHS

            if char == '=':
                # Ignore any '=' that follows '<', '>', '='
                state = comparison

            elif char in ['>', '<']:
                print("Error: only inequalities '>=', '<=' are allowed.", 
                        file=sys.stderr)
                return None

            elif is_num_dot | (char in ['+', '-']):
                state = rhs
                comptype = buff
                buff = char

            elif char == '\0':
                print("Error: cannot end a line with a comparison.",
                        file=sys.stderr)
                return None

            else:
                print("error at character (%c)" % char, file=sys.stderr)
                return None

        # Reading RHS
        elif state == rhs:
            # After a RHS we must end

            if is_num_dot | (char in ['e', 'E', '+', '-']):
                state = rhs
                buff += char

            elif char == '\0':
                state = final                
                try:
                    rhsval = float(buff)

                except ValueError:
                    print("Error converting right-hand-side '%s' to a number." 
                            % buff, file=sys.stderr)
                    return None

            else:
                print("error at character (%c)" % char, file=sys.stderr)
                return None

        # Finished reading
        elif state == final:
            pass

        # This should never happen
        else:
            print("Error: something impossible happened.", file=sys.stderr)
            return None

    if len(coeffs) == len(variables) == len(in_bracket):
        terms = zip(coeffs, variables, in_bracket)

    else:
        print("Error: number of coefficients and variables does not match.", 
                file=sys.stderr)
        return None

    return dict(label=label, terms=terms, comptype=comptype, rhsval=rhsval,
            bracket_state=bracket_state, bracket_factor=bracket_factor)


###############################################################################
def lp_parser(filename):
    '''
    This is a finite state machine to parse an lp file.

    Returns the problem specified as a dictionary.
    '''

    try:
        f = open(filename)

    except IOError:
        print("Error: file '%s' does not exist." % filename,
                file=sys.stderr)
        return None

    # States of the finite state machine
    initial = 'initial'
    objective = 'objective'
    constraints = 'constraints'
    binaries = 'binaries'
    final = 'final'

    bad_keywords = ['bound', 'bounds', 'general', 'generals', 'gen', 
            'semi-continuous', 'semi', 'semis', 'sos']

    # Initialize the finite state machine
    state = initial

    obj = dict(sense='', terms=[], label='', bracket_factor=1.0)
    constrs, binvars = [], []

    bracket_state = outside_brackets
    found_left_bracket = found_right_bracket = False

    line_number = 0

    for line in f:
        line_number += 1

        # Remove any comments
        ind = line.find('\\')
        if ind > -1:
            line = line[:ind]

        # Remove initial and trailing whitespace and newlines
        line = line.strip()

        # Detect brackets
        if '[' in line:
            if found_left_bracket:
                print("Error: cannot have two '[' in the same equation",
                        file=sys.stderr)
                print("error reading line %d:\n" % line_number, line, 
                        file=sys.stderr)
                return None

            else:
                found_left_bracket = True

        if ']' in line:
            if found_right_bracket:
                print("Error: cannot have two ']' in the same equation",
                        file=sys.stderr)
                print("error reading line %d:\n" % line_number, line, 
                        file=sys.stderr)
                return None
            
            else:
                found_right_bracket = True


        if state == initial:
            # We can start with 'minimize' or 'maximize'

            if line == '':
                pass

            elif line.lower() in ['minimize', 'min']:
                state = objective
                obj['sense'] = 'minimize'

            elif line.lower() in ['maximize', 'max']:
                state = objective
                obj['sense'] = 'maximize'

            else:
                print("Error: must begin problem with 'minimize' or "
                    + "'maximize'", file=sys.stderr)
                print("error reading line %d:\n" % line_number, line, 
                        file=sys.stderr)
                return None

        elif state == objective:
            # After objective we must have constraints or binaries

            if line == '':
                pass

            elif line.lower() in ['subject to', 'st', 's.t.']:
                state = constraints
                # Reset bracket detectors
                found_left_bracket = found_right_bracket = False

            elif line.lower() in ['binaries', 'binary', 'bin']:
                state = binaries

            elif line.lower() in bad_keywords:
                print("Unsupported keyword used on line %d:\n" 
                        % line_number, line, file=sys.stderr)
                return None

            else:
                res = parse_equation(line, bracket_state=bracket_state)

                if res != None:
                    bracket_state = res['bracket_state']

                    obj['terms'] += res['terms']
                    obj['bracket_factor'] *= res['bracket_factor']

                    if res['label'] != None:
                        obj['label'] = res['label']

                else:
                    print("error reading line %d:\n" % line_number, line, 
                            file=sys.stderr)
                    return None

        elif state == constraints:
            # After constraints we must have binaries

            if line == '':
                pass

            elif line.lower() in ['binaries', 'binary', 'bin']:
                state = binaries

            elif line.lower() in bad_keywords:
                print("Unsupported keyword used on line %d:\n" 
                        % line_number, line, file=sys.stderr)
                return None

            else:
                res = parse_equation(line, bracket_state=bracket_state)

                if res != None:
                    bracket_state = res['bracket_state']

                    if res['rhsval'] != None:  # End of constraint
                        # Reset bracket detectors
                        found_left_bracket = found_right_bracket = False

                    if len(constrs) == 0:
                        constrs.append(res)

                    elif constrs[-1]['rhsval'] != None:
                        # If the last constraint is complete, start a new
                        # constraint
                        constrs.append(res)

                    else:
                        # The last constraint is not complete, so add the
                        # terms, comparison type, and rhs value to it
                        constrs[-1]['terms'] += res['terms']
                        constrs[-1]['comptype'] = res['comptype']
                        constrs[-1]['rhsval'] = res['rhsval']
                        constrs[-1]['bracket_factor'] *= res['bracket_factor']

                        if res['label'] != None:
                            constrs[-1]['label'] = res['label']

                else:
                    print("error reading line %d:\n" % line_number, line, 
                            file=sys.stderr)
                    return None

        elif state == binaries:
            # After binaries we can have 'end' or end-of-file

            if line == '':
                pass

            elif line.lower() in ['end']:
                state = final

            elif line.lower() in bad_keywords:
                print("Unsupported keyword used on line %d:\n" 
                        % line_number, line, file=sys.stderr)
                return None

            else:
                binvars += line.split(' ')

        elif state == final:
            pass

    f.close()

    problem = dict(obj=obj, constrs=constrs, binvars=binvars)

    return problem


###############################################################################
def find_ijval(term, binvars, bracket_factor):
    '''
    Returns (ii, jj, val) or None if there is a problem.

    If the term is linear, jj = n + 1.
    '''

    n = len(binvars)

    (coeff, variable, bracket_state) = term

    if '*' in variable:
        product_list = variable.split('*')

        if len(product_list) == 2:
            var1, var2 = product_list

            try:
                ii = binvars.index(var1) + 1

            except ValueError:
                print("All variables must be binary: %s is not binary." 
                        % var1, file=sys.stderr)
                return None

            try:
                jj = binvars.index(var2) + 1

            except ValueError:
                print("All variables must be binary: %s is not binary." 
                        % var2, file=sys.stderr)
                return None

            # If term is inside the brackets, multiply the coeff by the
            # bracket_factor
            if bracket_state == inside_brackets:
                coeff *= bracket_factor

        else:
            # Error
            return None

    elif '^' in variable:
        square_term_list = variable.split('^')

        if len(square_term_list) == 2:
            var, exponent = square_term_list

            if exponent == '2':
                try:
                    ii = binvars.index(var) + 1

                except ValueError:
                    print("All variables must be binary: %s is not binary." 
                            % var, file=sys.stderr)
                    return None

                jj = ii

                # If term is inside the brackets, multiply the coeff by the
                # bracket_factor
                if bracket_state == inside_brackets:
                    coeff *= bracket_factor

            else:
                # Error
                return None

        else:
            # Error
            return None

    else:
        # Linear term
        try:
            ii = binvars.index(variable) + 1

        except ValueError:
            print("Error: %s is not binary." % variable, file=sys.stderr)
            return None

        # Put linear term in the last column
        jj = n + 1

        # If term is inside the brackets, report an error
        if bracket_state == inside_brackets:
            print("Error: linear term %s is inside brackets." %
                    variable, file=sys.stderr)
            return None

    # If term is not on the diagonal, divide coefficient by two
    if ii == jj:
        val = coeff

    else:
        val = coeff/2.                    

    return (ii, jj, val)


###############################################################################
def BC_file(prob):
    '''
    Takes a problem specified in the output of lp_parser and 
    prints the BC file.
    '''

    # Gather problem information
    obj = prob['obj']
    constrs = prob['constrs']
    binvars = prob['binvars']
    sense = obj['sense']
    comptypes = [constr['comptype'] for constr in constrs]
    rhsvals = [constr['rhsval'] for constr in constrs]
    m = len(constrs)
    n = len(binvars)

    # Output the variable names
    k = 0
    output = '# List of binary variables:'
    for var in binvars:
        k += 1
        output += '\n# %4d:  %s' % (k, var)

    # Output the optimization sense
    if sense == 'maximize':
        output += '\n1 = max problem'

    elif sense == 'minimize':
        output += '\n-1 = min problem'

    else:
        print("Error: sense must be 'minimize' or 'maximize'.",
                file=sys.stderr)
        return None

    # Output the number of constraints
    output += '\n%d = number of constraints' % m

    # Output the number of blocks and the size of the blocks
    num_ineqs = comptypes.count('<') + comptypes.count('>')
    if num_ineqs > 0:
        output += '\n2 = number of blocks'
        output += '\n%d, %d' % (n + 1, -num_ineqs)

    else:
        output += '\n1 = number of blocks'
        output += '\n%d' % (n + 1)

    # Output the RHS values
    if m > 0:
        output += '\n'
    for rhsval in rhsvals:
        output += '%.1f ' % rhsval

    # Output the objective function
    bracket_factor = obj['bracket_factor']
    for term in obj['terms']:
        res = find_ijval(term, binvars, bracket_factor)

        if res != None:
            (ii, jj, val) = res
            output += '\n0 1 %d %d %.1f' % (ii, jj, val)

        else:
            # Error
            return None

    # Output the constraints
    k = 0
    ineq_counter = 0
    for constr in constrs:
        k += 1  # Increment the constraint number

        bracket_factor = constr['bracket_factor']
        for term in constr['terms']:
            res = find_ijval(term, binvars, bracket_factor)

            if res != None:
                (ii, jj, val) = res
                output += '\n%d 1 %d %d %.1f' % (k, ii, jj, val)

            else:
                # Error
                return None

        if constr['comptype'] == '<':
            ineq_counter += 1
            output += '\n%d 2 %d %d 1.0' % (k, ineq_counter, ineq_counter)

        elif constr['comptype'] == '>':
            ineq_counter += 1
            output += '\n%d 2 %d %d -1.0' % (k, ineq_counter, ineq_counter)

        else:
            pass

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
            problem = lp_parser(FILENAME) 

            if problem != None:
                RES = BC_file(problem)

                if RES != None:
                    print(RES)

                else:
                    print("lp2bc: an error occured.", file=sys.stderr)

            else:
                print("lp2bc: an error occured.", file=sys.stderr)

    else:
        print('Usage:', sys.argv[0], 'file.lp > file.bc', file=sys.stderr)
        print('For detailed help:', sys.argv[0], 'help', file=sys.stderr)
