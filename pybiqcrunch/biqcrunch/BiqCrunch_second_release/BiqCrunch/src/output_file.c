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
#include <sys/stat.h>


extern FILE *output;


/*
 * Open the output file.
 * If the output file exist it add a number a counter at the end of the filename
 * @param instance: name of the instance to use as base for the output file name
 * @return 1 if succeed in the creation of the output file, 0 otherwise
 */
int createOutputFile(char *instance) 
{
    char output_path[200];

    sprintf(output_path, "%s.output", instance);
    struct stat buffer;
    int counter = 1;
    // Check if the file already exists, if so aappend _<NUMBER> to the end of the 
    // output file name
    while (stat(output_path, &buffer) == 0)
        sprintf(output_path, "%s.output_%d", instance, counter++);
    output = fopen(output_path, "w");

    // Display output file
    printf("Output file: %s\n", output_path);

    if (!output)
        return 0;
    return 1;
}


/*
 * Close the output file
 */
void closeOutputFile() 
{
    fclose(output);
}

