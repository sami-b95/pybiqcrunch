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
#include <string.h>

#include "bb.h"
#include "biqcrunch.h"

extern FILE *output;
extern BiqCrunchParameters params;

/*
 * Read the parameters contained in the file given by the argument..
 * @param path: the path to the .param file
 */
void readParameter(char *path) 
{
    FILE * paramfile;
    char s[128];
    char param_name[20];

    // Initialize every parameter with its default value
#define P(type, name, format, def_value)\
    params.name = def_value;
    PARAM_FIELDS;
#undef P

    if (strcmp(path, "") != 0) {
        // open parameter file
        if ((paramfile = fopen(path, "r")) == NULL) {
            fprintf(stderr,
                    "Error: parameter file %s not found.\n", path);
            exit(1);
        }
        printf("Parameter file: %s\n", path);
        while (!feof(paramfile)) {
            if ( fgets(s, 120, paramfile) != NULL ) {
                if (s[0] == '#' || s[0] == '\n' || s[0] == ' ') {
                    // it's a comment, do nothing
                    strcpy(param_name, "\0");
                } else {
                    sscanf(s, "%[^=^ ]", param_name);

                    // Compare all the parameter with the parameter read and change 
                    // the value only if they're the same
#define P(type, name, format, def_value)\
                    if(strcmp(#name, param_name) == 0)\
                    sscanf(s, "%*[^=]="format"\n", &(params.name));
                    PARAM_FIELDS;
#undef P
                }
            }
        }
        fclose(paramfile);
    }

    // Report the values of the parameters used
    fprintf(output, "BiqCrunch Parameters:\n");
#define P(type, name, format, def_value)\
    fprintf(output, "%20s = "format"\n", #name, params.name);
    PARAM_FIELDS;
#undef P
}

