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

// Convert Unconstrained Binary Quadratic Problems to BiqCrunch instances
// 
// 2012 Marco Casazza
// Part of BiqCrunch project : http://lipn.univ-paris13.fr/BiqCrunch/

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    FILE    *fid;
    int     ret;
    int     n, m, line;
    int     i, j, val;
    double  *Q;
    int     N;

    /*
       Check number of arguments
       */
    if(argc != 2) { 
        printf("Usage: qp2bc inputfile > outputfile\n");
        exit(1);
    }

    /*
       Open file
       */
    fid = fopen(argv[1],"r");
    if (fid == (FILE *) NULL) {
        printf("qp2bc: Couldn't open problem file for reading! \n");
        exit(1);
    }

    /*
       Read the number of nodes and edges from the first line of the file
       */
    ret = fscanf(fid,"%d %d",&n,&m);
    if (ret != 2) {
        printf("qp2bc: Incorrect input file. First line has to give n and m values.\n");
        fclose(fid);
        exit(1);
    }

    /*
       Initialize the objective matrix Q
       */
    N = n+1;
    Q = (double *) malloc(N*N*sizeof(double));
    for (i=0;i<N;i++)
        for (j=0;j<N;j++)
            Q[i+j*N] = 0;

    /*
       Read the matrix from the file and compute the objective function matrix for BiqCrunch
       */
    for (line=1;line<=m;line++) {

        ret = fscanf(fid,"%d %d %d",&i,&j,&val);

        if (ret != 3)
        {
            printf("qp2bc: Incorrect input file. Can't read values.\n");
            free(Q);
            fclose(fid);
            exit(1);
        }

        // if the value is on the diagonal of the matrix is moved on the last column and last row
        if(i == j){
            Q[(i-1)+(N-1)*N] = (double)val/2.0;
            Q[(N-1)+(i-1)*N] = (double)val/2.0;
        }else{
            Q[(i-1)+(j-1)*N] = (double)val;
            Q[(j-1)+(i-1)*N] = (double)val;
        }
    }
    fclose(fid);


    // Print the instance
    // # constraints, # blocks, size 
    printf("1 =max problem\n");
    printf("0\n1\n%d\n",N);
    // Q matrix    
    for (i=0;i<N;i++) {
        for (j=i+1;j<N;j++) {
            if (Q[i+j*N] != 0) 
                printf("0 1 %d %d %lf\n",i+1,j+1,Q[i+j*N]);
        }
    }

    free(Q);
    exit(0);
}
