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

////////////////////////////////////////////////////////////////////////////////
//									      //
//									      //
//									      //
//	Convert rudy files and dats file into BiqCrunch problem file	      //
//	for the k-cluster problem.					      //
//									      //
//									      //
// 	2012: Frederic Roupin (LIPN), Geoffrey Kozak			      //
// 	Part of BiqCrunch project : http://lipn.univ-paris13.fr/BiqCrunch/    //
//									      //
//									      //
//									      //
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TRUE 1
#define FALSE 0

//Compute adjacency matrix of the problem in file fileGraphe
double* getAdjacencyMatrix(char *fileGraphe);

//Main function to convert a graph instance into a BC file
void toBC(char *fileGraphe, char* outputname, int k, int weighted_graph);

//Convert dat file into the standard format use by rudy
int toStdFormat(char* fileGraphe);

//Check if the vertices of the graph are numbered from 0 to (n-1) or from 1 to n
int firstNodeIsZero(char* fileGraphe,char* inputfile);

void toBC(char *fileGraphe, char* outputname, int k, int weighted_graph)
{
    FILE *f_graphe, *f_out;
    int nb_node=0, nb_edge=0,i,j;

    f_out = fopen(outputname, "w");
    if (f_out != NULL)
    {
        f_graphe = fopen(fileGraphe, "r");
        if (f_graphe != NULL)
        {
            fscanf(f_graphe, "%d %d\n",&nb_node, &nb_edge);
            fclose(f_graphe);
        }
        else
        {
            perror("Reading error !!!\n");
            return;
        }

        fprintf(f_out,"1 =max problem\n");
        fprintf(f_out,"%d =mdim\n1 =nblocks\n{%d}\n",nb_node+1,nb_node+1);
        fprintf(f_out,"%lf ",(double)(2*k));
        for(i=0;i<nb_node;i++)
            fprintf(f_out,"%lf ",0.);
        fprintf(f_out,"\n");

        double *adjMatrix = getAdjacencyMatrix(fileGraphe);

        for(i=0 ; i < nb_node ; i++)
            for(j=i ; j<nb_node ; j++)
                if(adjMatrix[i*nb_node+j])
                {
                    if(weighted_graph)
                        fprintf(f_out,"0 1 %d %d %lf\n",i+1,j+1,(adjMatrix[i*nb_node+j]/2));
                    else
                        fprintf(f_out,"0 1 %d %d %lf\n",i+1,j+1,(1./2.));						
                }

        for(i=0 ; i < nb_node ; i++)
            fprintf(f_out,"1 1 %d %d %lf\n",i+1,nb_node+1,1.);

        for(i=0 ; i < nb_node ; i++)
        {
            for(j=0 ; j<nb_node ; j++)
            {
                if(i<j)
                    fprintf(f_out,"%d 1 %d %d %lf\n",i+2,i+1,j+1,1.);
                else if(i==j)
                    fprintf(f_out,"%d 1 %d %d %lf\n",i+2,i+1,i+1,2.);
                else
                    fprintf(f_out,"%d 1 %d %d %lf\n",i+2,j+1,i+1,1.);
            }
            fprintf(f_out,"%d 1 %d %d %lf\n",i+2,i+1,nb_node+1,-(double)k);
        }

        free(adjMatrix);

        fclose(f_out);
    }
    else
    {
        perror("Writting error !!!\n");
        return;
    }
}

int firstNodeIsZero(char* fileGraphe,char* inputfile)
{
    FILE *f_graphe;
    int nb_nodes,k,nb_edges,head,tail;
    double weight;

    f_graphe = fopen(fileGraphe, "r");
    if (f_graphe != NULL)
    {
        if(strcmp("txt",inputfile)==0)
            fscanf(f_graphe, "%d %d\n",&nb_nodes, &nb_edges);
        else if(strcmp("dat",inputfile)==0)
            fscanf(f_graphe, "%d %d %d\n",&nb_nodes, &k ,&nb_edges);

        while(fscanf(f_graphe,"%d %d %lf\n",&head, &tail,&weight)>0)
            if(tail==0 || head==0)
                return 1;

        fclose(f_graphe);
        return 0;
    }
    else
        return -1;
}

int toStdFormat(char* fileGraphe)
{
    FILE *f_graphe, *f_out;
    int nb_nodes, k, nb_edges,head,tail;
    double weight;
    int beginAtZero = firstNodeIsZero(fileGraphe,"dat");

    f_graphe = fopen(fileGraphe, "r");
    f_out = fopen("output.txt", "w");
    if (f_graphe != NULL)
    {
        if (f_out != NULL)
        {
            fscanf(f_graphe, "%d %d %d\n",&nb_nodes, &k ,&nb_edges);
            fprintf(f_out,"%d %d\n",nb_nodes,nb_edges/2);

            while(fscanf(f_graphe,"%d %d %lf\n",&head, &tail,&weight)>0)
            {
                if(tail>head && beginAtZero==1)
                    fprintf(f_out,"%d %d %lf\n",head+1,tail+1,weight*2);
                else if(tail>head && beginAtZero==0)
                    fprintf(f_out,"%d %d %lf\n",head,tail,weight*2);
            }
            fclose(f_out);
            return k;
        }	
        else
        {
            perror("Reading error !!!\n");
            return -1;
        }	

        fclose(f_graphe);
    }
    else
        perror("Reading error !!!\n");

    return -1;
}

double* getAdjacencyMatrix(char *fileGraphe)			//Compute the adjacency matrix of graph describe in the input file
{								//Note : it will build a good adjacency matrix if the input file have the same structure 
    //as the one use by rudy (graph generator made by Giovanni Rinaldi)
    int N=0,head,tail,i,j,nb_node=0,nb_edge=0;		
    double weight;
    double *Matrix = NULL;

    FILE *f_graphe;

    f_graphe = fopen(fileGraphe, "r");				//Open the file in reading mode

    if (f_graphe != NULL)						//If we succeed to open the file
    {
        fscanf(f_graphe, "%d %d\n",&nb_node, &nb_edge);		//We read the first line and we get the number of vertices(nb_nodes) and number of edge
        N=nb_node;

        Matrix=(double *) malloc(N * N * sizeof(double));	//Create the adjacency matrix,with a dynamic memory allocation	

        for(i=0 ; i < N ; i++)					//Initialization of the adjacency matrix
            for(j=0 ; j<N ; j++)
                Matrix[i*N+j]=0.;

        while(fscanf(f_graphe,"%d %d %lf\n",&head, &tail,&weight)>0)	//While we can read the file, we collect edge one by one
        {		
            Matrix[(head-1)*N+(tail-1)]=weight;			//We add the edge (head,tail) in the adjacency matrix
            Matrix[(tail-1)*N+(head-1)]=weight;			//complete the symmetrical part
        }
        fclose(f_graphe);						//close the file
    }
    else									//if we failed to open the file
        perror("Reading error !!!\n");					//error message

    return Matrix;								//return adjacency matrix compute
}


int main(int argc, char** argv)
{ 
    char inputfile[40], kval[10];
    char *ext,*outputname;
    int k;
    int weighted_graph;

    if(argc==1)
    {
        printf("\n\n[Use :] - ./KC2BC <inputfile.txt> <k> -w [0|1]\n\t=> Format .txt does not contain the k value\n\tso you have to specify it in the second parameter.\n\n\t- ./KC2BC <inputfile.dat> -w [0|1]\n\t=> Format .dat contains the k value, so you don't have to give it to the generator.\n\nOption \"-w\" : \n\t-> followed by 0 if you want to build the associated graph with\n\tunweighted edges in the output BC file.\n\t-> followed by 1 if you want to build the edge-weighted graph\n\tin the output BC file.\n\n");
        return -1;
    }

    strcpy (inputfile,argv[1]);
    outputname= strtok(inputfile,".");
    ext= strtok(NULL,".");


    if(strcmp("dat",ext)==0)
    {
        if(argc==4){
            if(strcmp("-w",argv[2])!=0)
            {
                printf("\n%s : Unknown option !\n",argv[2]);
                return -1;
            }
            weighted_graph=atoi(argv[3]);
            if(weighted_graph!=0 && weighted_graph!=1)
            {
                printf("\n\n[Option \"-w\"]\n\t-> followed by 0 if you want to build the associated graph with\n\tunweighted edges in the output BC file.\n\t-> followed by 1 if you want to build the edge-weighted graph\n\tin the output BC file.\n\n");
                return -1;
            }
            k=toStdFormat(argv[1]);
            sprintf(kval,"%d",k);
            outputname=strcat(outputname,"_k");
            outputname=strcat(outputname,kval);
            outputname=strcat(outputname,".bc");
            if(weighted_graph)
                toBC("output.txt",outputname,k,TRUE);
            else
                toBC("output.txt",outputname,k,FALSE);
            system("rm output.txt");
        }
        else if(argc<4){
            printf("\nSome parameters are missing !\n");
            return -1;
        }
        else{
            printf("\nToo many arguments !\n");
            return -1;
        }
    }
    else
    {
        if(argc==5)
        {
            if(strcmp("-w",argv[3])!=0)
            {
                printf("\n\n%s : Unknown option !!\n\n",argv[3]);
                return -1;
            }
            weighted_graph=atoi(argv[4]);
            if(weighted_graph!=0 && weighted_graph!=1)
            {
                printf("\n\n[Option \"-w\"]\n\t-> followed by 0 if you want to build the associated graph with\n\tunweighted edges in the output BC file.\n\t-> followed by 1 if you want to build the edge-weighted graph\n\tin the output BC file.\n\n");
                return -1;
            }
            k=atoi(argv[2]);
            outputname=strcat(outputname,"_k");
            outputname=strcat(outputname,argv[2]);
            outputname=strcat(outputname,".bc");
            if(weighted_graph)
                toBC(argv[1],outputname,k,TRUE);
            else
                toBC(argv[1],outputname,k,FALSE);
        }
        else if(argc<5){
            printf("\nSome parameters are missing !\n");
            return -1;
        }
        else{
            printf("\nToo many arguments !\n");
            return -1;
        }
    }

    return 0;
}
