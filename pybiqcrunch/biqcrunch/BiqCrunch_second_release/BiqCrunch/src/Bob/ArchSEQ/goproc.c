/*
 *        BOB version 1.0:  Branch and Bound Optimization LiBrary
 *                    PNN Team of PRiSM laboratory
 *             University of Versailles St-Quentin en Yvelines.
 *      Authors:  M. Benaichouche, V. Cung, S. Dowaji, B. Le Cun
 *                      T. Mautor, C. Roucairol.
 *                    (C) 1995 All Rights Reserved
 *
 *                              NOTICE
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted
 * provided that the above copyright notice appear in all copies and
 * that both the copyright notice and this permission notice appear in
 * supporting documentation.
 *
 * Neither the institutions (Versailles University, PRiSM Laboratory, 
 * the PNN Team), nor the Authors make any representations about the 
 * suitability of this software for any purpose.  This software is 
 * provided ``as is'' without express or implied warranty.
 *
 */

/* 
 *  File   : goproc.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: initializations for ARTP=SEQ.
 */

#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include "../Include/bb.h"

BobTPsSt  BobPsSt;

/*--------------------------------------------------------------*/
void Bob_PsStInit() {

  BobPsSt.Evl=0;
  BobPsSt.NoExp=0;
  BobPsSt.Ins=0;
  BobPsSt.Del=0;
  BobPsSt.Delg=0;
  Bob_PRINULL(BobPsSt.LastPri);
}

Bob_STEVL(){ BobPsSt.Evl++;}
int Bob_EVL(){ return BobPsSt.Evl;}
Bob_STNEXP(){ BobPsSt.NoExp++;}
int Bob_NEXP(){ return BobPsSt.NoExp;}

/*--------------------------------------------------------------*/


void Bob_PsStPrint(io)
FILE *io;
{
  fprintf(io,"--PROCESS ---Stat:----Bob_UpBd:%d\n",Bob_ULBGet());
  fprintf(io,"        Evl   NoExp     Ins     Dem    Delg ");
  Bob_PRISTR(io); fprintf(io,"\n");
  fprintf(io,"     %6d  %6d  %6d  %6d  %6d ",
           BobPsSt.Evl,BobPsSt.NoExp,BobPsSt.Ins,BobPsSt.Del,
           BobPsSt.Delg);
  Bob_PRIPRINT(io,BobPsSt.LastPri);
  fprintf(io,"\n");
}


/*--------------------- Data Structure STAT --------------------*/
/*--------------------------------------------------------------*/

BobTDsSt dst;

/*--------------------------------------------------------------*/

Bob_DsStIns() {
    dst.NbNode++;
    dst.MaxNode = ( dst.MaxNode<dst.NbNode ?
                    dst.NbNode : dst.MaxNode );
}

/*--------------------------------------------------------------*/

Bob_DsStDel(Pri)
BobTPri Pri;
{
    dst.LastPri=Pri;
    dst.NbNode--;
}
/*--------------------------------------------------------------*/

void Bob_DsStInit() {
int i;

  dst.NbNode  = 0;
  dst.MaxNode = 0;
  Bob_PRINULL(dst.LastPri);
}

/*--------------------------------------------------------------*/

Bob_DsStPrint(io)
FILE *io;
{
  int i=0;

  fprintf(io,"-- DSTR STAT ----------\n");
  fprintf(io,"  PQ  NbNode MaxNode ");
  Bob_PRISTR(io); fprintf(io,"\n");
  fprintf(io,"  %2d  %6d  %6d ",i,dst.NbNode,
          dst.MaxNode);
  Bob_PRIPRINT(io,dst.LastPri);
  fprintf(io,"\n");
}



/*--------------------------------------------------------------*/
Bob_Main(n,v,bt)
int n;
char **v;
BobTTime *bt;
{
  char **Para;
  int NbPara,ParaInd;


  ParaInd = Bob_PrsParam(n,v);
  Para = v+ParaInd;
  NbPara = n-ParaInd;

  Bob_PsStInit();
  Bob_DsStInit();

  Bob_Init(NbPara,Para);

//  if ( BobSTAT1 )
//     printf("Lower Bound :%d  Upper Bound :%d\n",
//             Bob_PRIEVAL(BobRoot->Pri),Bob_ULBGet());

  Bob_GPQAlloc();
  Bob_CtInit();

  Bob_TimeComp(bt);

  signal(SIGINT,Bob_StPrint);

  Bob_GoJob(Bob_Algo,v[0],&(BobCt.NbProcs),v+1);

  Bob_TimeEnd(bt);

  Bob_End();

  Bob_GPQFree();

}

/*--------------------------------------------------------------*/
int Bob_GoJob(fonct,binary,nbproc,v) 
int *nbproc;
char *binary;
int (*fonct)();
char **v;
{

  fonct(0);
}

