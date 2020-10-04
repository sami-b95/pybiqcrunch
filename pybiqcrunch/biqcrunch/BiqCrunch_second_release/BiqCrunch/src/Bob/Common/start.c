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
 *  File   : start.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Source file of startup and configuration functions.
 */

#include <stdio.h>
#include <math.h>
#include "../Include/bb.h"

VAR_SHARED BobNode     *BobRoot;
VAR_SHARED BobSolution *BobSol;
VAR_SHARED int          BobId=0;

#if ARTP==SHARED
VAR_SHARED BobTCt BobCtTab[MAXPROCS];
#else
VAR_SHARED BobTCt BobCt;
#endif

/*------------------------------------------------------------------*/
void Bob_StPrint() {
  
  Bob_PsStPrint(stderr);
  Bob_PrintSolution(BobSol);
  Bob_DsStPrint(stderr);
}


/*------------------------------------------------------------------*/
#if ARTP==SHARED
Bob_DataInit(d)
BobTData *d;
{
  d->Val = 0;
  ARCH_LCKINIT((*d));
}
#endif


/*------------------------------------------------------------------*/
int Bob_CalcMaxDepth()
{
  
#if ARTP!=SEQ
  if ( BobCt.ExpCt.Type&1 ) 
     BobCt.ExpCt.Max = BobId%BobCt.ExpCt.Max;
  if ( BobCt.ExpCt.Type&2 )  {
     if ( BobCt.ExpCt.Max > BobCt.NbProcs ) 
       BobCt.ExpCt.Max = (BobCt.ExpCt.Max/BobCt.NbProcs)*BobId;
     else  
       BobCt.ExpCt.Max = BobId%BobCt.ExpCt.Max;
  }
#endif

}


/*------------------------------------------------------------------*/

int Bob_ExpCtrl(iprof,ExpCt)
int iprof;
BobTExpCt *ExpCt;
{
#if ARTP==SHARED
  return (iprof<ExpCt->Max) && Bob_GPQTestNbNode();
#elif ARTP==DISTRIB || ARTP==MLTHR
  return (iprof<ExpCt->Max) && Bob_GPQTestNbNode();
#else
  return (iprof<ExpCt->Max);
#endif
}


/*--- Tache principale ---*/
void Bob_Algo(id)
int id;
{
  BobTPri Pri;
  BobNode *an;

  BobId = id;

  Bob_CalcMaxDepth(id);

  Bob_PQInitLocal(id);

  if ( BobId==0 ) {
     Pri = BobRoot->Pri;
     Bob_GPQIns(BobRoot);
  }

  while (1) {
    an = Bob_GPQDel();
    if ( an==NULL ) { return;}
    Pri = an->Pri;

    if ( Bob_ULBSup(Bob_PRIEVAL(an->Pri)) ) {
       Bob_GenChild(an,0,&(BobCt.ExpCt));
    } else {
       Bob_NodeFree(an);
       Bob_STNEXP();
    }
  }
}


/*------------------------------------------------------------------*/

Bob_PrsVerbose(v)
char *v;
{
  BobCt.Stat = atoi(v);
}

Bob_PrsDS(s)
char *s;
{
  if ( *s==':' ) {
     s++;
     if ( isdigit(*s) ) {
        BobCt.NbDS = atoi(s);
     } else if ( *s=='O') {
        BobCt.NbDS=-1;
     }
  } else BobCt.NbDS=1;
}


Bob_PrsDepth(v)
char *v;
{
char *s=v;

  BobCt.ExpCt.Max = atoi(v);
  while(*s && *s!=':' ) s++;
  if ( *s==':' ) {
     BobCt.ExpCt.Type = atoi(s+1);
  } else {
     BobCt.ExpCt.Type = 0;
  }
  if ( BobSTAT1 ) 
     printf("Max Depth=%d  Type Depth=%d\n",BobCt.ExpCt.Max,BobCt.ExpCt.Type);
}


Bob_PrsProcs(v)
char *v;
{ char *s=v;
  char *t;

  if ( BobSTAT1 ) printf("Parallel model\n");
  while(*s && *s!=':' ) s++;
  Bob_PrsDS(s);
  *s = 0;
  BobCt.NbProcs= atoi(v);
//  if ( BobSTAT1 )
//     printf("NbProcs=%d  NbDS=%d\n",BobCt.NbProcs,BobCt.NbDS);
}

/*------------------------------------------------------------------*/

Bob_PrsParam(n,v)
int n;
char **v;
{

#if ARTP==SEQ 
char *OptStr="d:v:";
char *UsageStr="Usage : %s -d n<:n> -v <n> ...\n";
#else
char *OptStr="p:d:v:";
char *UsageStr="Usage : %s -p n<:n> -d n<:n> -v <n> ...\n";
#endif
char c;
extern char *optarg;
extern int optind;

   BobCt.NbProcs  = 1;
   BobCt.NbDS     = 1;
   BobCt.PbSize   = 0;
   BobCt.ExpCt.Type= 0;
   BobCt.ExpCt.Max= 0;
   BobCt.Stat     = 0;

   while ((c = getopt(n, v, OptStr)) != GETOPTFAIL ) {
    switch(c) {
      case 'v' : Bob_PrsVerbose(optarg);
                 break;
      case 'd' : Bob_PrsDepth(optarg);
                 break;
#if ARTP!=SEQ 
      case 'p' : Bob_PrsProcs(optarg);
                 break;
#endif
      case '?' : fprintf(stderr,UsageStr,v[0]);
                 break;
    }
   }
   if ( BobSTAT1 ) {
       fprintf(stderr,"Procs:%d:%d  ExplCt:%d:%d Verbose:%d\n",
            BobCt.NbProcs,BobCt.NbDS,BobCt.ExpCt.Max, 
            BobCt.ExpCt.Type, BobCt.Stat);
   }
   return optind;
}

/*------------------------------------------------------------------*/
Bob_CtInit() {
int i;

#if ARTP==SHARED
  for (i=1;i<BobCt.NbProcs; i++ ) { 
    BobCtTab[i]= BobCtTab[0];
  }
#endif
}


/*------------------------------------------------------------------*/

Bob_TimeInit(bt)
BobTTime *bt;
{
  gettimeofday(&(bt->Init),&(bt->TimeZone));
}

Bob_TimeComp(bt)
BobTTime *bt;
{
  gettimeofday(&(bt->Comp),&(bt->TimeZone));
}

Bob_TimeEnd(bt)
BobTTime *bt;
{
  gettimeofday(&(bt->End),&(bt->TimeZone));
  bt->SecComp = bt->End.tv_sec-bt->Comp.tv_sec +
                (float)(bt->End.tv_usec-bt->Comp.tv_usec)/1e6;
  bt->SecInit = bt->Comp.tv_sec-bt->Init.tv_sec +
                (float)(bt->Comp.tv_usec-bt->Init.tv_usec)/1e6;

//  if ( BobSTAT1 ) {
//     printf("---PROCS: %d ---TEMPS : %f %f\n",BobCt.NbProcs,
//            bt->SecComp,bt->SecInit);
//     Bob_StPrint();
//  } else {
//     printf("%d  %f\n",BobCt.NbProcs, bt->SecComp);
//  }
}

/*------------------------------------------------------------------*/
#ifdef NEEDMAIN 
main(n,v)
int n;
char **v;
{
  BobTTime BobTime;

  Bob_TimeInit(&BobTime);

  Bob_Main(n,v,&BobTime);

  exit(0);
}
#endif

