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
 *  File   : futb.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995
 *  Comment: Source file of the Funnel-Table priority queue.
 */

#define FUNNELTABLE

#include <stdio.h>
#include <math.h>
#include STRINC


VAR_SHARED int gl_SizeFun;
VAR_PRIVATE int SizeFun,Locilb;
VAR_SHARED int glilb,gliub;

VAR_PRIVATE int    LastVal=0;
VAR_SHARED  BobPQ *pqs;

#if ORTP==MINIMISATION
/*--------------------------------------------------------------------*/
void Bob_PQAlloc(NbDS,Pilb,Piub)
{
  int i,j;

  glilb = Pilb;
  gliub = Piub;
  gl_SizeFun = gliub+2 - glilb;
  SizeFun = gl_SizeFun;

  if ( NbDS!=1 ) {fprintf(stderr,"Multi-Structure not supported\n");exit(1);}
  pqs = (BobPQ *)SH_MALLOC(sizeof(BobPQ)*NbDS);
  if ( pqs==NULL ) {fprintf(stderr,"Not enough Memory\n"); return;}
  for(j=0;j<NbDS; j++ ){
     pqs[j].funnel = (FIFO *)SH_MALLOC((SizeFun+1)*sizeof(FIFO));
     if ( pqs[j].funnel==NULL ) {fprintf(stderr,"Not enough Memory\n"); return;}

     for ( i=0 ; i<SizeFun+1 ; i++ ) {
       pqs[j].funnel[i].head = NULL;
       pqs[j].funnel[i].tail = NULL;
       pqs[j].funnel[i].NbProc = 0;
       Bob_PQLINIT( &(pqs[j].funnel[i]));
     }
     pqs[j].binf = SizeFun;
     Bob_PQLINIT(&(pqs[j]));
  }
  pqs[0].binf = 0;
  pqs[0].funnel[0].NbProc = 1;
  LastVal = Pilb;
}

/*--------------------------------------------------------------------*/
Bob_PQFree() 
{
int i;
  for (i=0;i<BobCt.NbDS;i++)
     SH_FREE(pqs[i].funnel);
  SH_FREE(pqs);
}

/*--------------------------------------------------------------------*/
Bob_PQInitLocal()
{
  SizeFun = gl_SizeFun;
  Locilb = glilb;
}


/*--------------------------------------------------------------------*/
/*------ ENQUEU ------*/
void enqueu(pq,n,leaf)
BobPQ *pq;
BobPQNode *n;
int leaf;
{
FIFO *funnel = pq->funnel;
#if PRTP==EVAL

  if (funnel[leaf].head==NULL )
    funnel[leaf].head = n;
  else
    funnel[leaf].tail->Next = n;

  funnel[leaf].tail = n;

#elif PRTP==EVDP
BobPQNode *p=funnel[leaf].head;
BobPQNode *pp=NULL;
  while(p!=NULL && p->an.Pri.Depth < n->an.Pri.Depth ) {
      pp=p;
      p= p->Next;
  }
  n->Next = p;     
  if ( pp==NULL ) {
    funnel[leaf].head = n;     
  } else {
    pp->Next = n;
  }
#endif
}




/*--------------------------------------------------------------------*/
/*------ DEQUEU ------*/
BobPQNode *dequeu(pq,leaf)
BobPQ *pq;
int leaf;
{
BobPQNode *n;
FIFO *funnel = pq->funnel;
  
  n=funnel[leaf].head;
  if ( funnel[leaf].head==NULL )
    return NULL;

  funnel[leaf].head = n->Next;
  if ( funnel[leaf].tail == n )
     funnel[leaf].tail = NULL;
  return n;
}


/*--------------------------------------------------------------------*/
/*------ RANGER ------*/
void Bob_PQIns(dsi,n)
int dsi;
BobPQNode *n;
{
BobTPri Val;
FIFO *funnel;
BobPQ *pq=pqs+dsi;
  
  Val = n->an.Pri;
  funnel = pq->funnel;
  n->Next = NULL;
  Bob_PQLON(&(funnel[Val-Locilb]));
  enqueu(pq,n,Bob_PRIEVAL(Val)-Locilb);
  Bob_PQLOFF(&(funnel[Val-Locilb]));
  if ( BobCt.NbDS>1 ) {
     Bob_PQLON(pq);
     if ( pq->binf > Bob_PRIEVAL(Val)-Locilb ) {
        pq->binf = Bob_PRIEVAL(Val)-Locilb;
     }
     Bob_PQLOFF(pq);
  }

}


/*--------------------------------------------------------------------*/
/*------ DEMPB ------*/
BobNode *Bob_PQDel(dsi)
int dsi;
{
BobPQNode *n=NULL;
int Ind,Val;
FIFO *funnel;
BobPQ *pq=pqs+dsi;

  funnel = pq->funnel;
  Bob_PQLON(&(funnel[LastVal-Locilb]));
  pq->funnel[LastVal-Locilb].NbProc--;
  Bob_PQLOFF(&(funnel[LastVal-Locilb]));
  Bob_PQLON(pq);
  Bob_PQLON(&(funnel[pq->binf]));
  while ( funnel[pq->binf].head==NULL && 
          funnel[pq->binf].NbProc ==0 &&
          pq->binf<SizeFun) {
        Bob_PQLON(&(funnel[pq->binf+1]));
        Bob_PQLOFF(&(funnel[pq->binf]));
        pq->binf++;
  }
  Ind = pq->binf;
  Bob_PQLOFF(pq);
  while ( 1 ) {
     while( Ind<SizeFun && funnel[Ind].head==NULL ) {
        Bob_PQLON(&(funnel[Ind+1]));
        Bob_PQLOFF(&(funnel[Ind]));
        Ind++;
     }
     if ( Ind<SizeFun ) {
        funnel[Ind].NbProc++;
        n = dequeu(pq,Ind);
        Bob_PQLOFF(&(funnel[Ind]));
        if ( n!=NULL ) break;
     } else {
        Bob_PQLOFF(&(funnel[Ind]));
        break;
     } 
     Bob_PQLON(pq);
     Ind = pq->binf;
     Bob_PQLOFF(pq);
     Bob_PQLON(&(funnel[Ind]));
  }
  if ( n==NULL ) {
      return NULL;
  }
  LastVal = Bob_PRIEVAL(n->an.Pri);
  return &(n->an);
}

/*--------------------------------------------------------------------*/
int Bob_PQDelG(dsi,Val)
int dsi;
int Val;
{
int i;
int cnt=0;
BobPQNode *n;
BobPQ *pq=pqs+dsi;
FIFO *funnel = pq->funnel;

   for (i=Val+1-Locilb ; i<SizeFun ; i++ ) {
      Bob_PQLON(&(funnel[i]));
      n = dequeu(pq,i);
      while( n!=NULL ) {
           Bob_NodeFree(n);
           cnt++;
           n = dequeu(pq,i);
      }
      Bob_PQLOFF(&(funnel[i]));
   }
   return cnt;
}

#endif


