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
 *  File   : ph.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Source file of Pairing Heap
 */

#define PAIRINGHEAP

/* #define PRELOCK */

#include <stdio.h>
#include <math.h>
#include STRINC

VAR_SHARED BobPQ *pqs;

/*------ Swap les deux arbres si besoin est ------*/
IsSwapP(n1,n2)
BobPQNode **n1,**n2;
{
BobPQNode *n3;

  if ( Bob_PRIL((*n1)->an.Pri,(*n2)->an.Pri) ) {
      n3 = *n2;
     *n2 = *n1;
     *n1 =  n3;
  }
}


/*----------------------------------------------------------------*/
void Bob_PQAlloc(NbDS)
int NbDS;
{
int i;

  pqs = (BobPQ *)SH_MALLOC(sizeof(BobPQ)*NbDS);
  if ( pqs ==NULL ) return;
  for (i=0;i<NbDS; i++ ) {
      Bob_PQLINIT(pqs[i]);
      pqs[i].Root = NULL;
  }
}

/*----------------------------------------------------------------*/
Bob_PQFree()
{
  SH_FREE(pqs);
}

/*----------------------------------------------------------------*/
Bob_PQInitLocal(){}



/*------ Pairing ------*/
BobPQNode *Pairing(Root)
BobPQNode *Root;
{
  BobPQNode *n,*Rest;

  if ( Root==NULL ) return NULL;
  if ( Root->ff == NULL ) return Root;

  n = Root->ff;
  Root->ff = NULL;

  Rest = n->ff;

  IsSwapP(&Root,&n);
  
  n->ff = Root->fs;
  Root->fs = n;

  Root->ff = Pairing(Rest);
  if ( Root->ff!=NULL ) {
     if ( Bob_PRIG(Root->ff->an.Pri,Root->an.Pri) ) {
        n = Root->ff; 
        Root->ff = n->ff;
        n->ff = Root;
        Root = n;
     }
  }
  return Root;
}

BobPQNode *Cons(Root)
BobPQNode *Root;
{
BobPQNode *n,*tmp;

  Root = Pairing(Root);

  n = Root->ff;    
  Root->ff = NULL;    
  for ( tmp = Root->fs ; tmp!=NULL && tmp->ff!=NULL ; tmp=tmp->ff);
  if ( tmp==NULL )
     Root->fs = n;
  else 
     tmp->ff = n;
  return Root;
}


BobNode *Bob_PQDel(dsi)
int dsi;
{
BobPQ *pq=pqs+dsi;
BobPQNode *First;

  Bob_EXLON(pq);
  if ( (First = pq->Root)==NULL ) return NULL;
  pq->Root = (pq->Root->fs!= NULL ? Cons(pq->Root->fs):NULL);
  Bob_EXLOFF(pq);
  return &(First->an);
}

void Bob_PQIns(dsi,n)
int dsi;
BobPQNode *n;
{
BobPQ *pq=pqs+dsi;

  n->fs=NULL;
  Bob_EXLON(pq);
  if ( pq->Root!=NULL ) {
     if ( Bob_PRIL(pq->Root->an.Pri,n->an.Pri) ) {
       n->fs = pq->Root;
       n->ff = NULL;
       pq->Root->ff = pq->Root->fs;
       pq->Root->fs = NULL;
       pq->Root = n;
     } else {
       n->ff=pq->Root->fs;
       pq->Root->fs = n;
     }
  }
  else {
     pq->Root = n;
     n->ff = NULL;
  }
  Bob_EXLOFF(pq);
}

/*----------------------------------------------------------------*/
int Bob_PQDelG(dsi,Val)
int dsi;
int Val;
{
BobPQ *pq=pqs+dsi;

   return 0;
}


