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
 *  File   : skh.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995
 *  Comment: Source file of Skew-Heap.
 *           First implementation B. Mans 1990.
 */

#define SKEWHEAP

/* #define PRELOCK */

#include <stdio.h>
#include <math.h>
#include STRINC

VAR_SHARED BobPQ *pqs;

/*----------------------------------------------------------------*/
void Bob_PQAlloc(nb)
int nb;
{
int i;
int j;
 
  pqs = (BobPQ *)SH_MALLOC(sizeof(BobPQ)*nb);
  if ( pqs==NULL ) return;
  for (i=0;i<nb;i++) {
    Bob_PQLINIT(pqs+i);
    (pqs+i)->Root = NULL;
  }
}

/*----------------------------------------------------------------*/
Bob_PQFree()
{ SH_FREE(pqs);
}

BobNode *Bob_GPQSearch() {}

/*----------------------------------------------------------------*/
Bob_PQCapa(dsi,nb)
int dsi;
int nb;
{
  int j;
  for (j=-1; nb ; nb=nb>>1 ) j++;
  return j/2;
}

/*----------------------------------------------------------------*/
Bob_PQInitLocal(i) { }

/*----------------------------------------------------------------*/
void swap(q1,q2)
BobPQNode **q1, **q2;
{
  BobPQNode *dum;
 
  if ( q2 == NULL ) return; 
  if (Bob_PRIG((*q2)->an.Pri,(*q1)->an.Pri) )
    {
      dum = *q1;
      *q1 = *q2;
      *q2 = dum;
    }
}

/*----------------------------------------------------------------*/
void meld(pq,q1,q2)
BobPQ *pq;
BobPQNode *q1,*q2;
{
  BobPQNode *p;
  
  if (q1==NULL)
    {
      pq->Root=q2;
      Bob_PQLOFF(pq);
      return;
    }
  if (q2==NULL)
    {
      pq->Root=q1;
      Bob_PQLOFF(pq);
      return;
    }

#ifdef PRELOCK
  Bob_NDLON(q1);
  Bob_NDLON(q2);
#endif
  swap(&q1, &q2);
  pq->Root=q1;
#ifndef PRELOCK
  Bob_NDLON(q1);
#endif
  Bob_PQLOFF(pq);
  p=q1;
  q1=p->right;
  p->right=p->left;

  while (q1!=NULL)
    {
#ifdef PRELOCK
      Bob_NDLON(q1);
#endif
      swap(&q1, &q2);
      p->left=q1;
#ifndef PRELOCK
      Bob_NDLON(q1);
#endif
      Bob_NDLOFF(p);
      p=q1;
      q1=p->right;
      p->right=p->left;
    }

  p->left=q2;
#ifdef PRELOCK
  Bob_NDLOFF(q2);
#endif
  Bob_NDLOFF(p);

}


/*----------------------------------------------------------------*/
void Bob_PQIns(dsi,p)
int dsi;
BobPQNode *p;
{
BobPQ *pq=pqs+dsi;

  p->left=NULL;
  p->right=NULL;
  Bob_NDLINIT(p);
      
  Bob_PQLON(pq);
  meld(pq,p,pq->Root);

}


/*----------------------------------------------------------------*/
BobNode *Bob_PQDel(dsi)
int dsi;
{
  BobPQ *pq=pqs+dsi;
  BobPQNode *p;
  int Val;

  Bob_PQLON(pq);
  if ( (p=pq->Root)==NULL ) { 
     Bob_PQLOFF(pq);
     return NULL;
  }
  Bob_NDLON(p);
  Bob_NDLOFF(p);
  meld(pq,p->right,p->left);
  return(&(p->an));
}
  

/*----------------------------------------------------------------*/
int Bob_PQDelG(dsi) 
int dsi;
{
   return 0;
}

