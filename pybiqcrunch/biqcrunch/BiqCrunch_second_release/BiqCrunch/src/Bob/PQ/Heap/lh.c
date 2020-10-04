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
 *  File   : lh.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995
 *  Comment: Source file of leftist heap.
 */

#define LEFTISTHEAP

#include <stdio.h>
#include <math.h>
#include STRINC

VAR_SHARED BobPQ *pqs;

#define Key(n) (n==NULL ? 0 : n->an.Pri)

Rank(n) 
BobPQNode *n;
{
return ((n)==NULL ? 0 : (n)->Rank);
}

/*----------------------------------------------------------------*/
void Bob_PQAlloc(NbDS)
int NbDS;
{
int i;
  
 pqs = (BobPQ *)SH_MALLOC(sizeof(BobPQ)*NbDS);
 if ( pqs ==NULL ) return;
 for (i=0;i<NbDS;i++) {
   pqs[i].Root = NULL;
   Bob_PQLINIT(pqs+i);
 }
}

/*----------------------------------------------------------------*/
Bob_PQFree()
{ 
   SH_FREE(pqs);
}

/*----------------------------------------------------------------*/
Bob_PQInitLocal() {}


/*----------------------------------------------------------------*/
BobPQNode *Mesh(i1,i2)
BobPQNode *i1,*i2;
{
BobPQNode *tmp;

   if ( i2!=NULL && Bob_PRIL(i1->an.Pri,i2->an.Pri) ) {
      tmp = i1;
      i1 = i2;
      i2 = tmp;
   } 
   i1->fd = ( i1->fd==NULL ? i2 : Mesh( i1->fd , i2 ));
   if ( Rank(i1->fg)<Rank(i1->fd) ) {
      tmp = i1->fg;
      i1->fg = i1->fd;
      i1->fd = tmp;
   }
   i1->Rank = Rank(i1->fd) + 1;
   return i1;
}


/*----------------------------------------------------------------*/
BobPQNode *Meld(i1,i2)
BobPQNode *i1,*i2;
{
   if ( i1==NULL ) return i2;
   if ( i2==NULL ) return i1;
   return Mesh(i1,i2);
}


/*----------------------------------------------------------------*/
BobNode *Bob_PQDel(dsi)
int dsi;
{
  BobPQNode *n;
  BobPQ *pq=pqs+dsi;

  Bob_PQLON(pq);
  n = pq->Root;
  if ( n!=NULL ) pq->Root = Meld(n->fg,n->fd);
  Bob_PQLOFF(pq);
  return &(n->an);
}

/*----------------------------------------------------------------*/
void Bob_PQIns(dsi,n)
int dsi;
BobPQNode *n;
{
  BobPQ *pq=pqs+dsi;

  n->Rank = 1;
  n->fg=NULL; n->fd=NULL;
  Bob_PQLON(pq);
  pq->Root = Meld(n,pq->Root);
  Bob_PQLOFF(pq);
}

/*----------------------------------------------------------------*/
int Bob_PQDelG(dsi,Val)
int dsi;
int Val;
{
   return 0;
}

