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
 *  File   : dh.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  File   : Source file of Rao and Kumar concurrent Heap.
 */

#define RKHEAP 
#define NBELEM 200000

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include STRINC


#define ABSENT  0
#define PRESENT 1
#define WANTED  2
#define PENDING 3

#define IN(i) (pq->q+i)
#define STATUS(i) pq->q[i].Status
#define VALUE(i)  (pq->q[i].n==NULL? MaxPri : pq->q[i].n->an.Pri)
#define APPNODE(i)  pq->q[i].n
#define LSON(i)   (i*2)
#define RSON(i)   ((i*2)+1)
/*#define MIN(i)    (Bob_PRIG(VALUE(LSON(i)),VALUE(RSON(i))) ? LSON(i) : RSON(i) )
#define MAX(i)    (Bob_PRIG(VALUE(LSON(i)),VALUE(RSON(i))) ? RSON(i) : LSON(i) )*/

VAR_SHARED BobPQ *pqs;

VAR_PRIVATE BobTPri MaxPri;

void Bob_PQAlloc(NbDS)
int NbDS;
{
INTPQNode *n;
INTPQNode *nSup;
int i,j;

  pqs = (BobPQ *)SH_MALLOC(sizeof(BobPQ)*NbDS);
  if ( pqs==NULL ) return;
  for (i=0;i<NbDS;i++) {
     Bob_PQLINIT(pqs+i);
     nSup= &(pqs[i].q[NBELEM-1]);
     for (j=0;j<NBELEM; j++ ) {
        n = &(pqs[i].q[j]);
        Bob_NDLINIT(n);
        n->n = NULL;
        n->Status = ABSENT;
     }
     pqs[i].lastelem=pqs[i].fulllevel=0;
  }
}

Bob_PQFree()
{
  SH_FREE(pqs);
}

Bob_PQInitLocal(){ 
   Bob_PRIINFI(MaxPri);
}


/****** inter-echange suivant priorite *****/
void Exchange(q1,q2)
BobPQNode **q1, **q2;
{
  BobPQNode *dum;
 
  dum = *q1;
  *q1 = *q2;
  *q2 = dum;
}

void Bob_PQIns(dsi,n) 
int dsi;
BobPQNode *n;
{
int target,i,j,k;
BobTPri pri;
BobPQ *pq=pqs+dsi;

   Bob_PQLON(IN(1));
   (pq->lastelem)++; 
   target = pq->lastelem;
   if ( pq->lastelem>=pq->fulllevel*2 ) {
      pq->fulllevel = pq->lastelem;
   } 
   i = target-pq->fulllevel;
   j = pq->fulllevel/2;
   k = 1;
   STATUS(target) = PENDING;
   
   while(j!=0 ) {
      if ( STATUS(target) == WANTED ) break;

      pri = VALUE(k);
      if ( Bob_PRIL(pri,n->an.Pri) ) Exchange(&n,&(APPNODE(k)));

      if ( i>=j ) {
         Bob_PQLON(IN(RSON(k)));
         Bob_PQLOFF(IN(k));
         k = RSON(k);
         i -=j;
      } else {
         Bob_PQLON(IN(LSON(k)));
         Bob_PQLOFF(IN(k));
         k = LSON(k);
      }
      j = j/2;
   }
   if ( STATUS(target)== WANTED ) {
      APPNODE(1) = n;
      STATUS(target)=ABSENT;
      STATUS(1) = PRESENT;
   } else {
      APPNODE(target)=n;
      STATUS(target)=PRESENT;
   }
   Bob_PQLOFF(IN(k));
}


/****** Preemption d'un probleme (un noeud de valeur minimale)  *****/
BobNode *Bob_PQDel(dsi)
int dsi;
{
int i,j,mi,ma;
BobPQ *pq =pqs+dsi;
BobPQNode *n;
BobTPri lp,rp,pi,pmi;

  Bob_PQLON(IN(1));
  if ( (pq->lastelem)==0) {
     Bob_PQLOFF(IN(1));
     return NULL;
  }

  n = APPNODE(1);
  i = 1;
  j = pq->lastelem;
  pq->lastelem--;
  if ( pq->lastelem < pq->fulllevel ) {
      pq->fulllevel = pq->fulllevel/2;
  }
  if ( j==1 ) {
      APPNODE(1) = NULL;
      STATUS(1) = ABSENT;
      Bob_PQLOFF(IN(1));
      return &(n->an);
  }
  Bob_PQLON(IN(j));
  if ( STATUS(j) == PRESENT ) {
      APPNODE(1) = APPNODE(j);
      STATUS(j) = ABSENT;
      APPNODE(j)  = NULL;
  } else {
      STATUS(1) = ABSENT;
      STATUS(j) = WANTED;
  }
  Bob_PQLOFF(IN(j));
  while( STATUS(i)==ABSENT );
   
  Bob_PQLON(IN(LSON(i)));
  Bob_PQLON(IN(RSON(i)));
  
  lp = VALUE(LSON(i));
  rp = VALUE(RSON(i));
  if ( Bob_PRIG(lp,rp) ) {mi = LSON(i);ma=RSON(i);}
  else { mi = RSON(i);ma=LSON(i);}
   
  pi= VALUE(i);
  pmi = VALUE(mi);
  while ( Bob_PRIL(pi,pmi) ) {
      Exchange(&(APPNODE(i)),&(APPNODE(mi)) );
      Bob_PQLOFF(IN(i));
      Bob_PQLOFF(IN(ma));
      i = mi;
      Bob_PQLON(IN(LSON(i)));
      Bob_PQLON(IN(RSON(i)));
      lp = VALUE(LSON(i));
      rp = VALUE(RSON(i));
      if ( Bob_PRIG(lp,rp) ) {mi = LSON(i);ma=RSON(i);}
      else { mi = RSON(i);ma=LSON(i);}
      pi= VALUE(i);
      pmi = VALUE(mi);
  }
  Bob_PQLOFF(IN(i));
  Bob_PQLOFF(IN(LSON(i)));
  Bob_PQLOFF(IN(RSON(i)));
  return &(n->an);
}
  

Bob_PQDelG(dsi,Val) 
int dsi;
int Val;
{
BobPQ *pq=pqs+dsi;

  return 0;
}



