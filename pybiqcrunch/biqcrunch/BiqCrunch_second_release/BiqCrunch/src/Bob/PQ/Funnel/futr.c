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
 *  File   : futr.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995
 *  Comment: Source file of the Funnel-Tree priority queue.
 *           First implementation B. Mans 1990.
 */

#define FUNNELTREE

#include <stdio.h>
#include <math.h>
#include STRINC

VAR_SHARED int fullev;
VAR_SHARED int glilb,gliub;
VAR_PRIVATE int Locilb;
VAR_PRIVATE int Lociub;

VAR_PRIVATE int Capa;

VAR_SHARED BobPQ *pqs;


#if ORTP==MINIMISATION

/*--------------------------------------------------------------------*/
void Bob_PQAlloc(NbDS,Pilb,Piub)
int NbDS,Pilb,Piub;
{
  int i,j,SizeFifo;
  int depth;

  glilb = Pilb; gliub = Piub;
  SizeFifo = gliub+2 - glilb;

  depth = (E_log2(gliub-glilb-1)+1);
  fullev = 1<<depth;
       /*  puissance de 2 superieure au gap */

  pqs = (BobPQ *)SH_MALLOC(sizeof(BobPQ)*NbDS);
  if ( pqs==NULL ) return;
  for(j=0;j<NbDS;j++){
    pqs[j].funnel = (BobFUNode *)SH_MALLOC(fullev*2*sizeof(BobFUNode));
    if ( pqs[j].funnel == NULL ) exit(1);
    for(i=1; i<fullev*2; i++) {
      pqs[j].funnel[i].Val=0;
      Bob_NDLINIT(&(pqs[j].funnel[i]));
    }
    pqs[j].queues = (FIFO *)SH_MALLOC(SizeFifo*sizeof(FIFO));
    if ( pqs[j].queues == NULL ) exit(1);

    for(i=0; i<SizeFifo; i++) {
        pqs[j].queues[i].head=NULL;
        pqs[j].queues[i].tail=NULL;
    }
    Bob_PQLINIT(&pqs[j]);
  }
}

/*--------------------------------------------------------------------*/
Bob_PQFree()
{
  int i;
  for (i=0;i<BobCt.NbDS;i++) {
     SH_FREE(pqs[i].funnel);
     SH_FREE(pqs[i].queues);
  }
  SH_FREE(pqs);
}

/*----------------------------------------------------------------*/
Bob_PQCapa(dsi,nb)
int dsi,nb;
{
  return Capa;
}


/*--------------------------------------------------------------------*/
Bob_PQInitLocal()
{ 

  Locilb = glilb; 
  Lociub = gliub; 
  Capa = (E_log2(Lociub-Locilb-1)+1);
}


/*--------------------------------------------------------------------*/
/*------ ENQUEU ------*/
void enqueu(pq,n,leaf)
BobPQ *pq;
BobPQNode *n;
int leaf;
{
#if PRTP==EVAL
  if (pq->queues[leaf].head==NULL )
    pq->queues[leaf].head = n;
  else
    pq->queues[leaf].tail->Next = n;

  pq->queues[leaf].tail = n;
#elif PRTP==EVDP
BobPQNode *p=pq->queues[leaf].head;
BobPQNode *pp=NULL;
  while(p!=NULL && p->an.Pri.Depth < n->an.Pri.Depth ) {
      pp=p;
      p= p->Next;
  }
  n->Next = p;     
  if ( pp==NULL ) {
    pq->queues[leaf].head = n;     
  } else {
    pp->Next = n;
  }
#endif
}



/*--------------------------------------------------------------------*/
/******* Insert thru the "sBobIde path" ******/
void metfil(pq,p)
BobPQ *pq;
BobPQNode *p;
{
  int leaf, e_targ, target, i, j;
  BobFUNode *funnel = pq->funnel;


  leaf=Bob_PRIEVAL(p->an.Pri)-Locilb;
  e_targ=leaf+fullev;
  target=1;
  i=leaf;
  j=fullev/2;
  Bob_PQLON(pq);
  Bob_NDLON(&funnel[target]);
  Bob_PQLOFF(pq);
  ++funnel[target].Val;
  while (e_targ>target) 
    {
      target*=2;
      if (i>=j)
	{
	  target++;
	  i-=j;
	}
      j/=2;
      Bob_NDLON(&funnel[target]);
      ++funnel[target].Val;
      Bob_NDLOFF(&funnel[target/2]);
    }
  enqueu(pq,p,leaf);
  Bob_NDLOFF(&funnel[target]);
}


/*--------------------------------------------------------------------*/
/*------ RANGER ------*/
void Bob_PQIns(dsi,n)
int dsi;
BobPQNode *n;
{
  n->Next = NULL;
  metfil(pqs+dsi,n);
}


/*--------------------------------------------------------------------*/
/*------ DEQUEU ------*/
BobPQNode *dequeu(pq,leaf)
BobPQ *pq;
int leaf;
{
BobPQNode *n;


  n=pq->queues[leaf].head;
  if ( pq->queues[leaf].head==NULL )
    return NULL;

  pq->queues[leaf].head = n->Next;
  if ( pq->queues[leaf].tail == n )
     pq->queues[leaf].tail = NULL;

  return n;
}


/*--------------------------------------------------------------------*/
BobNode *Bob_PQDel(dsi)
int dsi;
{
  int p, leaf, target;
  BobPQ *pq= pqs+dsi;
  BobPQNode *n;
  BobFUNode *funnel;

  Bob_PQLON(pq);
  Bob_NDLON((pq->funnel)+1);
  if ( pq->funnel[1].Val==0 ) {
    Bob_PQLOFF(pq);
    Bob_NDLOFF((pq->funnel)+1);
    return NULL;
  }
  
  funnel = pq->funnel;
  Bob_PQLOFF(pq);
  target=1;
  --funnel[target].Val;
  while (target<fullev) {
	  target*=2;
	  Bob_NDLON(&funnel[target]);
	  if (funnel[target].Val==0) {
	      Bob_NDLOFF(&funnel[target]);
	      ++target;
	      Bob_NDLON(&funnel[target]);
	  }
	  Bob_NDLOFF(&funnel[target/2]);
	  --funnel[target].Val;
  }
  leaf=target-fullev;
  n=dequeu(pq,leaf);
  Bob_NDLOFF(&funnel[target]);
  return(&(n->an));
}
  


/*--------------------------------------------------------------------*/
int Bob_PQDelG(dsi,valeur)
int dsi;
int valeur;
{
  int target, e_targ, leaf, nbelim, i, j, p;
  int SizeFifo;
  BobPQ *pq=pqs+dsi;
  BobPQNode *n;
  BobFUNode *funnel = pq->funnel;

  leaf=valeur-Locilb+1;
  e_targ=leaf+fullev;


/* Top-down Bob_NDLON of the whole path */
  target=1;
  i=leaf;
  j=fullev/2;
  nbelim=0;
  Bob_PQLON(pq);
  Bob_NDLON(&funnel[1]);
  Bob_PQLOFF(pq);
  while (target<e_targ)
    {
      target*=2;
      if (i>=j)
	{
	  ++target;
	  i-=j;
	}
      else
	{
	  Bob_NDLON(&funnel[target+1]);
	  nbelim+=funnel[target+1].Val;
	}
      j/=2;
      Bob_NDLON(&funnel[target]);
    }
  nbelim+=funnel[target].Val;


/* Top-down update-Bob_NDLOFF of the whole path */
  target=1;
  i=leaf;
  j=fullev/2;
  funnel[target].Val -= nbelim;

  while (target<e_targ)
    {
      target*=2;
      if (i>=j)
	{
	  ++target;
	  i-=j;
	}
      else
	{
	  nbelim-=funnel[target+1].Val;
	  funnel[target+1].Val=0;
	  Bob_NDLOFF(&funnel[target+1]);
	}      
      Bob_NDLOFF(&funnel[target/2]);
      j/=2;
      funnel[target].Val-=nbelim;
    }
  Bob_NDLOFF(&funnel[target]);
  SizeFifo = Lociub+1 - Locilb;

  j = 0;
  for(i=leaf; i<SizeFifo ; i++ ) {
     while ((n=dequeu(pq,i))!=NULL ) {Bob_NodeFree(n);j++;}
  }
  return j;
}

/*--------------------------------------------------------------------*/
int E_log2(p)
int p;
{
  unsigned  i;

  for(i=0; p!=0; i++, p>>=1);
  if (i)
    return(i-1);
  else
    return(0);
}

#endif
