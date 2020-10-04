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
 *  File   : Sp.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Source file of Semi-Splay priority queue.
 *           Sequentiel D. Browser.
 */

#define SPLAY


#include <stdio.h>
#include <math.h>
#include STRINC

VAR_SHARED BobPQ *pqs;

/*--------------------------------------------------------------------*/
void Bob_PQAlloc(NbDS)
{
int i;

  pqs = (BobPQ *)SH_MALLOC(sizeof(BobPQ)*NbDS);
  if ( pqs == NULL ) return;
  for (i=0;i<NbDS; i++ ) {
     pqs[i].Root = NULL;
     Bob_PQLINIT(pqs+i);
  }
}

void Bob_PQFree()
{
  SH_FREE(pqs);
}

Bob_PQInitLocal(){ }

/*--------------------------------------------------------------------*/
void Bob_PQIns(dsi,n)
int dsi;
BobPQNode *n;
{
BobPQ *pq=pqs+dsi;
BobPQNode *Top;          /* Noeusd de raccrochement */
BobPQNode Dummy;         /* Reference de Left et Right */
BobPQNode *Cas[2];       /* Left et Right */
BobPQNode *Next;         /* 1er Noeud */
BobPQNode *FarNext;      /* 2ieme Noeud */
BobPQNode *FarFarNext;   /* 3ieme Noeud */

BobTPri Key;

int ZigZag = FALSE;
int Sens=RIGHT;
int SensTop=RIGHT;
int Insere=FALSE;

    n->Link[LEFT] = NULL;
    n->Link[RIGHT] = NULL;
    Bob_NDLINIT(n);

    Key = n->an.Pri;       /* difficult enq */

    Top= NULL;
    Dummy.Link[LEFT] = NULL;
    Dummy.Link[RIGHT]= NULL;

    Cas[LEFT] = &Dummy;
    Cas[RIGHT]= &Dummy;

    Bob_PQLON(pq);
    Next = pq->Root;
/*-- Cas Top = pq->Root --*/

    if( Next == NULL ) {      /* trivial enq */
        pq->Root = n;
        Bob_PQLOFF(pq);
        return;
    }

    if ( Bob_PRIL(Next->an.Pri,Key) )
       Sens = LEFT;

    while ( Sens!=DONE ) {
          Bob_NDLON(Next);
          FarNext = Next->Link[Sens];
          if ( FarNext == NULL ) {
             Cas[1-Sens]->Link[Sens]=Next;
             Cas[Sens]->Link[1-Sens]=NULL;
             n->Link[Sens]   = Dummy.Link[1-Sens];
             n->Link[1-Sens] = Dummy.Link[Sens];
             
             if ( Top==NULL ) {
                pq->Root = n;
                Bob_PQLOFF(pq);
             }
             else {
                Top->Link[SensTop]=n;
                Bob_NDLOFF(Top);
             }
             Bob_NDLOFF(Next);
             
             Sens = DONE;    /* job done, entire tree split */
             continue;
          }
          if ( ( Sens==RIGHT && Bob_PRIL(FarNext->an.Pri,Key) ) ||
               ( Sens==LEFT  && Bob_PRIGE(FarNext->an.Pri,Key) ) ) {
             ZigZag=TRUE;
             Sens=1-Sens;
          }   
          Bob_NDLON(FarNext);
          FarFarNext = FarNext->Link[Sens];
          if ( FarFarNext == NULL ) {
             FarNext->Link[Sens] = n;
             FarFarNext = n;
             Insere = TRUE;
          }
          if ( ZigZag==TRUE ) {
             /* CAS ZIG-ZAG */
             ZigZag=FALSE;
             Cas[1-Sens]->Link[Sens] = FarNext;
             Cas[1-Sens] = FarNext;
             Bob_NDLOFF(FarNext);

             Cas[Sens]->Link[1-Sens] = Next;
             Cas[Sens] = Next;
             Bob_NDLOFF(Next);

             Next = FarFarNext;
             if ( Insere == TRUE ) {
                Bob_NDLON(Next);
                if ( Top == NULL ) {
                   pq->Root = Next;
                   Bob_PQLOFF(pq);
                } else {
                   Top->Link[SensTop]=Next;
                   Bob_NDLOFF(Top);
                }
                Next->Link[Sens] = Dummy.Link[1-Sens];
                Next->Link[1-Sens] = Dummy.Link[Sens];
                Cas[Sens]->Link[1-Sens] = NULL;
                Cas[1-Sens]->Link[Sens] = NULL;
                Bob_NDLOFF(Next);
                Sens = DONE;
                continue;
             }
          }
          else {
             /* CAS ZIG-ZIG */   
             if ( Top==NULL ) {
                pq->Root = FarNext;
                Bob_PQLOFF(pq);
             }
             else {
                Top->Link[SensTop]=FarNext;
                Bob_NDLOFF(Top);
             }


             Next->Link[Sens]= FarNext->Link[1-Sens];
             Cas[1-Sens]->Link[Sens] = Next;
             Bob_NDLOFF(Next);

             FarNext->Link[1-Sens] = Dummy.Link[Sens];
             FarNext->Link[Sens] = Dummy.Link[1-Sens];

             if ( Cas[Sens]==&Dummy ) {
                SensTop = Sens;
                Top = FarNext;
             }
             else {
                SensTop = 1-Sens;
                Top = Cas[Sens];
                Bob_NDLON(Top);
                Bob_NDLOFF(FarNext);
             }

             Dummy.Link[Sens] = Dummy.Link[1-Sens] = NULL;
             Cas[Sens] = Cas[1-Sens] = &Dummy;
          
             Next = FarFarNext;
             if ( Insere == TRUE ) {
                Top->Link[SensTop] = Next;
                Bob_NDLOFF(Top);
                Sens = DONE;
                continue;
             }
          }
          if ( ( Sens==RIGHT && Bob_PRIL(Next->an.Pri,Key) ) ||
               ( Sens==LEFT  && Bob_PRIGE(Next->an.Pri,Key) ) ) {
             Sens = 1-Sens;
          }
   } 
} /* Ranger */




/****** Preemption d'un probleme (un noeud de valeur minimale)  *****/
BobNode *Bob_PQDel(dsi)
int dsi;
{
BobPQ *pq=pqs+dsi;
register BobPQNode * p;
register BobPQNode * deq;         /* one to return */
register BobPQNode * Top;
register BobPQNode * Left;        /* the Left child of Top */
register BobPQNode * FarLeft;     /* the FarLeft child of Left */
register BobPQNode * FarFarLeft;  /* the FarFarLeft child of FarLeft */
int Val;

  Bob_PQLON(pq);
  Left = pq->Root;
  if ( Left == NULL ) {
        Bob_PQLOFF(pq);
        return NULL;
  }

  Top = NULL;
  for (;;) {
       /* Left is not it, FarLeft is not NULL, might be it */
       Bob_NDLON(Left);
       FarLeft = Left->Link[LEFT];
       if ( FarLeft == NULL ) {
          deq = Left;
          if ( Top==NULL ) {
             pq->Root = Left->Link[RIGHT];
             Bob_PQLOFF(pq);
          }
          else {
             Top->Link[LEFT] = Left->Link[RIGHT];
             Bob_NDLOFF(Top);
          }
          Bob_NDLOFF(Left);
          break;
       }
       Bob_NDLON(FarLeft);
       FarFarLeft = FarLeft->Link[LEFT];
       if ( FarFarLeft == NULL ) {
          deq = FarLeft;
          Left->Link[LEFT] = FarLeft->Link[RIGHT];
          if ( Top==NULL ) {
             Bob_PQLOFF(pq);
          } else {
             Bob_NDLOFF(Top);
          }
          Bob_NDLOFF(Left);
          Bob_NDLOFF(FarLeft);
          break;
       }
       /* Left, FarLeft, farFarLeft are not it, rotate */
       if ( Top==NULL ) {
          pq->Root = FarLeft;
          Bob_PQLOFF(pq);
       }
       else {
          Top->Link[LEFT] = FarLeft;
          Bob_NDLOFF(Top);
       }

       Left->Link[LEFT] = FarLeft->Link[RIGHT];
       FarLeft->Link[RIGHT] = Left;
       Bob_NDLOFF(Left);
       Top = FarLeft;
       Left = FarFarLeft;
  }
  return(&(deq->an));

}



/*--------------------------------------------------------------*/
int RecDel(sr,Val)
BobPQNode *sr;
int Val;
{
int nb;

   if ( sr==NULL ) return 0;
   nb = RecDel(sr->Link[LEFT],Val);
   nb += RecDel(sr->Link[RIGHT],Val);
   nb ++;
   Bob_NodeFree(sr);
   return nb;
}


int Bob_PQDelG(dsi,Val)
int dsi;
int Val;
{
BobPQ *pq=pqs+dsi;
BobPQNode *Top;          /* Noeud de raccrochement */
BobPQNode Dummy;         /* Reference de Left et Right */
BobPQNode *Cas[2];       /* Left et Right */
BobPQNode *Next;         /* 1er Noeud */

int Key;
int Sens;
int SensTop;
int Nb;

    Key = Val+1;       /* difficult enq */

    Top= NULL;
    Dummy.Link[LEFT] = NULL;
    Dummy.Link[RIGHT]= NULL;

    Cas[LEFT] = &Dummy;
    Cas[RIGHT]= &Dummy;

    Bob_PQLON(pq);
    Next = pq->Root;
/*-- Cas Top = pq->Root --*/


    while ( Next!=NULL ) {
          if ( Bob_EVALG(Bob_PRIEVAL(Next->an.Pri),Key) ) {
             Cas[LEFT]->Link[RIGHT] = Next;
             Cas[LEFT] = Next;
             Bob_NDLON(Next);
             if ( Top== NULL ) {
                pq->Root = Dummy.Link[RIGHT];
                Bob_PQLOFF(pq);
             } else {  
                Bob_NDLOFF(Top);
             }
             Top = Next;
             Next = Next->Link[RIGHT];
             continue;
          }
          if ( Bob_EVALLE(Bob_PRIEVAL(Next->an.Pri),Key) ) {
             Cas[RIGHT]->Link[LEFT] = Next;
             Cas[RIGHT] = Next;
             Bob_NDLON(Next);
             Bob_NDLOFF(Next);
             Next = Next->Link[LEFT];
          }
   }
   if ( Top== NULL ) {
      pq->Root = Dummy.Link[RIGHT];
      Bob_PQLOFF(pq);
   } else {  
      Top->Link[RIGHT] = NULL;
      Bob_NDLOFF(Top);
   }
   if ( Cas[RIGHT]!= &Dummy ) {
      Cas[RIGHT]->Link[LEFT] = NULL;
      Bob_NDLOFF(Cas[RIGHT]);
   }
   
   Nb = RecDel(Dummy.Link[LEFT],Val);
   return Nb;
}


