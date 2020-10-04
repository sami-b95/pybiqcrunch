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
 *  File   : SiSmSp.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Source file of Single-Simple-Semi-Splay Tree.
 */

#define SPLAYSINGLE 

#include <stdio.h>
#include <math.h>
#include STRINC

VAR_SHARED BobPQ *pqs;

void Bob_PQAlloc(NbDS)
int NbDS;
{
int i;

  pqs = (BobPQ *)SH_MALLOC(sizeof(BobPQ)*NbDS);
  if ( pqs==NULL ) return;
  for (i=0;i<NbDS;i++) {
     pqs[i].Root = NULL;
     Bob_PQLINIT(pqs+i);
  }
}

void Bob_PQFree()
{
  SH_FREE(pqs);
}

/*----------------------------------------------------------------*/
Bob_PQCapa()
{
  return 0;
}

Bob_PQInitLocal(){}

/*------ Simple Enqueue Dequeue -------*/
Enq(in,f)
INTPQNode *in;
BobPQNode *f;
{

    f->Next = NULL;
    in->nb++;

    if ( in->Tail==NULL ) {
       in->Head = f;
    } else {
       in->Tail->Next = f;
    }
    in->Tail = f;
}

int Deq(n,deq)
INTPQNode *n;
BobPQNode **deq;
{

    n->nb--;
    *(deq) = n->Head;
    n->Head = n->Head->Next;
    if ( n->Head == NULL ) {
       n->Tail=NULL;
       return TRUE;
    }
    /*printf("\n");*/
    return FALSE;
}
    

/*-- Simple SemiAdjusting Splay ---*/
void Bob_PQIns(dsi,f)
int dsi;
BobPQNode *f;
{
BobPQ *pq=pqs+dsi;
INTPQNode *n;            /* Noeud a inserer */
INTPQNode *Tmp;          /* Pointeur sur le Noeud a Modifier */
INTPQNode *Top;          /* Noeusd de raccrochement */
INTPQNode Dummy;         /* Reference de Left et Right */
INTPQNode *Cas[2];       /* Left et Right */
INTPQNode *Next;         /* 1er Noeud */
INTPQNode *FarNext;      /* 2ieme Noeud */
INTPQNode *FarFarNext;   /* 3ieme Noeud */

BobTPri Key;
int Alloc=FALSE;
int Sens=RIGHT;
int SensTop=RIGHT;
int Insere=FALSE;


    n = (INTPQNode *)SH_MALLOC(sizeof(INTPQNode));
    n->Key = f->an.Pri;
    n->nb = 0;
    n->Link[LEFT] = NULL;
    n->Link[RIGHT] = NULL;
    n->Head = n->Tail = NULL;
    Bob_NDLINIT(n);

    Key = f->an.Pri;       /* difficult enq */

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
        Enq(n,f);
        Bob_PQLOFF(pq);
        return;
    }

    if ( Bob_PRIE(Next->Key,Key) ) {
        Enq(Next,f);
        Bob_PQLOFF(pq);
        Bob_NDLDSTR(n);
        SH_FREE(n);
        return;
    }
    if ( Bob_PRIL(Next->Key,Key) )
       Sens = LEFT;

    while ( Sens!=DONE ) {
          Bob_NDLON(Next);
          FarNext = Next->Link[Sens];
          if ( FarNext == NULL || Bob_PRIE(FarNext->Key,Key) ) {
             Cas[1-Sens]->Link[Sens]=Next;
             if ( FarNext==NULL ) {
                Next->Link[Sens]=NULL;
                Tmp=n;
                Alloc = TRUE;
                Cas[Sens]->Link[1-Sens]=NULL;
             }
             else {
                Bob_NDLON(FarNext);
                Next->Link[Sens] = FarNext->Link[1-Sens];
                Cas[Sens]->Link[1-Sens] = FarNext->Link[Sens];
                Bob_NDLOFF(FarNext);
                Tmp = FarNext;
             }
             Tmp->Link[Sens]   = Dummy.Link[1-Sens];
             Tmp->Link[1-Sens] = Dummy.Link[Sens];
             Enq(Tmp,f);
             if ( Top==NULL ) {
                pq->Root = Tmp;
                Bob_PQLOFF(pq);
             }
             else {
                Top->Link[SensTop]=Tmp;
                Bob_NDLOFF(Top);
             }
             Bob_NDLOFF(Next);
             
             Sens = DONE;    /* job done, entire tree split */
             continue;
          }

          if ( ( Sens==RIGHT && Bob_PRIL(FarNext->Key,Key) ) ||
               ( Sens==LEFT  && Bob_PRIG(FarNext->Key,Key) ) ) {
             Bob_NDLOFF(Next);
             Cas[1-Sens]->Link[Sens]=Next;
             Cas[1-Sens] = Next;
             Next = FarNext;
             Sens=1-Sens;
             continue;
          }   
          Bob_NDLON(FarNext);
          FarFarNext = FarNext->Link[Sens];
          if ( FarFarNext == NULL || Bob_PRIE(FarFarNext->Key,Key) ) {
             Tmp = ( FarFarNext==NULL ? Alloc=TRUE,n :FarFarNext );
             FarNext->Link[Sens] = Tmp;
             FarFarNext = Tmp;
             Enq(Tmp,f);
             Insere = TRUE;
          }
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
          if ( ( Sens==RIGHT && Bob_PRIL(Next->Key,Key) ) ||
               ( Sens==LEFT  && Bob_PRIGE(Next->Key,Key) ) ) {
             Sens = 1-Sens;
          }
   } 
   if ( Alloc==FALSE ) {
      Bob_NDLDSTR(n);
      SH_FREE(n);
   }
} /* Ranger */



/****** Preemption d'un probleme (un noeud de valeur minimale)  *****/
BobNode *Bob_PQDel(dsi)
int dsi;
{
BobPQ *pq=pqs+dsi;
BobPQNode *deq;         /* one to return */
register INTPQNode * Top;
register INTPQNode * Left;        /* the Left child of Top */
register INTPQNode * FarLeft;     /* the FarLeft child of Left */
register INTPQNode * FarFarLeft;  /* the FarFarLeft child of FarLeft */
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
          if ( Deq(Left,&deq) == TRUE ) {
             if ( Top==NULL ) {
                pq->Root = Left->Link[RIGHT];
                Bob_PQLOFF(pq);
             } else {
                Top->Link[LEFT] = Left->Link[RIGHT];
                Bob_NDLOFF(Top);
             }
             Bob_NDLOFF(Left);
             Bob_NDLDSTR(Left);
             SH_FREE(Left);
             break;
          }
          if ( Top==NULL ) { Bob_PQLOFF(pq); }
          else { Bob_NDLOFF(Top); }

          Bob_NDLOFF(Left);
          break;
       }
       Bob_NDLON(FarLeft);
       FarFarLeft = FarLeft->Link[LEFT];
       if ( FarFarLeft == NULL ) {
          if ( Deq(FarLeft,&deq) == TRUE ) {
             Left->Link[LEFT] = FarLeft->Link[RIGHT];
             Bob_NDLOFF(FarLeft);
             Bob_NDLDSTR(FarLeft);
             SH_FREE(FarLeft);
          } else {
             Bob_NDLOFF(FarLeft);
          }
          if ( Top==NULL ) {
             Bob_PQLOFF(pq);
          } else {
             Bob_NDLOFF(Top);
          }
          Bob_NDLOFF(Left);
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
INTPQNode *sr;
int Val;
{
int nb;
BobPQNode *n;

   if ( sr==NULL ) return 0;
   nb = RecDel(sr->Link[LEFT],Val);
   nb += RecDel(sr->Link[RIGHT],Val);
   nb += sr->nb;
   while( sr->nb>0 ) {
     Deq(sr,&n);
     Bob_NodeFree(n);
   }
   Bob_NDLDSTR(sr);
   SH_FREE(sr);
   return nb;
}


int Bob_PQDelG(dsi,Val)
int dsi;
int Val;
{
BobPQ *pq=pqs+dsi;
INTPQNode *Top;          /* Noeud de raccrochement */
INTPQNode Dummy;         /* Reference de Left et Right */
INTPQNode *Cas[2];       /* Left et Right */
INTPQNode *Next;         /* 1er Noeud */

BobTPri Key;
int Sens;
int SensTop;
int Nb;

    Bob_PRIEVAL(Key) = Val+1;       /* difficult enq */

    Top= NULL;
    Dummy.Link[LEFT] = NULL;
    Dummy.Link[RIGHT]= NULL;

    Cas[LEFT] = &Dummy;
    Cas[RIGHT]= &Dummy;

    Bob_PQLON(pq);
    Next = pq->Root;
/*-- Cas Top = pq->Root --*/


    while ( Next!=NULL ) {
          if ( Bob_EVALG(Bob_PRIEVAL(Next->Key),Bob_PRIEVAL(Key)) ) {
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
          if ( Bob_EVALLE(Bob_PRIEVAL(Next->Key),Bob_PRIEVAL(Key)) ) {
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
   
   Nb= RecDel(Dummy.Link[LEFT],Val);
   return Nb;
}



