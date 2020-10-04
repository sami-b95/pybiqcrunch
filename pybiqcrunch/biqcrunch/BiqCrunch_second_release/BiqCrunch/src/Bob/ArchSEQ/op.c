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
 *  File   : op.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Global PQ for ARTP=SEQ.
 */

#include <stdio.h>
#include "../Include/bb.h"

/*-------------------------- TST SEQ --------------------------------*/
Bob_GPQAlloc()
{
int i,j;

  BobCt.NbDS=1;
  Bob_PQAlloc(1,Bob_PRIEVAL(BobRoot->Pri),Bob_ULBGet());
}

Bob_GPQFree() {
  Bob_PQFree();
}

void Bob_GPQIns(p)
BobNode *p;
{
  Bob_PQIns(0,p);
  BobPsSt.Ins++;
  Bob_DsStIns();

}

BobNode *Bob_GPQDel()
{
  BobNode *p;
  p = Bob_PQDel(0);
  if ( p!=NULL ) {
     BobPsSt.Del++;
     Bob_DsStDel(p->Pri);
     BobPsSt.LastPri = p->Pri;
  }
  return p;
}

void Bob_GPQDelG(ub)
int ub;
{
int nb;

  nb = Bob_PQDelG(0,ub);
  BobPsSt.Delg+=nb;
  dst.NbNode-=nb;
}


