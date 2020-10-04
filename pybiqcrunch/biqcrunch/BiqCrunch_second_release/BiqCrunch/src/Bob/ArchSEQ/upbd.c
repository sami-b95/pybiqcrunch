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
 *  File   : upbd.c
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: Global Upper/Lower Bound for ARTP=SEQ
 */

#include <stdio.h>
#include "../Include/bb.h"

BobTULB BobUpBd;

/*----------------------------------------------------------------------*/
int Bob_ULBGet() {
  return BobUpBd;
}

/*----------------------------------------------------------------------*/
void Bob_ULBInit(ub,bs)
int ub;
BobSolution *bs;
{
int i;

  BobSol = (BobSolution *)malloc(sizeof(BobSolution));
  if ( BobSol == NULL ) {
     fprintf(stderr,"Not enough memory for BobSolution\n");
     exit(1);
  }
  *BobSol = *bs;
  BobUpBd= ub;
}

/*----------------------------------------------------------------------*/
int Bob_ULBUpd(NewBks,bs) 
int NewBks;
BobSolution *bs;
{

  if ( Bob_EVALL(BobUpBd,NewBks) ) {
     BobUpBd = NewBks;
     *BobSol = *bs;
     Bob_GPQDelG(NewBks);
     return 1;
  } 
  return 0;
}

/*----------------------------------------------------------------------*/
int Bob_ULBSup(Pri) 
int Pri;
{

  return (Bob_EVALL(BobUpBd,Pri));
}

