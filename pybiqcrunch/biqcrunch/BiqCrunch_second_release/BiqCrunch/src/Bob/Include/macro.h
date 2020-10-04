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
 *  File   : macro.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: header file of configuration macros definitions.
 */


/*---- ORTP Definition -----*/
#define MINIMISATION 0
#define MAXIMISATION 1

/*---- ARTP Definitions ----*/
/*--- Architecture ----*/
#define SEQ      1
#define SHARED   2
#define DISTRIB  3
#define MLTHR    4

/*---- MDTP Definitions ----*/
/*---- parallelization methods MDTP for ARTP==DISTRIB||SHARED----*/
#define NONE 0
#define ASC 1
#define SYC 2
#define FRM 3

/*---- parallelization methods MDTP for ARTP==MLTHR ----*/
/*#define NONE 0*/
#define PM2 1
#define ATHAP 2

/*---- DSTP Definitions ----*/
/*---- Data structure type for shared memory architecture ----*/
/*#define NONE 0*/
#define ONE 1
#define FXD 2
#define MRG 3

/*---- PLTP Definitions ----*/
/*---- Partial locking protocol ----*/
/* the definition of SEQ is in ARTP section */
/*#define NONE 0*/
/*#define SEQ 1*/
#define TST 2
#define EXC 3
#define MRK 4
#define LCK 5

/*--- PRTP Definitions ----*/
/*---- Priority notions ----*/

/*#define NONE 0*/
#define POT  1
#define EVDP 2
#define EVAL 3



