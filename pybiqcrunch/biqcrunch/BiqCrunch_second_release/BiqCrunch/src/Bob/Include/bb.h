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
 *  File   : bb.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: header file of BOB Type definitions.
 */

/*-------------------------------------------------------------------------*/
/* Include architecture dependent definitions 
 * (Locks, massage interface,etc)
 */
#include "Define.h"

/*-------------------------------------------------------------------------*/
/* Priority definitions
 */

typedef int BobTEval;

#if PRTP==POT

typedef struct {
    int Dist,PtWk,NbAnc;
} BobTPri;

#elif PRTP==EVDP

typedef struct {
    BobTEval Eval;
    int Depth;
} BobTPri;

#elif PRTP==EVAL

typedef BobTEval BobTPri;

#endif


/*-------------------------------------------------------------------------*/
typedef struct {
   int Size;
   int Off;
} BobTNdInfo;

#define Bob_NDSIZE(n) n->BobNdInfo.Size
#define Bob_NDOFF(n) n->BobNdInfo.Off


/*     Application Include to find the BobNode and BobSolution definitions.
 *     define APPINC with the -D compilation option.
 */
#ifdef APPINC
#include APPINC
#endif

#define NOMORE -1

/*-------------------------------------------------------------------------*/
#include "priority.h"

/*-------------------------------------------------------------------------*/

#if ARTP==SHARED
typedef struct {
    ARCHLock Lock;
    VAR_ALIGN int Val;
} BobTData;
#endif

/*-------------------------------------------------------------------------*/
#if ARTP==SHARED
typedef BobTData BobTULB;
#else
typedef int BobTULB;
#endif

/*-------------------------------------------------------------------------*/
typedef struct {
#if ARTP==SHARED
    ARCHLock Lock;
#endif
    VAR_ALIGN int NbNode;
    int MaxNode;
    BobTPri LastPri;
} BobTDsSt;

#if ARTP==SEQ
extern BobTDsSt dst;
#elif ARTP==SHARED
extern BobTDsSt dst[MAXPROCS];
#elif ARTP==DISTRIB || ARTP==MLTHR
extern BobTDsSt *dstTab;
extern BobTDsSt dst;
#endif

typedef struct {
    VAR_ALIGN int Evl;
    int NoExp;
    int Ins;
    int Del;
    int Delg;
    BobTPri LastPri;
} BobTPsSt;

#if ARTP==SEQ
extern BobTPsSt BobPsSt;
#elif ARTP==SHARED
extern BobTPsSt BobPsSt[MAXPROCS];
#elif ARTP==DISTRIB || ARTP==MLTHR
extern BobTPsSt *BobPsStTab;
extern BobTPsSt BobPsSt;
#endif


/*-------------------------------------------------------------------------*/
typedef struct {
  int Type;
  int Max;
} BobTExpCt;


/*-------------------------------------------------------------------------*/
typedef struct {
    VAR_ALIGN int NbProcs;
    int PbSize;
    BobTExpCt ExpCt;
    int NbDS;
    int Stat;
} BobTCt;

#if ARTP==SHARED
extern VAR_SHARED  BobTCt BobCtTab[MAXPROCS];
#define BobCt BobCtTab[BobId]
#else
extern BobTCt BobCt;
#endif

#define BobPbSize BobCt.PbSize
#define BobSTAT1  BobCt.Stat&1
#define BobSTAT2  BobCt.Stat&2
#define BobSTAT4  BobCt.Stat&4
#define BobSTAT8  BobCt.Stat&8
#define BobSTAT16 BobCt.Stat&16


/*-------------------------------------------------------------------------*/

#include <sys/time.h>

typedef struct {
  struct timeval Init,Comp,End;
  struct timezone TimeZone;
  float SecInit,SecComp;
} BobTTime;


/*-------------------------------------------------------------------------*/
#include "datastruct.h"

/*-------------------------------------------------------------------------*/
BobNode *Bob_GPQDel();
void     Bob_GPQIns();
void     Bob_GPQDelG();

/*-------------------------------------------------------------------------*/

extern VAR_PRIVATE BobTid BobId;
extern VAR_SHARED  BobNode   *BobRoot;
extern VAR_SHARED  BobSolution *BobSol;


extern void Bob_Algo();
extern void Bob_StPrint();

#if ARTP==DISTRIB || ARTP==MLTHR
extern int *tids;
#endif




