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
 *  File   : datastruct.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: header file of PQ Types definitions.
 *           - Pairing,Leftist,Skew,Rao&Kumar Heap;
 *           - Semi-Splay,Simple-Semi-Splay, 
 *             Single-Semi-Splay,Single-Simple-Semi-Splay Tree. 
 *           - Funnel-Table, Funnel-Tree.
 *           - without PQ.
 */


/*-------------------------------------------------------------------------*/
#ifdef PAIRINGHEAP
#define DATASTRUCT

typedef struct __PQNode {
  VAR_ALIGN BobNode an;
  struct __PQNode *ff,*fs;
} BobPQNode;
          
typedef struct __Heap { 
   BobPQLock Lock; 
   VAR_ALIGN BobPQNode *Root;
} BobPQ;

#endif
/*-------------------------------------------------------------------------*/
#ifdef LEFTISTHEAP
#define DATASTRUCT

typedef struct __PQNode {
   VAR_ALIGN BobNode an;
   int Rank;
   struct __PQNode *fg,*fd;
} BobPQNode;
          
typedef struct __Heap { 
   BobPQLock Lock; 
   VAR_ALIGN BobPQNode *Root;
} BobPQ;

#endif

/*-------------------------------------------------------------------------*/
#ifdef SKEWHEAP
#define DATASTRUCT

typedef struct __PQNode {
   VAR_ALIGN BobNode an;
   BobNDLock Lock; 
   struct __PQNode *left,*right;
} BobPQNode;
          
typedef struct __Heap { 
   VAR_ALIGN BobPQLock Lock; 
   BobPQNode *Root;
} BobPQ;

#endif

/*-------------------------------------------------------------------------*/
#ifdef RKHEAP
#define DATASTRUCT

typedef struct {
   VAR_ALIGN BobNode an;
} BobPQNode;

typedef struct { 
   VAR_ALIGN BobPQNode *n;
   BobNDLock Lock;
   int Status;
} INTPQNode;


typedef struct { 
   BobPQLock Lock;
   int lastelem;
   int fulllevel;
   VAR_ALIGN INTPQNode q[NBELEM]; 
} BobPQ;

#endif

/*-------------------------------------------------------------------------*/
#ifdef SPLAYSINGLE
#define DATASTRUCT

#define DONE -1 
#define LEFT  0 
#define RIGHT 1 
#define NONE  0 

typedef struct __PQNode {
    VAR_ALIGN BobNode an;
    struct __PQNode *Next;
} BobPQNode; 
          
typedef struct __IntNode {                                           
    BobNDLock Lock;
    VAR_ALIGN BobPQNode *Head,*Tail;                                                
    int nb;                                                          
    struct __IntNode *Link[2];                                       
    BobTPri Key;
} INTPQNode;                                                           
          
typedef struct __Heap {                                              
   BobPQLock Lock;                                                  
   VAR_ALIGN INTPQNode *Root;              /* root node */                      
} BobPQ;                                                              

#endif

/*-------------------------------------------------------------------------*/
#ifdef SPLAY
#define DATASTRUCT

#define DONE -1
#define LEFT  0
#define RIGHT 1
#define NONE  0

typedef struct __PQNode {
    VAR_ALIGN BobNode an;
    BobNDLock Lock;
    struct __PQNode *Link[2];
} BobPQNode;
 
 
typedef struct { 
   BobPQLock Lock;
   VAR_ALIGN BobPQNode *Root;         /* root node */
} BobPQ;

#endif

/*-------------------------------------------------------------------------*/
#ifdef FUNNELTABLE
#define DATASTRUCT

typedef struct __PQNode { 
  VAR_ALIGN BobNode an;
  struct __PQNode  *Next;
} BobPQNode;

 
typedef struct {
  BobPQLock Lock;
  VAR_ALIGN BobPQNode *head, *tail;
  int NbProc;
} FIFO;

typedef struct {
  BobPQLock Lock;
  VAR_ALIGN int binf;
  FIFO *funnel;
} BobPQ;

#endif

/*-------------------------------------------------------------------------*/
#ifdef FUNNELTREE
#define DATASTRUCT

typedef struct __PQNode {
  VAR_ALIGN BobNode an;
  struct __PQNode *Next;
} BobPQNode;


typedef struct __Fifo {
  VAR_ALIGN BobPQNode *head, *tail;
} FIFO;


typedef struct __FUNode { 
  BobNDLock Lock;
  VAR_ALIGN int Val;
} BobFUNode;
 
typedef struct {
  BobPQLock Lock;
  VAR_ALIGN int Val;
  BobFUNode *funnel;
  FIFO *queues;
} BobPQ;


#endif

/*-------------------------------------------------------------------------*/

#ifdef NOPQ
#define DATASTRUCT

typedef struct {
  BobNode an;
} BobPQNode;
#endif

#define TRUE 1
#define FALSE 0

/*-------------------------------------------------------------------------*/
#ifdef DATASTRUCT

BobNode *Bob_NodeAlloc(n)
int n;
{
BobPQNode *an;

  an = (BobPQNode *)SH_MALLOC(sizeof(BobPQNode)+n);
  an->an.BobNdInfo.Off = sizeof(BobPQNode);
  an->an.BobNdInfo.Size = n+sizeof(BobPQNode);
  return &(an->an);
}

void Bob_NodeFree(n)
BobPQNode *n;
{
#if defined(SPLAY) || defined(SKEWHEAP)
  Bob_NDLDSTR(n);
#endif
  SH_FREE(n);
}

#else
BobNode *Bob_NodeAlloc();
void     Bob_NodeFree();

void    *Bob_PQAlloc();
void     Bob_PQFree();

BobNode *Bob_PQDel();
void     Bob_PQIns();
int      Bob_PQDelG();

void     Bob_PQSpecial();
void     Bob_PQInitLocal();
#endif

