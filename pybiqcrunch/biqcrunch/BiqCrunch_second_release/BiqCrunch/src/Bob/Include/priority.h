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
 *  File   : priority.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: header file of macros used for priority operations.
 */

/*------------------------ Evaluation Define -------------------*/

#define BobEvalINFI 200000
#define BobEvalNULL 0
#ifndef ORTP
#define ORTP MINIMISATION
#endif

#if ORTP==MINIMISATION
#define Bob_EVALE(e1,e2)  (e1==e2)
#define Bob_EVALL(e1,e2)  (e1>e2)
#define Bob_EVALLE(e1,e2) (e1>=e2)
#define Bob_EVALG(e1,e2)  (e1<e2)
#define Bob_EVALGE(e1,e2) (e1<=e2)
#elif ORTP==MAXIMISATION
#define Bob_EVALE(e1,e2)  (e1==e2)
#define Bob_EVALL(e1,e2)  (e1<e2)
#define Bob_EVALLE(e1,e2) (e1<=e2)
#define Bob_EVALG(e1,e2)  (e1>e2)
#define Bob_EVALGE(e1,e2) (e1>=e2)
#endif

/*------------------------ Priority Define -------------------*/
#if PRTP==POT

#define Bob_PRIE(p1,p2)  (p1.Dist ==p2.Dist) && \
              (p1.PtWk ==p2.PtWk) && \
              (p1.NbAnc==p2.NbAnc )
#define Bob_PRIL(p1,p2)  (p1.Dist > p2.Dist) || \
              ((p1.Dist==p2.Dist) && (p1.PtWk<p2.PtWk)) || \
              ((p1.Dist==p2.Dist) && (p1.PtWk==p2.PtWk) && (p1.NbAnc<p2.NbAnc))
#define Bob_PRILE(p1,p2) (p1.Dist>p2.Dist) || \
              ((p1.Dist==p2.Dist) && (p1.PtWk<p2.PtWk)) || \
              ((p1.Dist==p2.Dist) && (p1.PtWk==p2.PtWk) && (p1.NbAnc<=p2.NbAnc))
#define Bob_PRIG(p1,p2)  (p1.Dist<p2.Dist) || \
              ((p1.Dist==p2.Dist) && (p1.PtWk>p2.PtWk)) || \
              ((p1.Dist==p2.Dist) && (p1.PtWk==p2.PtWk) && (p1.NbAnc>p2.NbAnc))
#define Bob_PRIGE(p1,p2) (p1.Dist<p2.Dist) || \
              ((p1.Dist==p2.Dist) && (p1.PtWk>p2.PtWk)) || \
              ((p1.Dist==p2.Dist) && (p1.PtWk==p2.PtWk) && (p1.NbAnc>=p2.NbAnc))

#define Bob_PRICREAT(Pri,eval,Wk,lvl) Pri.Dist  = eval; \
                                         Pri.PtWk = Wk; \
                                         Pri.NbAnc = lvl
#define Bob_PRINULL(Pri) Pri.Dist=0;Pri.PtWk=0;Pri.NbAnc=0
#define Bob_PRIINFI(Pri) Pri.Dist=0;Pri.PtWk=0;Pri.NbAnc=0

#define Bob_PRISTR(io)       fprintf(io,"( Dist   PtWk NbAn)");
#define Bob_PRIPRINT(io,Pri) fprintf(io,"%6d %6d %3d",\
                            Pri.Dist,Pri.PtWk,Pri.NbAnc);

#if ARTP==DISTRIB || ARTP==MLTHR
#define Bob_pkPRI(Pri) Bob_pkint(&(Pri.Dist),1,1);\
                       Bob_pkint(&(Pri.PtWk),1,1);\
                       Bob_pkint(&(Pri.NbAnc),1,1)
#define Bob_upkPRI(Pri) Bob_upkint(&(Pri.Dist),1,1);\
                       Bob_upkint(&(Pri.PtWk),1,1);\
                       Bob_upkint(&(Pri.NbAnc),1,1)
#endif

/*----------------------------------*/
#elif PRTP==EVDP

#define Bob_PRI1E(p1,p2)  (Bob_EVALE(p1.Eval,p2.Eval)) 
#define Bob_PRI2E(p1,p2)  (p1.Depth==p2.Depth) 
#define Bob_PRIE(p1,p2)   Bob_PRI1E(p1,p2) && Bob_PRI2E(p1,p2)

#define Bob_PRI1L(p1,p2)  (Bob_EVALL(p1.Eval,p2.Eval)) 
#define Bob_PRI2L(p1,p2)  (p1.Depth<p2.Depth)
#define Bob_PRIL(p1,p2)   Bob_PRI1L(p1,p2) || (Bob_PRI1E(p1,p2) && Bob_PRI2L(p1,p2))

#define Bob_PRI1LE(p1,p2) (Bob_EVALLE(p1.Eval,p2.Eval))
#define Bob_PRI2LE(p1,p2) (p1.Depth<=p2.Depth)
#define Bob_PRILE(p1,p2)  Bob_PRI1LE(p1,p2) || (Bob_PRI1E(p1,p2) && Bob_PRI2LE(p1,p2))

#define Bob_PRI1G(p1,p2)  (Bob_EVALG(p1.Eval,p2.Eval))
#define Bob_PRI2G(p1,p2)  (p1.Depth>p2.Depth)
#define Bob_PRIG(p1,p2)   Bob_PRI1G(p1,p2) || (Bob_PRI1E(p1,p2) && Bob_PRI2G(p1,p2))

#define Bob_PRI1GE(p1,p2) (Bob_EVALG(p1.Eval,p2.Eval))
#define Bob_PRI2GE(p1,p2) (p1.Depth>=p2.Depth)
#define Bob_PRIGE(p1,p2)  Bob_PRI1G(p1,p2) || (Bob_PRI1E(p1,p2) && Bob_PRI2GE(p1,p2))

#define Bob_PRICREAT(Pri,eval,Wk,lvl) Pri.Eval=eval;Pri.Depth=lvl

#define Bob_PRINULL(Pri) Pri.Eval=BobEvalNULL;Pri.Depth=0
#define Bob_PRIINFI(Pri) Pri.Eval=BobEvalINFI;Pri.Depth=0
#define Bob_PRIEVAL(Pri) Pri.Eval

#define Bob_PRISTR(io)       fprintf(io,"( Eval Depth)");
#define Bob_PRIPRINT(io,Pri) fprintf(io,"%6d %4d",Pri.Eval,Pri.Depth);

#if ARTP==DISTRIB || ARTP==MLTHR
#define Bob_pkPRI(Pri) Bob_pkint(&(Pri.Eval),1,1);\
                       Bob_pkint(&(Pri.Depth),1,1)
#define Bob_upkPRI(Pri) Bob_upkint(&(Pri.Eval),1,1);\
                       Bob_upkint(&(Pri.Depth),1,1)
#endif

/*----------------------------------*/
#elif PRTP==EVAL

#define Bob_PRI1E(p1,p2)  Bob_EVALE(p1,p2)
#define Bob_PRIE(p1,p2)   Bob_EVALE(p1,p2)
#define Bob_PRI1L(p1,p2)  Bob_EVALL(p1,p2)
#define Bob_PRIL(p1,p2)   Bob_EVALL(p1,p2)
#define Bob_PRI1LE(p1,p2) Bob_EVALLE(p1,p2)
#define Bob_PRILE(p1,p2)  Bob_EVALLE(p1,p2)
#define Bob_PRI1G(p1,p2)  Bob_EVALG(p1,p2)
#define Bob_PRIG(p1,p2)   Bob_EVALG(p1,p2)
#define Bob_PRI1GE(p1,p2) Bob_EVALGE(p1,p2)
#define Bob_PRIGE(p1,p2)  Bob_EVALGE(p1,p2)
       
#define Bob_PRICREAT(Pri,eval,PtWk,lvl) Pri = eval

#define Bob_PRINULL(Pri) Pri=BobEvalNULL
#define Bob_PRIINFI(Pri) Pri=BobEvalINFI
#define Bob_PRIEVAL(Pri) Pri

#define Bob_PRISTR(io)       fprintf(io,"( Eval)");
#define Bob_PRIPRINT(io,Pri) fprintf(io,"%6d",Pri);

#if ARTP==DISTRIB || ARTP==MLTHR
#define Bob_pkPRI(Pri) Bob_pkint(&(Pri),1,1)
#define Bob_upkPRI(Pri) Bob_upkint(&(Pri),1,1)
#endif

#endif

