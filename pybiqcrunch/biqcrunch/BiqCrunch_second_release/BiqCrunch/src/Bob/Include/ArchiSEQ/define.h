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
 *  File   : define.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: sequential machine types and defines.
 */

typedef int BobTid;

/*-------- Message defines --------*/
#define Bob_exit
#define Bob_mytid
#define Bob_spawn
#define Bob_config
#define Bob_pstat
#define Bob_kill

#define Bob_joingroup
#define Bob_barrier
#define Bob_recv
#define Bob_mcast 

#define Bob_initsend
#define Bob_pkbyte
#define Bob_pkcplx
#define Bob_pkdcplx
#define Bob_pkdouble
#define Bob_pkfloat
#define Bob_pkint
#define Bob_pklong
#define Bob_pkshort 
#define Bob_pkstr
#define Bob_send

#define Bob_probe
#define Bob_bufinfo
#define Bob_nrecv
#define Bob_upkbyte
#define Bob_upkcplx
#define Bob_upkdcplx
#define Bob_upkdouble
#define Bob_upkfloat
#define Bob_upkint
#define Bob_upklong
#define Bob_upkshort
#define Bob_upkstr


/*-------- Locks defines --------*/
typedef char  ARCHLock;
#define ARCH_LCKON(n)
#define ARCH_LCKOFF(n)
#define ARCH_LCKINIT(n)
#define ARCH_LCKDSTR(n)

typedef ARCHLock BobPQLock;
#define Bob_PQLINIT(p)
#define Bob_PQLDSTR(p)
#define Bob_PQLON(p)
#define Bob_PQLOFF(p)
#define Bob_EXLON(p)
#define Bob_EXLOFF(p)

typedef char BobNDLock;
#define Bob_NDLINIT(p)
#define Bob_NDLDSTR(p)
#define Bob_NDLON(p)
#define Bob_NDLOFF(p)


#define MAXPROCS 1

/*--- memory functions ---*/

#define INIT_SHMEM

#define SH_MALLOC malloc

#define SH_FREE   free


/*--- Shared variables defines ---*/

#define VAR_ALIGN

#define VAR_SHARED 

#define VAR_PRIVATE 

/*--- Get ID ---*/
#define ARCH_GETID(notask) 

