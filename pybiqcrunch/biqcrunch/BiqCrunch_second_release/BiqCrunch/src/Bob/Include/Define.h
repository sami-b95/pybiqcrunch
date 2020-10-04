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
 *  File   : Define.h
 *  Author : B. Le Cun.
 *  Date   : May 19 1995.
 *  Comment: header file of architecture dependent definitions.
 */

#include "macro.h"

#define NEEDMAIN

#if ARTP==SHARED
#include "ArchiSHARED/define.h"
#elif ARTP==DISTRIB
#include "ArchiDISTRIB/define.h"
#elif ARTP==MLTHR
#if MDTP==ATHAP
#undef NEEDMAIN
#endif
#include "ArchiDISTRIB/define.h"
#else
#include "ArchiSEQ/define.h"
#endif

#ifdef RS6000
#define GETOPTFAIL 255
#else
#define GETOPTFAIL -1
#endif
