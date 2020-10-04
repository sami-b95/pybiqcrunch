# 
# 
#         BOB version 1.0:  Branch and Bound Optimization LiBrary
#                     PNN Team of PRiSM laboratory
#              University of Versailles St-Quentin en Yvelines.
#       Authors:  M. Benaichouche, V. Cung, S. Dowaji, B. Le Cun
#                       T. Mautor, C. Roucairol.
#                     (C) 1995 All Rights Reserved
# 
#                               NOTICE
# 
#  Permission to use, copy, modify, and distribute this software and
#  its documentation for any purpose and without fee is hereby granted
#  provided that the above copyright notice appear in all copies and
#  that both the copyright notice and this permission notice appear in
#  supporting documentation.
# 
#  Neither the institutions (Versailles University, PRiSM Laboratory, 
#  the PNN Team), nor the Authors make any representations about the 
#  suitability of this software for any purpose.  This software is 
#  provided ``as is'' without express or implied warranty.
# 
# 
# 
#  
#   File   : makefile.pq
#   Author : B. Le Cun.
#   Date   : May 19 1995.
#   Comment: makefile for PQ compiling and linking.
# 

include Bob/makefile.inc

PQDIR = $(BOBDIR)/PQ
SPLAY = $(PQDIR)/Splay/
FUN   = $(PQDIR)/Funnel/
HEAP  = $(PQDIR)/Heap/
BIN   = Bin/$(PRTP)$(MDTP)$(DSTP)$(LBTP)$(PLTP)
OBJFILES = $(OBJAPP) 

SSPLAYBIN = $(BIN)SiSmSp $(BIN)SiSp 
NSPLAYBIN = $(BIN)SmSp $(BIN)Sp 
SKEWBIN  = $(BIN)skh 
FUNBIN   = $(BIN)futb $(BIN)futr 
HEAPBIN  = $(BIN)dh
OTHERSBIN= $(BIN)lh $(BIN)ph

INCAPP   = -DAPPINC=\"$(INC)\" -DSTRINC=\"../../Include/bb.h\"

PARALLEL = $(LIB) 

DTP = $(BOBCFLG) -DPLTP=$(PLTP) -DDSTP=$(DSTP) -DLBTP=$(LBTP)

DSOPT = $(INCAPP) $(DEBUG) $(PARALLEL) $(DTP)

ALL : $(PLTP) $(PQTP)


#LCK : SPLAY SKEW FUNNEL HEAP
LCK : SKEW HEAP FUNNEL SPLAY

TST : SKEW HEAP FUNNEL SPLAY OTHERS

SEQ : SPLAY SKEW FUNNEL HEAP OTHERS
#SEQ : OTHERS

SPG : SPLAY SKEW FUNNEL HEAP

MRK : SPLAY SKEW FUNNEL 
#MRK : SKEW 


SPLAY : $(NSPLAYBIN) $(SSPLAYBIN)

NSPLAY : $(NSPLAYBIN) 

SSPLAY : $(SSPLAYBIN)

SKEW  : $(SKEWBIN)

FUNNEL : $(FUNBIN)

HEAP  : $(HEAPBIN)

OTHERS : $(OTHERSBIN)

NOPQ :  $(BIN)nopq

SEQCmd =
SHAREDCmd =
DISTRIBCmd = cp $@ $(PVM_ROOT)/bin/$(PVM_ARCH)
MLTHRCmd = cp $@ $(PVM_ROOT)/bin/$(PVM_ARCH)



$(BIN)nopq : $(OBJFILES)
	$(CC) $(HEAP)/nopq.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)SiSmSp  : $(SPLAY)SiSmSp.c $(OBJFILES)
	$(CC) $(SPLAY)SiSmSp.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)
	

$(BIN)SiSp    : $(SPLAY)SiSp.c $(OBJFILES)
	$(CC) $(SPLAY)SiSp.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)SmSp    : $(SPLAY)SmSp.c $(OBJFILES)
	$(CC) $(SPLAY)SmSp.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)Sp      : $(SPLAY)Sp.c $(OBJFILES)
	$(CC) $(SPLAY)Sp.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)nSp      : $(SPLAY)nSp.c $(OBJFILES)
	$(CC) $(SPLAY)nSp.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)skh    : $(HEAP)skh.c $(OBJFILES)
	$(CC) $(HEAP)skh.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)futb    : $(FUN)futb.c $(OBJFILES)
	$(CC) $(FUN)futb.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)futr    : $(FUN)futr.c $(OBJFILES)
	$(CC) $(FUN)futr.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)dh    : $(HEAP)dh.c $(OBJFILES)
	$(CC) $(HEAP)dh.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)lh    : $(HEAP)lh.c $(OBJFILES)
	$(CC) $(HEAP)lh.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

$(BIN)ph    : $(HEAP)ph.c $(OBJFILES)
	$(CC) $(HEAP)ph.c -o $@ $(OBJFILES) $(DSOPT)
	$($(ARTP)Cmd)

