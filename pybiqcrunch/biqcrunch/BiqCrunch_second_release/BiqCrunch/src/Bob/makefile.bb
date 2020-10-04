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
#   File   : makefile.bb
#   Author : B. Le Cun.
#   Date   : May 19 1995.
#   Comment: makefile for the all library, branching on machine architecture.
# 
# 
LDIR  = Bob

PQFLG = PQTP="$(PQTP)"
LBFLG = LBTP="$(LBTP)"
PRFLG = PRTP=$(PRTP)

FLG = $(LBFLG) $(PQFLG)


SEQLIB     = ARTP=SEQ $(PQFLG) $(PRFLG)
SHAREDLIB  = ARTP=SHARED  PLTP="$(PLTP)" $(FLG) $(PRFLG)
DISTRIBLIB = ARTP=DISTRIB PLTP=SEQ $(FLG) $(PRFLG)
MLTHRLIB   = ARTP=MLTHR PLTP=SEQ $(FLG) $(PRFLG)

DEBUG = -g

ADIR = $(LDIR)/Arch$(ARTP)

ALL : $(HOSTTYPE)
	
# Sequentiel Machine
sun4 : SEQ DISTRIB MLTHR

sun3 : SEQ DISTRIB

# Sequent Balance 
bal : SEQ SHARED DISTRIB

# KSR 1
ksr1 : SEQ SHARED DISTRIB
 

SEQ : 
	rm -f *.o
	$(MAKE) -f $(LDIR)/ArchSEQ/makefile.arch $(PLTP) $(SEQLIB)

SHARED : 
	rm -f *.o
	$(MAKE) -f $(LDIR)/ArchSHARED/makefile.arch $(MDTP) $(SHAREDLIB) 

DISTRIB : 
	rm -f *.o
	$(MAKE) -f $(LDIR)/ArchDISTRIB/makefile.arch $(MDTP) $(DISTRIBLIB) 

MLTHR : 
	rm -f *.o
	$(MAKE) -f $(LDIR)/ArchMLTHR/makefile.arch $(MDTP) $(MLTHRLIB) 


