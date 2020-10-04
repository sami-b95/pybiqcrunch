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
#   File   : makefile.com
#   Author : B. Le Cun.
#   Date   : May 19 1995.
#   Comment: makefile for common file compiling.
# 
# 
include Bob/makefile.inc

PQDIR = $(BOBDIR)/PQ
OBJ   = Obj/$(ARTP)

COMOBJ = $(OBJ)start.o 
COMSRC = $(BOBDIR)/Common/

OBJFILES = $(OBJAPP) $(COMOBJ)

FLG = -DAPPINC=\"$(INC)\" -DPLTP=$(PLTP) $(DEBUG) $(BOBCFLG) -DMDTP=$(MDTP)

all : $(COMOBJ)
	$(MAKE) -f $(PQDIR)/makefile.pq $(PQTP) OBJAPP="$(OBJFILES)"
	

$(OBJ)start.o : $(COMSRC)start.c $(LIBINC)
	$(CC) -c $(COMSRC)start.c $(FLG)
	mv start.o $(OBJ)start.o

