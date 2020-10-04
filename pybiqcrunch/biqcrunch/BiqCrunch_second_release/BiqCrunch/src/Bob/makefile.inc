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
#   File   : makefile.inc
#   Author : B. Le Cun.
#   Date   : May 19 1995.
#   Comment: makefile to include in all other makefile.
#            
# 
# 

BOBDIR  = Bob
BOBOBJ   = Obj/$(PRTP)$(ARTP)
BOBINC = $(BOBDIR)/Include/macro.h $(BOBDIR)/Include/bb.h  \
         $(BOBDIR)/Include/datastruct.h $(BOBDIR)/Include/Define.h

ARCHCFLG    = -I$(BOBDIR)/Include

SEQCFLG     = 
SHAREDCFLG  = 
DISTRIBCFLG = -I$(HOME)/pvm3/include
MLTHRCFLG   = -I$(HOME)/marcel/include  -I$(HOME)/pm2/include \
              -I$(HOME)/pvm3/include 
#MLTHRCFLG   = -I$(ATHAPASCAN)/include  -I$(ATHAPASCAN)/include/Ath0a \
#              -I$(HOME)/pvm3/include -D$(ATH_ARCH)

MAKEAPPFLG = $(ARTP) ARTP=$(ARTP) PQTP=$(PQTP) PRTP=$(PRTP) DEBUG=$(DEBUG)
MAKEARCHFLG = $(MDTP) MDTP=$(MDTP) DEBUG=$(DEBUG)

BOBCFLG = -DPRTP=$(PRTP) -DARTP=$(ARTP) $($(ARTP)CFLG) $(ARCHCFLG)
