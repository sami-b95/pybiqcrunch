# *****************************************************************************
# *                                                                           *
# *  BiqCrunch is a semidefinite-based solver for binary quadratic problems.  *
# *  It uses a branch-and-bound method featuring an improved semidefinite     *
# *  bounding procedure, mixed with a polyhedral approach. BiqCrunch uses     *
# *  particular input files format to describe the combinatorial problems,    *
# *  similar to the SDPA format.                                              *
# *                                                                           *
# *  Copyright (C) 2010-2016 Nathan Krislock, Jérôme Malick, Frédéric Roupin  *
# *                                                                           *
# *****************************************************************************
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# *****************************************************************************

#### Definitions ####

# Directories
#MKLROOT	=  /opt/intel/mkl
BCDIR  = BiqCrunch
BOBDIR  = Bob
LBFGSBDIR = Lbfgsb.3.0

# Compilers
CC = gcc
FC = gfortran

# Bob makefile parameters
ARTP=SEQ
PRTP=EVDP

BOBOBJ   = Obj/$(PRTP)$(ARTP)

# Get list of problems
PROBLEMLIST = $(shell echo `ls ../problems/`)
BINS = $(addprefix ../problems/, $(addsuffix /biqcrunch, $(PROBLEMLIST)))

# Libraries
LIBS = lib/libbiqcrunch.a

# BiqCrunch objects
BCOBJS = $(BOBOBJ)redirect_output.o $(BOBOBJ)bound.o $(BOBOBJ)build_problem.o \
         $(BOBOBJ)read_parameters.o $(BOBOBJ)output_file.o $(BOBOBJ)update_solution.o \
         $(BOBOBJ)sparse_matrix.o $(BOBOBJ)bounding_procedure.o $(BOBOBJ)inequalities.o \
         $(BOBOBJ)alloc_free.o $(BOBOBJ)bob_functions.o $(BOBOBJ)evaluate_node.o \
	 $(BOBOBJ)rounding.o

# LBFGSB objects
LBFGSBOBJS =  $(BOBOBJ)blas.o $(BOBOBJ)lbfgsb.o $(BOBOBJ)timer.o $(BOBOBJ)linpack.o

# Bob objects
BOBOBJS = $(BOBOBJ)op.o $(BOBOBJ)upbd.o $(BOBOBJ)goproc.o Obj/$(ARTP)start.o $(BOBOBJ)skh.o

# All objects
OBJS = $(BCOBJS) $(LBFGSBOBJS) $(BOBOBJS)

# Determine the machine architecture
ifeq ($(shell uname -m),x86_64)
	ARCH = -m64
else
	ARCH = -m32
endif

# Determine if MKL is installed
ifeq ($(MKLROOT),)
	# LAPACK
	LINALG = -llapack -lblas
else
	# MKL
	ifeq ($(ARCH),-m64)
		MKLARCH = 64
	else
		MKLARCH = 
	endif
	LINALG = -lgfortran -I$(MKLROOT)/include -L$(MKLROOT)/lib -lmkl_intel_lp$(MKLARCH) -lmkl_sequential -lmkl_core -pthread
endif

OPTI    = -O3

BCSTDLIB = -Llib -lbiqcrunch -lm $(OPTI) $(LINALG)
BOBCFLG = -DAPPINC=\"bob.h\" -DPRTP=$(PRTP) -DARTP=$(ARTP) -I$(BOBDIR)/Include -I.

CFLAGS = -Wall $(OPTI) $(ARCH) $(BOBCFLG)
FFLAGS = -fno-automatic -fomit-frame-pointer -fPIC $(OPTI) $(ARCH)


#### Rules ####

.PHONY : all clean


# Default rule is to create all binaries #
all : $(BINS)


# Rules for the final binaries #
../problems/%/biqcrunch : $(BOBOBJ)%-heur.o $(LIBS)
	$(FC) -o $@ $< $(BCSTDLIB)


# heur.c rule #
$(BOBOBJ)%-heur.o : ../problems/%/heur.c
	$(CC) $(CFLAGS) -c -o $@ $<


# BOBOBJS rules #
Obj/$(ARTP)start.o : Bob/Common/start.c
	$(CC) -w $(BOBCFLG) -c -o $@ $<

$(BOBOBJ)skh.o : Bob/PQ/Heap/skh.c
	$(CC) -w $(BOBCFLG) -DSTRINC=\"bb.h\" -c -o $@ $<

$(BOBOBJ)%.o : Bob/ArchSEQ/%.c
	$(CC) -w $(BOBCFLG) -c -o $@ $<


# BiqCrunch code rules #
$(BOBOBJ)%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(BOBOBJ)%.o : %.f
	$(FC) $(FFLAGS) -c -o $@ $<

lib/libbiqcrunch.a : $(OBJS)
	$(AR) -cr $@ $^


# L-BFGS-B rule #
$(BOBOBJ)%.o : $(LBFGSBDIR)/%.f
	$(FC) $(FFLAGS) -c -o $@ $<


# Clean rule #
clean :
	rm $(BINS) $(OBJS) $(LIBS)
