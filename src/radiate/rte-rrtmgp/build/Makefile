#!/usr/bin/env make
include ../../../../../include.mk

COMPILE_DIR = $(PWD)
RTE_DIR = ../rte
RRTMGP_DIR = ../rrtmgp
#
# Choose kernels depending on environment variable
#
ifeq ($(RTE_KERNELS), openacc)
  RTE_KERNEL_DIR    = ../rte/kernels-openacc
  RRTMGP_KERNEL_DIR = ../rrtmgp/kernels-openacc
endif
RTE_KERNEL_DIR += ../rte/kernels
RRTMGP_KERNEL_DIR += ../rrtmgp/kernels
VPATH = $(RTE_DIR):$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR)
#
# Compiler variables FC, FCFLAGS must be set in the environment
#
# Make all the libraries though we'll only use the interface + kernels
all: for_rams
#all:           librte.a                           librrtmgp.a
separate-libs: librtekernels.a librteinterfacef.a librrtmgpkernels.a librrtmgpinterfacef.a

COMPILE = $(F_COMP) $(F_OPTS3) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<

include $(RTE_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends

for_rams: $(RTE_FORTRAN_KERNELS) $(RRTMGP_FORTRAN_KERNELS) $(RRTMGP_FORTRAN_INTERFACE) $(RTE_FORTRAN_INTERFACE)
	$(ARCH) $(COMPILE_DIR)/rams-$(RAMS_VERSION).a $(RTE_FORTRAN_KERNELS) $(RRTMGP_FORTRAN_KERNELS) \
	$(RRTMGP_FORTRAN_INTERFACE) $(RTE_FORTRAN_INTERFACE) 
	rm *.o

clean:
	rm -f *.mod *.o 
