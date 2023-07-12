include ../../../../../include.mk

RTE_RRTMG_OBJS = \
  rte_rrtmgp_rams.o \
  mo_simple_netcdf.o \
  rams_rrtmgp_driver.o \
  mo_cloud_optics.o \
  mo_load_coefficients.o \
  mo_load_cloud_coefficients.o \
  mo_heating_rates.o
#
# Location of RTE+RRTMGP libraries, module files.
#
RRTMGP_ROOT = ..
COMPILE_DIR = $(PWD)
RRTMGP_BUILD = $(RRTMGP_ROOT)/build
#
# RRTMGP module files
#
FCINCLUDE += -I$(RRTMGP_BUILD)
FCINCLUDE += -I$(COMPILE_DIR)

#
# netcdf library, module files
# Environment variables NCHOME and NFHOME point to root of C and Fortran interfaces respectively -
#
FCINCLUDE += -I$(NETCDF_FORTRAN_ROOT)/include
LDFLAGS   += -L$(NETCDF_FORTRAN_ROOT)/lib -L$(NETCDF_C_ROOT)/lib
LIBS      += -lnetcdff -lnetcdf

VPATH = $(RRTMGP_ROOT)/extensions/cloud_optics:$(RRTMGP_ROOT)/extensions

# Compilation rules
%.o: %.F90
	$(F_COMP) $(F_OPTS3) $(FCINCLUDE) -c $<
%: %.o
	$(F_COMP) $(F_OPTS3) -o $@ $^ $(LDFLAGS) $(LIBS)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_simple_netcdf.o mo_cloud_optics.o mo_load_coefficients.o mo_load_cloud_coefficients.o
ADDITIONS += mo_heating_rates.o rte_rrtmgp_rams.o $(COMPILE_DIR)/micro_prm.mod

#
# Targets
#
all: rams_rrtmgp_driver.o
	$(ARCH) $(COMPILE_DIR)/rams-$(RAMS_VERSION).a $(RTE_RRTMG_OBJS) 
	cp rte_rrtmgp_rams.mod $(COMPILE_DIR)/.
	rm *.o

rams_rrtmgp_driver.o: $(ADDITIONS) rams_rrtmgp_driver.F90

mo_heating_rates.o: mo_heating_rates.F90
mo_cloud_optics.o: mo_cloud_optics.F90
mo_load_coefficients.o: mo_simple_netcdf.o mo_load_coefficients.F90
mo_load_cloud_coefficients.o: mo_simple_netcdf.o mo_cloud_optics.o mo_load_cloud_coefficients.F90
rte_rrtmgp_rams.o: $(COMPILE_DIR)/micro_prm.mod mo_cloud_optics.mod rte_rrtmgp_rams.F90

clean:
	-rm *.o *.mod
