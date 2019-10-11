# Makefile for the 3D KPP boundary-layer ocean model
# Written by Nicholas P. Klingaman
# University of Reading, Department of Meteorology and NCAS-Climate
# Substantially revised on 03 Sept. 2009 to remove separate dot files for each target

# NOTE: You *must* run <make clean> if you wish to recompile the model for another target
# (e.g., if you have compiled for OASIS2 and wish to recompile for OASIS3, you must run <make clean> first,
# then run <make oasis3_coupled>.)

# Makefile for HECToR - Linux with Pathscale Fortran compiler
# NOTE: The CFS/GFS targets do not work on HECToR as the CFS/GFS libraries are not installed here.

# Fortran compiler
F90=gfortran-mp-4.8
#F90=ifort

# Flags for the Fortran compiler
#F90_FLAGS=-x -fno-second-underscore -I. -pedantic
F90_FLAGS=-x f77-cpp-input -O3 -fdefault-real-8 
#F90_FLAGS=-pedantic -fno-second-underscore 
#F90_FLAGS=-fpp -shared-intel -mcmodel=large -O3 -fast -r8 -fp-model precise
#F90_FLAGS=-fpp -O3 -r8 -traceback -fp-model precise -mcmodel=large
#F90_FLAGS=-fpp -O3 -r8 -traceback -mcmodel=large

# NetCDF dynamic libraries (not required on HECToR; use `module load netcdf` instead)
NCDF_LIB=-L/opt/local/lib -lnetcdf -lnetcdff
#NCDF_LIB=-lnetcdf -lnetcdff
# NetCDF include files (not required on HECToR ; use `module load netcdf` instead)
NCDF_INC=-I/opt/local/include

# Name of the KPP executable
EXECUTABLE=KPP_ocean
# -- No user modifications beyond this point --

FORCED_FILES = byjas_1D_ocn.f ocn.f vmix.f kpp.f subs1D.f steves_fluxes.f ncdf_in_out.f  

FORCED_OBJS = $(FORCED_FILES:.f=.o)

FORCED_LINK = $(F90) $(FORCED_OBJS) -o ${EXECUTABLE} ${NCDF_INC} ${NCDF_LIB}

.SUFFIXES:
.SUFFIXES: .f .f90 .o

.f.o :
	$(F90) $(F90_FLAGS) $(NCDF_INC) -c $< -o $@
.f90.0 :
	$(F90) $(F90_FLAGS) $(NCDF_INC) -c $< -o $@

forced : $(FORCED_OBJS)
	$(FORCED_LINK)

clean : 
	rm -rf *.i *.o *.mod $(EXECUTABLE)

obj_clean : 
	rm -rf *.i *.o *.mod

%:%.o
	$(F90) $(LDFLAGS) -o $@ $^ $(LDLIBS)
