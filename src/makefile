#-------------------------------------------------------------------------------
#
# include configuration
#
ifeq ($(wildcard make.config),make.config)
$(info )
$(info Using customized "$(wildcard make.config)" file.)
$(info )
include make.config
else
$(info )
$(info Could not find the customized "make.config" file!)
$(info )
$(info Copy "make.default" to "make.config" and customize it, if you do not want the default configuration.)
$(info )
$(shell sleep 5)
include make.default
endif

#-------------------------------------------------------------------------------
#
# host variables
#
HOST    := $(shell hostname -f)
OS      := $(shell uname)
DATE    := $(shell date +%Y%m%d%H%M)

include hosts/default
ifeq ($(wildcard hosts/$(HOST)),hosts/$(HOST))
$(info Using compiler setup from "hosts/$(HOST)" file.)
$(info )
include hosts/$(HOST)
endif

#-------------------------------------------------------------------------------
#
# pass additional flags to compiler
#
ifeq ($(DBLPREC),Y)
FFLAGS  += -r8 -DDBLPREC
endif

ifeq ($(HDF5),Y)
FFLAGS  += -DHDF5
ifeq ($(HDF5COMPRESS),GZIP)
FFLAGS  += -DHDF5GZIP
endif
FFLAGS  += -I${HDF5DIR}/lib
LIBS    += -L${HDF5DIR}/lib -lhdf5_fortran -lhdf5 -lz
else
FFLAGS  += -I${HDF4DIR}/include
LIBS    += -L${HDF4DIR}/lib -lmfhdf -ldf -ljpeg -lz
endif

ifeq ($(MHD),Y)
FFLAGS  += -DMHD
else
VPOT       = N
endif

ifeq ($(VPOT),Y)
FFLAGS  += -DVPOT
endif

ifeq ($(ISO),Y)
FFLAGS  += -DISO
endif

ifeq ($(DIMS),D2)
FFLAGS  += -DR2D
endif

ifeq ($(DIMS),D2.5)
FFLAGS  += -DR2D5
endif

ifeq ($(DIMS),D3)
FFLAGS  += -DR3D
endif

ifeq ($(TIME),RK3)
FFLAGS  += -DRK3
else
ifeq ($(TIME),EULER)
FFLAGS  += -DEULER
else
FFLAGS  += -DRK2
endif
endif

ifeq ($(FLUX),HLL)
FFLAGS  += -DHLL
endif

ifeq ($(FLUX),HLLC)
ifeq ($(MHD),Y)
FFLAGS  += -DHLLD
else
ifeq ($(ISO),Y)
FFLAGS  += -DHLL
else
FFLAGS  += -DHLLC
endif
endif
endif

ifeq ($(FLUX),HLLD)
ifeq ($(MHD),Y)
FFLAGS  += -DHLLD
else
FFLAGS  += -DHLLC
endif
endif

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = godunov-amr

default: $(name).x

sources = blocks.F90 config.F90 driver.F90 error.F90 evolution.F90 io.F90 \
          mesh.F90 problem.F90 timer.F90
objects = blocks.o   config.o   driver.o   error.o   evolution.o   io.o   \
          mesh.o   problem.o   timer.o
files   = $(sources) makefile make.default config.in license.txt hosts

$(name).x: $(objects)
	$(LD) $(LDFLAGS) $(objects) $(LIBS) -o $(name).x

arch: $(files)
	tar czvf $(name)-$(DATE).tar.gz $(files)

clean:
	rm -rf *.x *.o *.mod *.il

clean-all:
	rm -rf *.x *.o *.mod *.il *.out *.dat *.bin *.h5 *.hdf *~ gprof.txt

#-------------------------------------------------------------------------------

blocks.o        : blocks.F90 config.o
config.o        : config.F90 error.o
driver.o        : driver.F90 config.o evolution.o io.o mesh.o timer.o
error.o         : error.F90
evolution.o     : evolution.F90 blocks.o
io.o            : io.F90 blocks.o error.o
mesh.o          : mesh.F90 blocks.o config.o error.o problem.o
problem.o       : problem.F90 blocks.o
timer.o         : timer.F90

#-------------------------------------------------------------------------------
