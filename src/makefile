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
# pass additional flags to the compiler
#
# number of dimensions
#
FFLAGS  += -DNDIMS=${NDIMS}

# output data format
#
FFLAGS  += -D${OUTPUT}

# precision
#
FFLAGS  += -DPREC=${PREC}

# time integration method
#
FFLAGS  += -D${TIME}

# set of equations
#
FFLAGS  += -D${EQS}

# equation of state
#
FFLAGS  += -D${EOS}

# flux approximation method
#
FFLAGS  += -D${FLUX}

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = godunov-amr

default: $(name).x

sources = blocks.F90 boundaries.F90 config.F90 driver.F90 error.F90            \
          evolution.F90 interpolation.F90 io.F90 mesh.F90 problem.F90          \
          scheme.F90 timer.F90
objects = blocks.o   boundaries.o   config.o   driver.o   error.o              \
          evolution.o   interpolation.o   io.o   mesh.o   problem.o            \
          scheme.o timer.o
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
boundaries.o    : boundaries.F90 blocks.o config.o error.o interpolation.o
config.o        : config.F90 error.o
driver.o        : driver.F90 config.o evolution.o io.o mesh.o timer.o
error.o         : error.F90
evolution.o     : evolution.F90 blocks.o boundaries.o config.o mesh.o scheme.o
interpolation.o : interpolation.F90
io.o            : io.F90 blocks.o error.o
mesh.o          : mesh.F90 blocks.o config.o error.o interpolation.o problem.o
problem.o       : problem.F90 blocks.o
scheme.o        : scheme.F90 blocks.o config.o interpolation.o
timer.o         : timer.F90

#-------------------------------------------------------------------------------
