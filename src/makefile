#-------------------------------------------------------------------------------
#
# include configuration
#
$(info )
ifeq ($(wildcard make.config),make.config)
$(info Using customized file 'make.config'.)
$(info )
include make.config
else
$(warning Could not find customized file 'make.config'!)
$(info )
$(info File 'make.default' is an example with all available compilation time   \
       options.  You can modify it, but it is better to copy it to a new file  \
       'make.config' and then customize it.  This will also make disappear     \
       the following 15 second pause every time you compile.)
$(info )
$(shell sleep 15)
include make.default
endif

#-------------------------------------------------------------------------------
#
# host name
#
HOST := $(shell hostname)

include hosts/default
ifeq ($(wildcard hosts/$(HOST)),hosts/$(HOST))
$(info Using customized compiler setup from 'hosts/$(HOST)' file.)
$(info )
include hosts/$(HOST)
else
$(info Using default compiler setup from file 'hosts/default'.                 \
       Do not modify it!)
$(info Instead, copy this file to 'hosts/$(HOST)' and customize compilers and  \
       compilation options there.  This will also get rid of the following     \
       15 second pause every time you compile.)
$(info )
$(shell sleep 15)
endif

#-------------------------------------------------------------------------------
#
# check flag conditions
#

#-------------------------------------------------------------------------------
#
# pass additional flags to the compiler
#
# compiler
#
FFLAGS  += ${CPPPREFIX}-D${COMPILER}

# number of dimensions
#
FFLAGS  += ${CPPPREFIX}-DNDIMS=${NDIMS}

# output data format
#
FFLAGS  += ${CPPPREFIX}-D${OUTPUT}

# compression
#
ifneq ($(COMPRESS),NONE)
FFLAGS  += ${CPPPREFIX}-DCOMPRESS ${CPPPREFIX}-D${COMPRESS}
endif

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = amun

default: $(name).x

sources = algebra.F90 blocks.F90 boundaries.F90 constants.F90 coordinates.F90  \
          domains.F90 driver.F90 equations.F90 error.F90 evolution.F90         \
          integrals.F90 interpolations.F90 io.F90 mesh.F90 mpitools.F90        \
          operators.F90 parameters.F90 problems.F90 random.F90 refinement.F90  \
          schemes.F90 shapes.F90 sources.F90 timers.F90 utils.F90
objects = algebra.o blocks.o boundaries.o constants.o coordinates.o domains.o  \
          driver.o equations.o error.o evolution.o integrals.o                 \
          interpolations.o io.o mesh.o mpitools.o operators.o parameters.o     \
          problems.o random.o refinement.o schemes.o shapes.o sources.o        \
          timers.o utils.o

$(name).x: $(objects)
	$(LD) $(LDFLAGS) $(objects) $(LIBS) -o $(name).x

clean-bak:
	rm -rf *.bak *~

clean-data:
	rm -rf *.dat *.h5

clean-exec:
	rm -rf *.x

clean-logs:
	rm -rf *.log

clean-modules:
	rm -rf *.mod

clean-objects:
	rm -rf *.o

clean: clean-bak clean-modules clean-objects

clean-all: clean-bak clean-data clean-exec clean-logs clean-modules            \
           clean-objects

#-------------------------------------------------------------------------------

algebra.o        : algebra.F90 constants.o
blocks.o         : blocks.F90 error.o timers.o
boundaries.o     : boundaries.F90 blocks.o coordinates.o equations.o error.o   \
                   interpolations.o mpitools.o timers.o
constants.o      : constants.F90
coordinates.o    : coordinates.F90 parameters.o
driver.o         : driver.F90 blocks.o coordinates.o equations.o evolution.o   \
                   integrals.o interpolations.o io.o mesh.o mpitools.o         \
                   operators.o parameters.o problems.o random.o refinement.o   \
                   schemes.o shapes.o sources.o
equations.o      : equations.F90 algebra.o coordinates.o parameters.o timers.o
error.o          : error.F90
evolution.o      : evolution.F90 blocks.o boundaries.o coordinates.o           \
                   equations.o mesh.o mpitools.o parameters.o schemes.o        \
                   shapes.o sources.o
domains.o        : domains.F90 blocks.o boundaries.o coordinates.o parameters.o
integrals.o      : integrals.F90 blocks.o coordinates.o equations.o error.o    \
                   evolution.o mpitools.o parameters.o timers.o
interpolations.o : interpolations.F90 blocks.o coordinates.o error.o           \
                   parameters.o timers.o
io.o             : io.F90 blocks.o coordinates.o equations.o error.o           \
                   evolution.o mesh.o mpitools.o random.o refinement.o timers.o
mesh.o           : mesh.F90 blocks.o coordinates.o domains.o equations.o       \
                   error.o interpolations.o mpitools.o problems.o refinement.o \
                   timers.o
mpitools.o       : mpitools.F90 error.o timers.o
operators.o      : operators.F90 timers.o
parameters.o     : parameters.F90 mpitools.o
problems.o       : problems.F90 blocks.o constants.o coordinates.o equations.o \
                   error.o operators.o parameters.o random.o timers.o
random.o         : random.F90 mpitools.o parameters.o
refinement.o     : refinement.F90 blocks.o coordinates.o equations.o           \
                   operators.o parameters.o timers.o
schemes.o        : schemes.F90 algebra.o coordinates.o equations.o             \
                   interpolations.o timers.o
shapes.o         : shapes.F90 blocks.o constants.o coordinates.o equations.o   \
                   parameters.o timers.o
sources.o        : sources.F90 blocks.o coordinates.o equations.o operators.o  \
                   parameters.o timers.o
timers.o         : timers.F90
utils.o          : utils.F90 error.o

#-------------------------------------------------------------------------------
