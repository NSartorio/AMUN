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
# check flag conditions
#
ifeq ($(EQS),HYDRO)
ifeq ($(FLUX),HLLD)
FLUX = HLLC
endif
ifeq ($(EOS),ISO)
FLUX = HLL
endif
endif
ifeq ($(EQS),MHD)
ifeq ($(FLUX),HLLC)
FLUX = HLLD
endif
endif

#-------------------------------------------------------------------------------
#
# pass additional flags to the compiler
#
# compiler
#
FFLAGS  += -D${COMPILER}

# number of dimensions
#
FFLAGS  += -DNDIMS=${NDIMS}
ifeq ($(SHAPE),Y)
FFLAGS  += -DSHAPE
endif
ifeq ($(REFINE),Y)
FFLAGS  += -DREFINE
endif

# check if reconstruction is TVD
#
ifeq ($(SPACE),MINMOD)
FFLAGS  += -DTVD
endif
ifeq ($(SPACE),LF)
FFLAGS  += -DTVD
endif
ifeq ($(SPACE),MP5)
FFLAGS  += -DMP
endif
ifeq ($(SPACE),MP7)
FFLAGS  += -DMP
endif
ifeq ($(SPACE),MP9)
FFLAGS  += -DMP
endif

# output data format
#
FFLAGS  += -D${OUTPUT}

# compression
#
ifneq ($(COMPRESS),NONE)
FFLAGS  += -DCOMPRESS -D${COMPRESS}
endif

# precision
#
FFLAGS  += -DPREC=${PREC}

# time integration method
#
FFLAGS  += -D${TIME}

# spacial reconstruction method
#
FFLAGS  += -D${SPACE}

# set of equations
#
FFLAGS  += -D${EQS}

# equation of state
#
FFLAGS  += -D${EOS}

# fully conservative scheme
#
ifeq ($(CONSERVATIVE),Y)
FFLAGS  += -DCONSERVATIVE
endif

# flux approximation method
#
FFLAGS  += -D${FLUX}

# electromotive force approximation method
#
ifeq ($(EQS),MHD)
FFLAGS  += -D${FLUXEMF}
endif

# source terms
#
ifeq ($(FORCING),Y)
FFLAGS  += -DFORCE
endif
ifeq ($(VISCOSITY),Y)
FFLAGS  += -DVISCOSITY
endif
ifeq ($(EQS),MHD)
ifeq ($(RESISTIVITY),Y)
FFLAGS  += -DRESISTIVITY
endif
endif

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = amun

default: $(name).x

sources = blocks.F90 boundaries.F90 config.F90 constants.F90 coordinates.F90   \
          driver.F90 equations.F90 error.F90 evolution.F90 forcing.F90         \
          integrals.F90 interpolations.F90 io.F90 mesh.F90 mpitools.F90        \
          parameters.F90 problem.F90 random.F90 scheme.F90 timers.F90          \
          variables.F90
objects = blocks.o boundaries.o config.o constants.o coordinates.o driver.o    \
          equations.o error.o evolution.o forcing.o integrals.o                \
          interpolations.o io.o mesh.o mpitools.o parameters.o problem.o       \
          random.o scheme.o timers.o variables.o
files   = $(sources) makefile make.default config.in license.txt hosts

$(name).x: $(objects)
	$(LD) $(LDFLAGS) $(objects) $(LIBS) -o $(name).x

arch: $(files)
	tar czvf $(name)-$(DATE).tar.gz $(files)

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

blocks.o         : blocks.F90 error.o
boundaries.o     : boundaries.F90 blocks.o config.o error.o interpolations.o   \
                   mpitools.o variables.o
config.o         : config.F90 error.o
constants.o      : constants.F90
coordinates.o    : coordinates.F90 parameters.o
driver.o         : driver.F90 blocks.o config.o coordinates.o equations.o      \
                   evolution.o forcing.o integrals.o interpolations.o io.o     \
                   mesh.o mpitools.o parameters.o random.o
equations.o      : equations.F90 parameters.o variables.o
error.o          : error.F90
evolution.o      : evolution.F90 blocks.o boundaries.o config.o coordinates.o  \
                   forcing.o interpolations.o mesh.o mpitools.o problem.o      \
                   scheme.o variables.o
forcing.o        : forcing.F90 config.o constants.o coordinates.o mpitools.o   \
                   random.o variables.o
integrals.o      : integrals.F90 blocks.o config.o coordinates.o evolution.o   \
                   mpitools.o variables.o
interpolations.o : interpolations.F90 blocks.o coordinates.o parameters.o      \
                   variables.o
io.o             : io.F90 blocks.o coordinates.o error.o evolution.o           \
                   mpitools.o random.o scheme.o variables.o
mesh.o           : mesh.F90 blocks.o coordinates.o error.o interpolations.o    \
                   mpitools.o problem.o variables.o
mpitools.o       : mpitools.F90 timers.o
parameters.o     : parameters.F90 mpitools.o
problem.o        : problem.F90 blocks.o constants.o coordinates.o mpitools.o   \
                   random.o scheme.o variables.o
scheme.o         : scheme.F90 blocks.o coordinates.o interpolations.o          \
                   variables.o
random.o         : random.F90 mpitools.o parameters.o
timers.o         : timers.F90
variables.o      : variables.F90

#-------------------------------------------------------------------------------
