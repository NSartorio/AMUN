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
ifeq ($(DEBUG),Y)
FFLAGS  = -g -O3 -DDEBUG
endif

ifeq ($(PROFILE),Y)
FFLAGS  += -pg -DPROFILE
endif

ifeq ($(OMP),Y)
ifeq ($(PROFILE),Y)
FFLAGS += -openmp-profile -openmp_report=2 -auto
else
FFLAGS += -openmp -openmp_report=2 -auto
endif
endif # OMP == Y

ifeq ($(SIGNALS),Y)
FFLAGS  += -DSIGNALS
endif

ifeq ($(EXCEPTIONS),Y)
FFLAGS  += -DFPE
endif

ifeq ($(DBLPREC),Y)
FFLAGS  += -r8 -DDBLPREC
endif

ifeq ($(MPI),Y)
FFLAGS  += -DMPI
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

ifeq ($(DIMS),D3)
ifeq ($(SLICES),Y)
FFLAGS  += -DSLICES
endif
endif

ifeq ($(KMHD),Y)
FFLAGS  += -DKMHD
ifeq ($(KMHDEN),Y)
FFLAGS  += -DKMHDEN
endif
MHD        = Y
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

ifeq ($(FORCING),Y)
FFLAGS  += -DFORCE

ifeq ($(FSTATS),Y)
ifeq ($(FORCERK),N)
FFLAGS  += -DFSTATS
endif
endif

ifeq ($(FCORR),Y)
FFLAGS  += -DFCORR
endif

ifeq ($(FORCERK),Y)
FFLAGS  += -DFORCERK
endif

ifeq ($(SUBVOL),Y)
FFLAGS  += -DSUBVOL
endif

else
LARGESCALE = N
endif

ifeq ($(LARGESCALE),Y)
FFLAGS  += -DLARGESCALE
endif

ifeq ($(CHOFOR),Y)
FFLAGS  += -DCHO
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

ifeq ($(SHAPE),Y)
FFLAGS  += -DSHAPE
endif

ifeq ($(MOVINGFRAME),Y)
FFLAGS  += -DMOVEFRAME

ifeq ($(ROTATEFRAME),Y)
FFLAGS  += -DROTATEFRAME
endif
endif

ifeq ($(ORD5TH),Y)
FFLAGS  += -DHIORDER -DORD5TH
ORD3RD   = N
endif

ifeq ($(ORD3RD),Y)
FFLAGS  += -DHIORDER -DORD3RD
endif

ifeq ($(LIMITER),LF)
FFLAGS  += -DLFLIM
endif

ifeq ($(POSPRE),Y)
FFLAGS  += -DPOSPRE
endif

ifeq ($(LIMITMAX),Y)
FFLAGS  += -DLIMITMAX
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

ifeq ($(CONSRECON),Y)
FFLAGS  += -DCONREC
endif

ifeq ($(SOURCES),Y)
FFLAGS  += -DSOURCES

ifeq ($(SOURCERK),Y)
FFLAGS  += -DSOURCERK
endif

ifeq ($(GRAVITY),Y)
FFLAGS  += -DGRAV
endif

ifeq ($(SELFGRAVITY),Y)
FFLAGS  += -DSELFGRAV
ifeq ($(DBLPREC),Y)
LIBS    += -lfftw3_threads -lfftw3
else
LIBS    += -lfftw3f_threads -lfftw3f
endif

endif

ifeq ($(VISCOSITY),Y)
FFLAGS  += -DVISCOSITY
endif

ifeq ($(MHD),Y)
ifeq ($(RESISTIVITY),Y)
FFLAGS  += -DRESISTIVITY

ifeq ($(RESIST_ANOM),Y)
FFLAGS  += -DRESIST_ANOM
endif
endif
endif
endif

ifeq ($(PARTICLES),Y)
FFLAGS  += -DPARTICLES
ifeq ($(NEAREST),Y)
FFLAGS  += -DNEAREST
endif

ifeq ($(PACCEL),Y)
FFLAGS  += -DPACCEL
endif
endif

ifeq ($(CVENDOR),PGI)
FFLAGS  += -DPGI
endif

ifeq ($(USERINT),Y)
ifeq ($(MHD),Y)
FFLAGS  += -DUSERINT
endif
endif

ifeq ($(STATIC),Y)
LDFLAGS  += -static
endif

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = godunov-amr

default: $(name).x

sources = blocks.F90 config.F90 driver.F90 error.F90
objects = blocks.o   config.o   driver.o   error.o
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

blocks.o        : blocks.F90
config.o        : config.F90 error.o
driver.o        : driver.F90 blocks.o config.o
error.o         : error.F90

#-------------------------------------------------------------------------------
