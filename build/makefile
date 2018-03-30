#-------------------------------------------------------------------------------
#
# directories indicating where the source files are, where the object files
# should be created and where to save the resulting executable file

SRCSDIR := ../src
OBJSDIR := ./obj
DESTDIR := .

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
# generate new makedeps file if does not exists
#
$(info Generating dependencies.)
$(shell ./mkdeps.sh $(SRCSDIR) $(OBJSDIR) > makedeps)
$(shell sleep 1)

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

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

#-------------------------------------------------------------------------------

name    = amun

modules = algebra blocks boundaries constants coordinates domains driver       \
          equations error evolution gravity integrals interpolations io mesh   \
          mpitools operators parameters problems random refinement schemes     \
          shapes sources timers user_problem utils

sources := $(addprefix $(SRCSDIR)/,$(addsuffix .F90, $(modules)))
objects := $(addprefix $(OBJSDIR)/,$(addsuffix .o, $(modules)))

all: $(name).x

$(name).x: $(deps) $(objects) | $(DESTDIR)
	$(LD) $(LDFLAGS) $(objects) $(LIBS) -o $(DESTDIR)/$(name).x

$(OBJSDIR)/%.o : $(SRCSDIR)/%.F90
	$(FC) -c $(FFLAGS) -J $(OBJSDIR) $< -o $@

$(objects): | $(OBJSDIR)

$(OBJSDIR):
	mkdir -p $(OBJSDIR)

$(DESTDIR):
	mkdir -p $(DESTDIR)

clean:
	rm -rf $(OBJSDIR) $(DESTDIR)/$(name).x

#-------------------------------------------------------------------------------
#
# object files dependencies (generated using mkdeps.sh script)
#
include makedeps

#-------------------------------------------------------------------------------
