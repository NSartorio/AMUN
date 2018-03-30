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
# directories indicating where the source files are, where the object files
# should be created and where to save the resulting executable file

SRCSDIR := ../src
OBJSDIR := ./obj
DESTDIR := .

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

name = amun

modules = algebra blocks boundaries constants coordinates domains driver       \
          equations error evolution gravity integrals interpolations io mesh   \
          mpitools operators parameters problems random refinement schemes     \
          shapes sources timers user_problem utils

sources := $(addprefix $(SRCSDIR)/,$(addsuffix .F90, $(modules)))
objects := $(addprefix $(OBJSDIR)/,$(addsuffix .o, $(modules)))

all: $(name).x

clean:
	rm -rf $(OBJSDIR) $(DESTDIR)/$(name).x

$(OBJSDIR)/%.o : $(SRCSDIR)/%.F90 makefile
	$(FC) -c $(FFLAGS) -J $(OBJSDIR) $< -o $@

$(name).x: $(objects) | $(DESTDIR)
	$(LD) $(LDFLAGS) $(objects) $(LIBS) -o $(DESTDIR)/$(name).x

$(objects): | $(OBJSDIR)

$(OBJSDIR):
	mkdir -p $(OBJSDIR)

$(DESTDIR):
	mkdir -p $(DESTDIR)

#-------------------------------------------------------------------------------
#
# object files dependencies
#
algebra_deps        = constants error
blocks_deps         = error timers
boundaries_deps     = blocks coordinates equations error gravity               \
                      interpolations mpitools timers user_problem
constants_deps      =
coordinates_deps    = parameters
driver_deps         = blocks coordinates equations evolution gravity integrals \
                      interpolations io mesh mpitools operators parameters     \
                      problems random refinement schemes shapes sources        \
                      user_problem
equations_deps      = algebra coordinates error parameters timers
error_deps          =
evolution_deps      = blocks boundaries coordinates equations mesh mpitools    \
                      parameters schemes shapes sources
domains_deps        = blocks boundaries coordinates parameters
gravity_deps        = parameters timers user_problem
integrals_deps      = blocks coordinates equations error evolution mpitools    \
                      parameters timers
interpolations_deps = algebra blocks coordinates error parameters timers
io_deps             = blocks coordinates equations error evolution mesh        \
                      mpitools random refinement timers
mesh_deps           = blocks coordinates domains equations error               \
                      interpolations mpitools problems refinement timers
mpitools_deps       = error timers
operators_deps      = timers
parameters_deps     = mpitools
problems_deps       = blocks constants coordinates equations error parameters  \
                      random timers user_problem
random_deps         = mpitools parameters
refinement_deps     = blocks coordinates equations operators parameters timers
schemes_deps        = algebra coordinates equations interpolations timers
shapes_deps         = blocks constants coordinates equations parameters timers \
                      user_problem
sources_deps        = blocks coordinates equations gravity operators           \
                      parameters timers user_problem
timers_deps         =
user_problem_deps   = blocks constants coordinates equations error operators   \
                      parameters random timers
utils_deps          = error

$(OBJSDIR)/algebra.o        : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(algebra_deps)))
$(OBJSDIR)/blocks.o         : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(blocks_deps)))
$(OBJSDIR)/boundaries.o     : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(boundaries_deps)))
$(OBJSDIR)/constants.o      : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(constants_deps)))
$(OBJSDIR)/coordinates.o    : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(coordinates_deps)))
$(OBJSDIR)/driver.o         : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(driver_deps)))
$(OBJSDIR)/equations.o      : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(equations_deps)))
$(OBJSDIR)/error.o          : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(error_deps)))
$(OBJSDIR)/evolution.o      : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(evolution_deps)))
$(OBJSDIR)/domains.o        : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(domains_deps)))
$(OBJSDIR)/gravity.o        : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(gravity_deps)))
$(OBJSDIR)/integrals.o      : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(integrals_deps)))
$(OBJSDIR)/interpolations.o : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(interpolations_deps)))
$(OBJSDIR)/io.o             : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(io_deps)))
$(OBJSDIR)/mesh.o           : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(mesh_deps)))
$(OBJSDIR)/mpitools.o       : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(mpitools_deps)))
$(OBJSDIR)/operators.o      : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(operators_deps)))
$(OBJSDIR)/parameters.o     : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(parameters_deps)))
$(OBJSDIR)/problems.o       : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(problems_deps)))
$(OBJSDIR)/random.o         : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(random_deps)))
$(OBJSDIR)/refinement.o     : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(refinement_deps)))
$(OBJSDIR)/schemes.o        : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(schemes_deps)))
$(OBJSDIR)/shapes.o         : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(shapes_deps)))
$(OBJSDIR)/sources.o        : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(sources_deps)))
$(OBJSDIR)/timers.o         : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(timers_deps)))
$(OBJSDIR)/user_problem.o   : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(user_problem_deps)))
$(OBJSDIR)/utils.o          : $(addprefix $(OBJSDIR)/,$(addsuffix .o,          \
                                $(utils_deps)))

#-------------------------------------------------------------------------------
