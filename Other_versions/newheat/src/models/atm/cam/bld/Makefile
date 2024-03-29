#----------------------------------------------------------------------
# Makefile for building CAM on various platforms.
#
# This makefile assumes the existence the file:
#
# Filepath      The directories, listed one per line, that contain the source
#               code required to build CAM.  This list is used to set the
#               VPATH variable which is used by GNU make to search for
#               dependencies (after it looks in the directory from which
#               it was invoked).  This list of directories, prepended with ./,
#               is also used to construct the list of search directories 
#               used by the preprocessor (as specified by -I command-line options).
#
# The following macros may be set in the user's environment:
#

# ROOTDIR       Root directory for the CAM source distribution.  If ROOTDIR is not
#               set the makefile searches for it in a file called "Rootdir".
# EXENAME       Name to call the executable.  Default: atm
# MODEL_EXEDIR  Directory to build the executable in.  Default: ./
# INC_NETCDF    Directory containing the NetCDF include files.  Default: /usr/local/include
# LIB_NETCDF    Directory containing the NetCDF library.  Default: /usr/local/lib
# MOD_NETCDF    Directory containing the NetCDF module files.  Default: /usr/local/lib
# INC_MPI       Directory containing the MPI include files.  Default: /usr/local/include
# LIB_MPI       Directory containing the MPI library.  Default: /usr/local/lib
# ESMF_ROOT     Root directory of ESMF source tree.  Default: $(ROOTDIR)/models/utils/esmf
# ESMF_BLD      Upper level build directory for ESMF
#               The ESMF .o files are built in $(ESMF_BLD)/obj
#               The ESMF .m files are built in arch dependent subdirs of $(ESMF_BLD)/mod
#               The ESMF library is built in architecture and optimization
#               dependent subdirectories of $(ESMF_BLD)/lib
#               Default: ./esmf
# CARMA_ROOT    Root directory of CARMA source tree.  Default: $(ROOTDIR)/models/utils/carma  EJL - changed for Titan
# CARMA_MAKE    Name of CARMA makefile.  Default: Makefile.darwin.intel
# DEBUG         Set to TRUE to turn on compiler debugging options.  Default: FALSE 
# SPMD          Whether to build in SPMD mode or not. 	[values TRUE FALSE]
#               Default: read from ./misc.h
#               N.B. To ensure consistency between the Makefile and misc.h the SPMD
#                    macro should not be set.  This forces it to be read from misc.h.
# SMP           Set to TRUE to enable building in SMP mode (uses OpenMP).  
#               Currently implemented for IBM, SGI, linux-pgf90. (default is TRUE on IBM and 
#               linux-pgf90, and depends on SPMD setting on SGI).
# NESTED_OMP    Set to TRUE to enable nested OpenMP support.  AIX only.
# USER_FC       Allow user to override the default Fortran compiler specified in Makefile.
# USER_CC       Allow user to override the default C compiler specified in Makefile (linux only).
# USER_CPPDEFS  Additional CPP defines.
# USER_CFLAGS   Additional C compiler flags that the user wishes to set.
# USER_FFLAGS   Additional Fortran compiler flags that the user wishes to set.
# USER_LDLAGS   Additional load flags that the user wishes to set.
# F_OPTIMIZATION_OVERRIDE
#               Use this to replace default Fortran compiler optimization.  
# NO_SWITCH     On Compaq if the hardware switch is not available
#               set this env variable to "TRUE".
#
# Note: The ESMF library is included in the CAM distribution in 
#       $ROOTDIR/models/utils/esmf and is built using this makefile.  
# Note: The CARMA library is included in the CAM distribution in 
#       $ROOTDIR/models/utils/carma and is built using this makefile.  
#------------------------------------------------------------------------

# Set up special characters
null  :=
space := $(null) $(null)
comma := $(null),$(null)

# Determine distribution root directory
ifeq ($(ROOTDIR),$(null))
ROOTDIR := /home/larsonej/titancarma/src
endif

# Check for override of default Fortran compiler optimizations
ifneq ($(F_OPTIMIZATION_OVERRIDE),$(null))
FORTRAN_OPTIMIZATION := $(F_OPTIMIZATION_OVERRIDE)
endif

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /Volumes/Data/Libraries/netcdf/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /Volumes/Data/Libraries/netcdf/include
endif

ifeq ($(MOD_NETCDF),$(null))
MOD_NETCDF := $(LIB_NETCDF)
endif

# Check for the MPI library and include directories 
ifeq ($(LIB_MPI),$(null))
LIB_MPI := /Volumes/Data/Libraries/mpich/lib
endif

ifeq ($(INC_MPI),$(null))
INC_MPI := /Volumes/Data/Libraries/mpich/include
endif

# Build the ESMF library
ifeq ($(ESMF_ROOT),$(null))
ESMF_ROOT := ${ROOTDIR}/models/utils/esmf
endif

ifeq ($(ESMF_BLD),$(null))
ESMF_BLD := $(shell (mkdir ./esmf 2>/dev/null;cd ./esmf;pwd))
endif

ESMF_BOPT := O
ifeq ($(DEBUG),TRUE)
  ESMF_BOPT := g
endif
ESMF_MOD := $(ESMF_BLD)/mod/mod$(ESMF_BOPT)
ESMF_LIB := $(ESMF_BLD)/lib/lib$(ESMF_BOPT)

# Build the CARMA libraries
#ifeq ($(CARMA_ROOT),$(null))
CARMA_ROOT := ${ROOTDIR}/models/utils/carma
#endif

ifeq ($(CARMA_ROOT),$(null))
CARMA_MAKE := Makefile.linux.intel
#Makefile.linux.mpif90
endif

ifeq ($(CARMA_BLD),$(null))
CARMA_BLD := $(shell (mkdir ./carma 2>/dev/null;cd ./carma;pwd))
endif

ifeq ($(INC_CARMA),$(null))
INC_CARMA := $(CARMA_ROOT)/aer $(CARMA_ROOT)/rad $(CARMA_ROOT)/include $(CARMA_ROOT)/include/cam
endif

CARMA_BOPT := O
ifeq ($(DEBUG),TRUE)
  CARMA_BOPT := g
endif

TRICK := $(shell rm $(CARMA_BLD)/aer/libaer.a $(CARMA_BLD)/rad/librad.a)

# Check for directory in which to put executable
ifeq ($(MODEL_EXEDIR),$(null))
MODEL_EXEDIR := .
endif

# Check for name of executable
ifeq ($(EXENAME),$(null))
EXENAME := atm
endif

# Check if SPMD is defined in "misc.h"
# Ensure that it is defined and not just "undef SPMD" set in file
ifeq ($(SPMD),$(null))
	SPMDSET := $(shell grep SPMD misc.h)
	ifneq (,$(findstring define,$(SPMDSET)))
		SPMD := TRUE
	else
		SPMD := FALSE
	endif
endif

# Load dependency search path.
dirs := . $(shell cat Filepath)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs) $(INC_NETCDF) $(INC_MPI) $(INC_CARMA)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
VPATH    := $(subst $(space),:,$(VPATH))               

#------------------------------------------------------------------------
# Primary target: build the model
#------------------------------------------------------------------------
all: $(MODEL_EXEDIR)/$(EXENAME)

# Get list of files and build dependency file for all .o files
#   using perl scripts mkSrcfiles and mkDepends

SOURCES   := $(shell cat Srcfiles)

# Newer makes set the CURDIR variable.
CURDIR := $(shell pwd)

$(CURDIR)/Depends: $(CURDIR)/Srcfiles $(CURDIR)/Filepath
	$(ROOTDIR)/models/atm/cam/bld/mkDepends Filepath Srcfiles > $@

$(CURDIR)/Srcfiles: $(CURDIR)/Filepath
	$(ROOTDIR)/models/atm/cam/bld/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES)))

CPPDEF    := -DCAM -DNO_SHR_VMATH -DHIDE_SHR_MSG $(USER_CPPDEFS)
 
$(MODEL_EXEDIR)/$(EXENAME): $(OBJS) $(CARMA_BLD)/aer/libaer.a $(CARMA_BLD)/rad/librad.a 
	$(FC) -o $@ $(OBJS) -L$(LIB_HDF5) -lhdf5_hl -lhdf5 -L$(LIB_HDF4) -lmfhdf -ldf -L$(LIB_NETCDF) -lnetcdff -lnetcdf -L$(ESMF_LIB)/$(ESMF_ARCH) -lesmf -L$(CARMA_BLD)/aer -laer -L$(CARMA_BLD)/rad -lrad $(LDFLAGS)

debug: $(OBJS)
	echo "FFLAGS: $(FFLAGS)"
	echo "LDFLAGS: $(LDFLAGS)"
	echo "OBJS: $(OBJS)"

test_fc: test_fc.o
	$(FC) -o $@ test_fc.o $(LDFLAGS)
test_nc: test_nc.o
	$(FC) -o $@ test_nc.o -L$(LIB_NETCDF) -lnetcdff -lnetcdf $(LDFLAGS)
test_mpi: test_mpi.o
	$(FC) -o $@ test_mpi.o $(LDFLAGS)

# Architecture-specific flags and rules
#
# Determine platform 
ifeq ($(UNAMES),$(null))
  UNAMES := $(shell uname -s)
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)

ifeq ($(SMP),$(null))
  SMP := TRUE
endif

#ADDRESS    := Q32
 ADDRESS    := Q64
ifeq ($(ADDRESS),Q32)
  QCMP      :=
  QLDR      := -bmaxdata:0x80000000 -bmaxstack:0x10000000
  ESMF_ARCH := rs6000_sp
endif
ifeq ($(ADDRESS),Q64)
  QCMP      := -q64
  QLDR      := -q64
  ESMF_ARCH := rs6000_64
endif

CC         := xlc_r
CPPDEF     += -DAIX
ifeq ($(NESTED_OMP),TRUE)
  CPPDEF   += -DNESTED_PAR
endif
CFLAGS     := $(cpp_path) -O2 $(CPPDEF) $(QCMP)

cpre = $(null)-WF,-D$(null)
FPPFLAGS   := $(patsubst -D%,$(cpre)%,$(CPPDEF))

mod_path   := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
FFLAGS     := $(cpp_path) $(mod_path) $(FPPFLAGS) -qarch=auto -qrealsize=8 -qdpc=e \
		-qspillsize=2500 -g -qfullpath $(QCMP)
FREEFLAGS  := -qsuffix=f=f90:cpp=F90
FIXEDFLAGS := -qfixed=132
#
# As of Apr/15/2003 cam2_0_2_dev18 xlfrte8.1.0.3 using -lmass causes the code to core-dump
# when using DEBUG compiler options.
#
# ... hence the following ugly expedient ...
ifeq ($(DEBUG),TRUE)
  LDFLAGS := $(QLDR)
else
  LDFLAGS := $(QLDR)
  LDFLAGS += -lmass
endif

# Use the thread-safe compilers even when disabling openMP.
ifeq ($(SPMD),TRUE)
  FC := mpxlf90_r
else
  FC := xlf90_r
  FFLAGS += -WF,-DHIDE_MPI
endif

ifeq ($(SMP),TRUE)
  # THREADED_PTHREADS is used by the timing library
  CFLAGS += -DTHREADED_PTHREADS
endif

ifeq ($(DEBUG),TRUE)
#
# Bounds checking is unreliable on the IBM.
# Sometimes you can get it to go if you turn threading off (by deleting -qsmp=omp)
# Only turn float-trapping on for debug mode as it's a 20% performance hit.
#
  FFLAGS += -qinitauto=FF911299 -qflttrap=ov:zero:inv:en
  ifeq ($(SMP),TRUE)
    ifeq ($(NESTED_OMP),TRUE)
      FFLAGS += -qsmp=omp:noopt:nested_par
      LDFLAGS += -qsmp=omp:noopt:nested_par
    else
      FFLAGS += -qsmp=omp:noopt
      LDFLAGS += -qsmp=omp:noopt
    endif
  endif
else

# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
# Inline when not debugging
    FORTRAN_OPTIMIZATION := -O3 -qstrict -Q
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
  ifeq ($(SMP),TRUE)
    ifeq ($(NESTED_OMP),TRUE)
      FFLAGS += -qsmp=omp:nested_par
      LDFLAGS += -qsmp=omp:nested_par
    else
      FFLAGS += -qsmp=omp
      LDFLAGS += -qsmp=omp
    endif
  endif
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<

.c.o:
	cc -c  $(CFLAGS) $<


endif

#------------------------------------------------------------------------
# DARWIN
#------------------------------------------------------------------------

ifeq ($(UNAMES),Darwin)

  ifeq ($(SMP),$(null))
    SMP := TRUE
  endif

  # Set the Fortran compilerq
  ifeq ($(USER_FC),$(null))
    FC := xlf90_r
  else
    FC := $(USER_FC)
  endif

  ifeq ($(USER_CC),$(null))
    CC := gcc
  else
    CC := $(USER_FC)
  endif

  ifeq ($(findstring ifort,$(FC)),ifort)
    ADDRESS    := Q32
#    ADDRESS    := Q64
    ifeq ($(ADDRESS),Q32)
      QCMP      :=
      ESMF_ARCH := Darwin_intel
      CARMA_MAKE := Makefile.darwin.intel
    endif
  
    ifeq ($(ADDRESS),Q64)
      QCMP      := -m64
      QLDR      := -m64
      ESMF_ARCH := Darwin_intel_64 
    endif

    CPPDEF      += -DFORTRANUNDERSCORE -DLINUX -DDARWIN
    mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
    
    # -mp is needed for NAN detection of physics buffer in the CARMA init code
    FFLAGS      := $(cpp_path) $(mod_path) $(CPPDEF) -132 -autodouble -ftz -mp
    SPEC_FFLAGS := $(FFLAGS)
    LDFLAGS     :=
    FREEFLAGS   := -FR

    ifeq ($(DEBUG),TRUE)
      FFLAGS += -CB -O0 -g
    else
      # Check for override of default Fortran compiler optimizations
      ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
	      FORTRAN_OPTIMIZATION := -O2
      endif
      FFLAGS += $(FORTRAN_OPTIMIZATION)
    endif

    ifeq ($(SMP),TRUE)
      FFLAGS += -openmp
      LDFLAGS += -openmp
    endif

    ifeq ($(SPMD),TRUE)
      LDFLAGS += -L$(LIB_MPI) -l$(MPI_LIB_NAME)f90 -lp$(MPI_LIB_NAME) -l$(MPI_LIB_NAME)
    else
      FFLAGS  += -DHIDE_MPI
    endif
  endif

  ifeq ($(findstring xlf,$(FC)),xlf)
    CARMA_MAKE := Makefile.darwin.xlf
    ESMF_ARCH       := Darwin_xlf
    ESMF_C_COMPILER := gcc

    ADDRESS    := Q32
    CPPDEF     += -DAIX -DDARWIN

    FPPFLAGS   := $(patsubst -D%,$(cpre)%,$(CPPDEF))

    mod_path   := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
    FFLAGS     := $(cpp_path) $(mod_path) $(FPPFLAGS) -qarch=auto -qrealsize=8 -qdpc=e \
 		  -qspillsize=2500 -g -qfullpath $(QCMP)
    FREEFLAGS  := -qsuffix=f=f90:cpp=F90
    FIXEDFLAGS := -qfixed=132

    # Use the thread-safe compilers even when disabling openMP.
    ifeq ($(SPMD),TRUE)
      FC := mpxlf90_r
    else
      FC := xlf90_r
      FFLAGS += -WF,-DHIDE_MPI
    endif

    ifeq ($(DEBUG),TRUE)
    #
    # Bounds checking is unreliable on the IBM.
    # Sometimes you can get it to go if you turn threading off (by deleting -qsmp=omp)
    # Only turn float-trapping on for debug mode as it's a 20% performance hit.
    #
      FFLAGS += -qinitauto=FF911299 -qflttrap=ov:zero:inv:en
      ifeq ($(SMP),TRUE)
	ifeq ($(NESTED_OMP),TRUE)
	  FFLAGS += -qsmp=omp:noopt:nested_par
	  LDFLAGS += -qsmp=omp:noopt:nested_par
	else
	  FFLAGS += -qsmp=omp:noopt
	  LDFLAGS += -qsmp=omp:noopt
	endif
      endif
    else

    # Check for override of default Fortran compiler optimizations
      ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
    # Inline when not debugging
	FORTRAN_OPTIMIZATION := -O3 -qstrict -Q
      endif
      FFLAGS += $(FORTRAN_OPTIMIZATION)
      ifeq ($(SMP),TRUE)
	ifeq ($(NESTED_OMP),TRUE)
	  FFLAGS += -qsmp=omp:nested_par
	  LDFLAGS += -qsmp=omp:nested_par
	else
	  FFLAGS += -qsmp=omp
	  LDFLAGS += -qsmp=omp
	endif
      endif
    endif
  endif

  ifeq ($(NESTED_OMP),TRUE)
    CPPDEF   += -DNESTED_PAR
  endif
  CFLAGS     := $(cpp_path) -O2 $(CPPDEF) $(QCMP)

  cpre = $(null)-WF,-D$(null)
  LDFLAGS += -lSystemStubs

  ifeq ($(SMP),TRUE)
    # THREADED_PTHREADS is used by the timing library
    CFLAGS += -DTHREADED_PTHREADS
  endif


.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<

.c.o:
	cc -c  $(CFLAGS) $<

endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)

# default is either pure SPMD or pure SMP parallelism
ifeq ($(SMP),$(null))
  ifeq ($(SPMD),TRUE)
    SMP := FALSE
  else
    SMP := TRUE
  endif
endif

ESMF_ARCH   := IRIX64
FC          := f90
CPP         := /lib/cpp

CPPDEF      += -DIRIX64
CFLAGS      := $(cpp_path) -64 -O2 $(CPPDEF)
mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
FFLAGS      := $(cpp_path) $(mod_path) -64 -r8 -i4 -c -cpp -extend_source \
		$(CPPDEF)
LDFLAGS     = -64 -mp

ifeq ($(SMP),TRUE)
  # THREADED_OMP is used by the timing library
  CFLAGS += -DTHREADED_OMP
# Only enable OpenMP, not SGI specific parallelism
# Note that -lmp should precede -lmpi when running hybrid MPI/OpenMP
  FFLAGS  += -mp -MP:dsm=OFF -MP:old_mp=OFF
  LDFLAGS += -mp -lmp
endif

ifeq ($(SPMD),TRUE)
  FFLAGS  += -I$(INC_MPI)
  LDFLAGS += -L$(LIB_MPI) -lmpi
else
  FFLAGS  += -DHIDE_MPI
endif

# WARNING: -mp and -g together has been known to cause wrong answers
ifeq ($(DEBUG),TRUE)
  FFLAGS  += -g -DEBUG:trap_uninitialized=ON -C
else
# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
    FORTRAN_OPTIMIZATION := -O2
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) $(FFLAGS) $<

.F90.o:
	$(FC) $(FFLAGS) $<
.c.o:
	cc -c  $(cpp_path) $(CFLAGS) $<


endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
ESMF_ARCH = solaris
FC	:= f90
FC77	:= f77
CPP     := /usr/ccs/lib/cpp
CFLAGS  := $(cpp_path) -DSUNOS $(CPPDEF)
mod_path:= -M$(ESMF_MOD)/$(ESMF_ARCH) -M$(MOD_NETCDF)
# Don't use OpenMP on Solaris as it currently causes problems
FFLAGS  := $(cpp_path) $(mod_path) -xs -stackvar -Qoption f90comp -r8const \
		-xtypemap=real:64,integer:32 -e -DSUNOS -DHIDE_MPI \
		$(CPPDEF)
LDFLAGS := -L/opt/SUNWspro/lib -lf77compat -openmp -fast
SPEC_FFLAGS := $(cpp_path) -r8 -i4 -c

ifeq ($(DEBUG),TRUE)
  FFLAGS += -g -dalign
  SPEC_FFLAGS += -g -dalign
else
# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
# Inline code when not debugging -inline
# Use symbol table (-g) as make execution stable and can find out where problems are.
# Turn off aggressive optimization options (fsimple=1, fns=no) so the climate
# will be similar to climate on the IBM
    FORTRAN_OPTIMIZATION := -inline=%auto -fast -g -fsimple=1 -fns=no -O4
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
endif

ifeq ($(SPMD),TRUE)
  FFLAGS  += -I$(INC_MPI)
  LDFLAGS += -L$(LIB_MPI) -lmpich -lnsl -lsocket
else
  FFLAGS  += -DHIDE_MPI
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F90.o:
	$(FC) -c $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)

CPPDEF += -DLINUX

# Set the Fortran compiler (default: PGI pgf90)
ifeq ($(USER_FC),$(null))
  FC := pgf90
else
  FC := $(USER_FC)
endif

# Decide whether or not a PGI Fortran compiler is being used by looking
# for a match of the string 'pgf' in the compiler name.
PGI_FC := FALSE
ifeq ($(findstring pgf,$(FC)),pgf)
  PGI_FC := TRUE
endif

# Set the C compiler (default: pgcc if using PGI Fortran compiler, cc otherwise)
# and cc if not
ifeq ($(USER_CC),$(null))
  ifeq ($(PGI_FC),TRUE)
    CC := pgcc
  else
    CC := cc
  endif
else
  CC := $(USER_CC)
endif

ifeq ($(SMP),TRUE)
  # THREADED_OMP and THREADED_PTHREADS are used by the timing library
  ifeq ($(CC),pgcc)
    CPPDEF += -DTHREADED_OMP
  else
    CPPDEF += -DTHREADED_PTHREADS
  endif
endif

CFLAGS = $(cpp_path) $(CPPDEF)
ifeq ($(CC),pgcc)
  CFLAGS += -fast
endif

# pgf90
# -DPGF90 is for phcs and gauaw which normally use r16 arithmetic but is unavailable under pgf90

ifeq ($(PGI_FC),TRUE)

  ifeq ($(CC),pgcc)
    ESMF_ARCH = linux_pgi
  else
    ESMF_ARCH = linux_gnupgf90
  endif

  ifeq ($(SMP),$(null))
    SMP := TRUE
  endif

  mod_path  := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
  FFLAGS    := $(cpp_path) $(mod_path) -r8 -i4 $(CPPDEF) -Mdalign -Mextend -DPGF90 -DNO_R16 -byteswapio
  FREEFLAGS := -Mfree
  LDFLAGS :=

  ifeq ($(DEBUG),TRUE)
    FFLAGS      += -g -Ktrap=fp -Mrecursive -Mbounds
    SPEC_FFLAGS := $(FFLAGS) 
  else
    SPEC_FFLAGS := $(FFLAGS) 
    # Check for override of default Fortran compiler optimizations
    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
      FORTRAN_OPTIMIZATION := -O1
    endif
    FFLAGS += $(FORTRAN_OPTIMIZATION)
  endif

  ifeq ($(SMP),TRUE)
    FFLAGS += -mp
    LDFLAGS += -mp
  endif

endif

# lf95
#
# Note that as of lf95 version 6.1 threading does NOT work because of
# ridiculously small per thread stacksize limits.
#
#  -CcdRR8 is an undocumented flag which promotes only vars declared "real", not "real*4"
# --trace produces a call traceback on abort
# --trap causes code to stop on divide by zero or overflow exceptions
# --pca prevents overwriting constant arguments
# --chk for basic compiler checking (a,e,s,u,x)

ifeq ($(FC),lf95)
  ESMF_ARCH   := linux_lf95
  mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
  FFLAGS      := $(cpp_path) $(mod_path) -CcdRR8 $(CPPDEF) --trace --trap --wide
  SPEC_FFLAGS := $(FFLAGS)
  LDFLAGS     :=
  ifeq ($(DEBUG),TRUE)
    #TBH:  this works    FFLAGS += -g --chk --pca
    #TBH:  this FAILS    FFLAGS += -g --chk a,e,s,u,x --pca
    FFLAGS += -g --chk a,e,s,u --pca
  else
    # Check for override of default Fortran compiler optimizations
    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
      FORTRAN_OPTIMIZATION := -O
    endif
    FFLAGS += $(FORTRAN_OPTIMIZATION)
  endif
endif

#ifort
ifeq ($(findstring ifort,$(FC)),ifort)
  ESMF_ARCH   := linux_intel
  mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
  FFLAGS      := $(cpp_path) $(mod_path) $(CPPDEF) -132 -autodouble -ftz -g -ip -no-prec-div
  SPEC_FFLAGS := $(FFLAGS)
  LDFLAGS     :=
  FREEFLAGS   := -FR
  MPILDFLAGS  := -lmpi_f90 -lmpi_f77 -lopen-rte -lopen-pal -ldl -Wl,-rpath -lnsl -lutil  
ifeq ($(DEBUG),TRUE)
    FFLAGS += -CB
  else
    # Check for override of default Fortran compiler optimizations
    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
      FORTRAN_OPTIMIZATION := -O2
    endif
    FFLAGS += $(FORTRAN_OPTIMIZATION)
  endif
  ifeq ($(SMP),TRUE)
    FFLAGS += -openmp
    LDFLAGS += -openmp
  endif
    CARMA_MAKE := Makefile.linux.intel
endif

# mpif90
   ifeq ($(findstring mpif90,$(FC)),mpif90)
     CPPFLAGS    := -D TITAN_FAO -D _DIML -D _TITAN
     FPPFLAGS    := -D TITAN_FAO -D _DIML -D _TITAN
     ESMF_ARCH   := linux_mpif90
     mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
     FFLAGS      := $(cpp_path) $(mod_path) $(CPPDEF) -132 -autodouble -ftz -ip -no-prec-div
     SPEC_FFLAGS := $(FFLAGS)
     LDFLAGS     := 
     MPILDFLAGS  := -lmpi_f90 -lmpi_f77 -lopen-rte -lopen-pal -ldl -Wl,-rpath -lnsl -lutil
#     MPILDFLAGS  := -lmpichf90 -lmpichf77 -lopen-rte -lopen-pal -ldl -Wl,-rpath -lnsl -lutil
     FREEFLAGS   := -FR
     ifeq ($(DEBUG),TRUE)
       FFLAGS += -CB -check all -check bounds
     else
       # Check for override of default Fortran compiler optimizations
       ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
         FORTRAN_OPTIMIZATION := -O2
       endif
       FFLAGS += $(FORTRAN_OPTIMIZATION)
     endif
     ifeq ($(SMP),TRUE)
       FFLAGS += -openmp
       LDFLAGS += -openmp
     endif
   	CARMA_MAKE := Makefile.linux.mpif90
   endif

# pathf90

ifeq ($(FC),pathf90)

  CC        := pathcc
  ESMF_ARCH := linux_pathscale

  CARMA_MAKE := Makefile.linux.pathscale
  
  mod_path  := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
  FFLAGS    := $(cpp_path) $(mod_path) $(CPPDEF) -r8 -extend_source -fno-second-underscore -DNO_R16
  FREEFLAGS :=
  LDFLAGS :=

  ifeq ($(DEBUG),TRUE)
    FFLAGS      += -g -O1 -C
    SPEC_FFLAGS := $(FFLAGS) 
  else
    SPEC_FFLAGS := $(FFLAGS) 
    # Check for override of default Fortran compiler optimizations
    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
      FORTRAN_OPTIMIZATION := -O
    endif
    FFLAGS += $(FORTRAN_OPTIMIZATION)
  endif

  ifeq ($(SMP),TRUE)
    FFLAGS +=
    LDFLAGS +=
  endif

endif


#ifeq ($(findstring ifort,$(FC)),ifort)
#  ESMF_ARCH   := linux_intel
#  mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
#  FFLAGS      := $(cpp_path) $(mod_path) $(CPPDEF) -132 -autodouble -ftz -g
#  SPEC_FFLAGS := $(FFLAGS)
#  LDFLAGS     :=
 # FREEFLAGS   := -FR
 # ifeq ($(DEBUG),TRUE)
  #  FFLAGS += -CB
#  else
#    # Check for override of default Fortran compiler optimizations
#    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
#      FORTRAN_OPTIMIZATION := -O2
#    endif
#    FFLAGS += $(FORTRAN_OPTIMIZATION)
#  endif
#  ifeq ($(SMP),TRUE)
#    FFLAGS += -openmp
#    LDFLAGS += -openmp
#  endif
#endif

# Flags common to all compilers

ifeq ($(SPMD),TRUE)
  LDFLAGS += -L$(LIB_MPI) -l$(MPI_LIB_NAME) $(MPILDFLAGS)
else
  FFLAGS  += -DHIDE_MPI
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

ifeq ($(FC),pgf90)
#
# To fix hanging problem when using sld dynamics, compile sgexx without "-fast"
#

sgexx.o: sgexx.F
	$(FC) -c $(SPEC_FFLAGS) $<
endif

ifeq ($(FC),lf95)
# lahey fails on binary_io due to writing wrap areas
binary_io.o: binary_io.F90
	$(FC) -c $(SPEC_FFLAGS) $<
wrap_nf.o: wrap_nf.F90
	$(FC) -c $(SPEC_FFLAGS) $<
wrap_mpi.o: wrap_mpi.F90
	$(FC) -c $(SPEC_FFLAGS) $<
endif

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<
.F.o:
	$(FC) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)

ESMF_ARCH     := alpha
CFLAGS        := $(cpp_path) -DOSF1 -O2 -omp $(CPPDEF)
FC            := f90
mod_path      := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
FFLAGS        := $(cpp_path) $(mod_path) $(CPPDEF) -r8 -i4 -c -omp -automatic -fpe3 \
		 -check omp_bindings
FFLAGS_DOTF90 := -DOSF1 -free -fpe3
FFLAGS_DOTF   := -extend_source -omp -automatic
LDFLAGS       := -omp -lcxml

ifeq ($(SPMD),TRUE)
  FFLAGS  += -I$(INC_MPI) 
  LDFLAGS += -lmpi
else
  FFLAGS  += -DHIDE_MPI
endif

ifeq ($(SMP),TRUE)
  # THREADED_OMP is used by the timing library
  CFLAGS += -DTHREADED_OMP
endif

ifeq ($(DEBUG),TRUE)
  FFLAGS += -g3 -C
else
# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
# Inline when not debugging
    FORTRAN_OPTIMIZATION := -O2 -inline speed
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
endif

ifeq ($(NO_SWITCH),$(null))
  NO_SWITCH := FALSE
endif
ifneq ($(NO_SWITCH),TRUE)
  LDFLAGS += -lelan
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) $(FFLAGS) $(FFLAGS_DOTF) $<

.F90.o:
	$(FC) $(FFLAGS) $(FFLAGS_DOTF90) $<

.c.o:
	cc -c  $(CFLAGS) $<


endif

#------------------------------------------------------------------------
# NEC SX-6
#------------------------------------------------------------------------

ifeq ($(UNAMES),SUPER-UX)

ESMF_ARCH = SX6

VER = inst
VER = rev285

FC        := sxf90 -Yf,/SX/opt/sxf90/$(VER)/lib -V
CC        := sxc++
AS        := sxas
CPP       := cpp -traditional
PROF      := -ftrace
CPPDEF    += -USX -DFORTRANUNDERSCORE -DNEC_SX -D_SX -DDISABLE_TIMERS
cpp_path  += -I$(ESMF_MOD)/$(ESMF_ARCH)
CFLAGS    := $(PROF)
FFLAGS    := $(PROF) -Wf,"-A idbl4 -pvctl fullmsg loopcnt=1000000 noassume noloopchg" \
             -Wf,"-L fmtlist map transform -ptr byte"
LDFLAGS   := $(PROF) -Wl"-Z 5000M -ZL 5000M"

ifeq ($(DEBUG),TRUE)
	FFLAGS  += -g -Cvsafe
else
# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
    FORTRAN_OPTIMIZATION := -Cvopt
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
endif

ifeq ($(SMP),$(null))
  SMP := TRUE
endif
ifeq ($(SPMD),TRUE)
  cpp_path += -I$(INC_MPI)
  LDFLAGS  += -L$(LIB_MPI) -lmpi
else
  CPPDEF   += -DHIDE_MPI
endif

ifeq ($(SMP),TRUE)
  CPPDEF += -DTHREADED_OMP
  CFLAGS += -P openmp
  FFLAGS += -P openmp
  LDFLAGS += -lcpp -P openmp
endif

CPPFLAGS  := $(CPPDEF) $(cpp_path)
FPPFLAGS  := $(CPPFLAGS)

.SUFFIXES:
.SUFFIXES: .F .F90 .c .s .o
.F.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $<

.c.o:
	$(CC) -c  $(CPPFLAGS) $(CFLAGS) $<
.s.o:
	$(AS) -m  $<

ifeq ($(SPMD),TRUE)
#...jps....added for MPI Global Memory usage
phys_grid.o:        phys_grid.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -gmalloc $<

spmd_dyn.o: spmd_dyn.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -gmalloc $<
endif

radae.o:        radae.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi rexp=phi,psi,fh2oself $<

wv_saturation.o:        wv_saturation.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi rexp=estblf $<

wv_saturation.f90:	wv_saturation.F90
	$(CPP)	$(CPPFLAGS) $< > wv_saturation.f90

cldwat.o:	cldwat.F90 wv_saturation.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi expin=wv_saturation.f90 rexp=estblf $<

FrictionVelocityMod.o:	FrictionVelocityMod.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi $<

shr_orb_mod.f90:        shr_orb_mod.F90
	$(CPP)  $(CPPFLAGS) $< > shr_orb_mod.f90

zenith.o:       zenith.F90 shr_orb_mod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi expin=shr_orb_mod.f90 rexp=shr_orb_cosz $<

# Temporary workaround for error in compiler optimization
courlim.o: courlim.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi -Cvsafe $<

# inlining for CLM

QSatMod.f90:	QSatMod.F90
	$(CPP) $(CPPFLAGS) $< > QSatMod.f90

CanopyFluxesMod.o:	CanopyFluxesMod.F90 QSatMod.f90  FrictionVelocityMod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $< -pi exp=qsat,moninobukini expin=QSatMod.f90,FrictionVelocityMod.f90

Biogeophysics1Mod.o:	Biogeophysics1Mod.F90 QSatMod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $< -pi exp=qsat expin=QSatMod.f90

FrictionVelocityMod.f90	: FrictionVelocityMod.F90
	$(CPP) $(CPPFLAGS) $< > FrictionVelocityMod.f90

BareGroundFluxesMod.o	: BareGroundFluxesMod.F90 FrictionVelocityMod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $< -pi exp=moninobukini expin=FrictionVelocityMod.f90

BiogeophysicsLakeMod.o	: BiogeophysicsLakeMod.F90 FrictionVelocityMod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $< -pi exp=moninobukini expin=FrictionVelocityMod.f90

endif

#------------------------------------------------------------------------
# Earth Simulator
#------------------------------------------------------------------------

ifeq ($(UNAMES),ES)

ESMF_ARCH = ES

FC        := esf90 -V
CC        := esc++
AS        := esas
CPP       := cpp -traditional
PROF      :=
CPPDEF    += -USX -DFORTRANUNDERSCORE -DNEC_SX -D_SX -DDISABLE_TIMERS
cpp_path  += -I$(ESMF_MOD)/$(ESMF_ARCH)
CFLAGS    := $(PROF)
FFLAGS    := $(PROF) -Wf,"-A idbl4 -pvctl fullmsg loopcnt=1000000 noassume noloopchg" \
             -Wf,"-L fmtlist map transform -ptr byte"
LDFLAGS   := $(PROF) -Wl"-Z 5000M -ZL 5000M"

ifeq ($(DEBUG),TRUE)
	FFLAGS  += -g -Cvsafe
else
# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
    FORTRAN_OPTIMIZATION := -Cvopt
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
endif

ifeq ($(SMP),$(null))
  SMP := FALSE
endif
ifeq ($(SPMD),TRUE)
  cpp_path += -I$(INC_MPI)
  LDFLAGS  += -L$(LIB_MPI) -lmpi
else
  CPPDEF   += -DHIDE_MPI
endif

ifeq ($(SMP),TRUE)
  CPPDEF += -DTHREADED_OMP
  CFLAGS += -P openmp
  FFLAGS += -P openmp
  LDFLAGS += -lcpp -P openmp
endif

CPPFLAGS  := $(CPPDEF) $(cpp_path)
FPPFLAGS  := $(CPPFLAGS)

.SUFFIXES:
.SUFFIXES: .F .F90 .c .s .o
.F.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $<

.c.o:
	$(CC) -c  $(CPPFLAGS) $(CFLAGS) $<
.s.o:
	$(AS) -m  $<

ifeq ($(SPMD),TRUE)
#...jps....added for MPI Global Memory usage
phys_grid.o:        phys_grid.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -gmalloc $<

# (mitsuo) add "-gmalloc"
spmd_dyn.o: spmd_dyn.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -gmalloc $<

endif

radae.o:        radae.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi rexp=phi,psi,fh2oself $<

wv_saturation.o:        wv_saturation.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi rexp=estblf $<

wv_saturation.f90:	wv_saturation.F90
	$(CPP)	$(CPPFLAGS) $< > wv_saturation.f90

cldwat.o:	cldwat.F90 wv_saturation.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi expin=wv_saturation.f90 rexp=estblf $<

FrictionVelocityMod.o:	FrictionVelocityMod.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi $<

shr_orb_mod.f90:        shr_orb_mod.F90
	$(CPP)  $(CPPFLAGS) $< > shr_orb_mod.f90

zenith.o:       zenith.F90 shr_orb_mod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi expin=shr_orb_mod.f90 rexp=shr_orb_cosz $<

# Temporary workaround for error in compiler optimization
courlim.o: courlim.F90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) -pi -Cvsafe $<

# inlining for CLM

QSatMod.f90:	QSatMod.F90
	$(CPP) $(CPPFLAGS) $< > QSatMod.f90

CanopyFluxesMod.o:	CanopyFluxesMod.F90 QSatMod.f90  FrictionVelocityMod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $< -pi exp=qsat,moninobukini expin=QSatMod.f90,FrictionVelocityMod.f90

Biogeophysics1Mod.o:	Biogeophysics1Mod.F90 QSatMod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $< -pi exp=qsat expin=QSatMod.f90

FrictionVelocityMod.f90	: FrictionVelocityMod.F90
	$(CPP) $(CPPFLAGS) $< > FrictionVelocityMod.f90

BareGroundFluxesMod.o	: BareGroundFluxesMod.F90 FrictionVelocityMod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $< -pi exp=moninobukini expin=FrictionVelocityMod.f90

BiogeophysicsLakeMod.o	: BiogeophysicsLakeMod.F90 FrictionVelocityMod.f90
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $< -pi exp=moninobukini expin=FrictionVelocityMod.f90

endif

#------------------------------------------------------------------------
# Cray X1
#------------------------------------------------------------------------

ifeq ($(UNAMES),UNICOS/mp)

ifeq ($(SMP),$(null))
  SMP := TRUE
endif

ifeq ($(SPMD),$(null))
  SPMD := TRUE
endif

ifeq ($(SSP),$(null))
  SSP := FALSE
endif

ifeq ($(SSP),FALSE)
  ESMF_ARCH := cray_x1
else
  ESMF_ARCH := cray_x1_ssp
endif

CC := cc
CPPFLAGS := -UCRAY -DUNICOSMP -DCPP_VECTOR $(CPPDEF)
ifeq ($(SMP),FALSE)
  COPT := -h task0
else
  COPT := -h task1
  # THREADED_OMP is used by the timing library
  CPPFLAGS += -DTHREADED_OMP
endif
ifeq ($(SSP),FALSE)
  COPT += -h gen_private_callee
else
  COPT += -h ssp
  CPPFLAGS += -DSSP
endif
CFLAGS := $(cpp_path) $(CPPFLAGS) $(COPT) -h list=m 

FC := ftn
FPPFLAGS := -F -UCRAY -DUNICOSMP -DCPP_VECTOR $(CPPDEF)
FOPT := -Ofp1,nomodinline
# Check for override of default Fortran compiler optimizations
ifneq ($(F_OPTIMIZATION_OVERRIDE),$(null))
  FOPT += $F_OPTIMIZATION_OVERRIDE
endif
ifeq ($(SPMD),FALSE)
  FPPFLAGS += -DHIDE_MPI
endif
ifeq ($(SMP),FALSE)
  FOPT += -Otask0 -xomp
else
  FOPT += -Otask1 -xcsd
endif
ifeq ($(SSP),FALSE)
  FOPT += -Ogen_private_callee
else
  FOPT += -Ossp
  FPPFLAGS += -DSSP
endif
mod_path := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(INC_NETCDF) \
  -p$(ESMF_MOD)/$(ESMF_ARCH) -p$(MOD_NETCDF)
FFLAGS := $(cpp_path) $(mod_path) $(FPPFLAGS) $(FOPT) -s real64 -r a -dp  -dy -rm -em
FREEFLAGS := -f free
FIXEDFLAGS := -f fixed -N 132

LDFLAGS := -L. -p$(ESMF_MOD)/$(ESMF_ARCH) -lsci
ifeq ($(SSP),TRUE)
  LDFLAGS += -Ossp
endif

ifeq ($(DEBUG),TRUE)
  CFLAGS += -Gp
  FFLAGS += -G1
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS)  $<

.c.o:
	$(CC) -c  $(CFLAGS) $<

# Specific targets for chunking reproducibility
vertical_diffusion.o : 
	$(FC) -c $(FREEFLAGS) $(FFLAGS) -Ofp0 $<

# Specific targets for CAM and CLM2 for good performance
wv_saturation.o FrictionVelocityMod.o QSatMod.o shr_orb_mod.o :
	$(FC) -c $(FREEFLAGS) $(FFLAGS) -Omodinline $<
sulchem.o radsw.o radlw.o radae.o :
	$(FC) -c $(FREEFLAGS) $(FFLAGS) -Oscalar3,vector3  $<

endif



#------------------------------------------------------------------------
# Targets/rules that depend on architecture specific variables.
#------------------------------------------------------------------------

# The ESMF library is not made to be built in parallel, so specify only one job is to run.
$(ESMF_LIB)/$(ESMF_ARCH)/libesmf.a:
	cd $(ESMF_ROOT); \
           echo "Build the ESMF library."; \
           echo "ESMF is NOT supported by the CCSM project, but by the ESMF core team in NCAR/SCD"; \
           echo "See http://www.esmf.ucar.edu"; \
           $(MAKE) -j 1 BOPT=$(ESMF_BOPT) ESMF_BUILD=$(ESMF_BLD) ESMF_DIR=$(ESMF_ROOT) ESMF_ARCH=$(ESMF_ARCH)

time_manager.o : $(ESMF_LIB)/$(ESMF_ARCH)/libesmf.a

# Build the CARMA libraries.
$(CARMA_BLD)/aer/libaer.a :
	mkdir $(CARMA_BLD)/aer; \
	cd $(CARMA_BLD); \
	echo "Build the CARMA aerosol library."; \
	echo "CARMA is NOT supported by the CCSM project, but by Ackerman and McKie at NASA Ames."; \
	$(MAKE) -j 1 -f $(CARMA_ROOT)/$(CARMA_MAKE) BOPT=$(CARMA_BOPT) CARMA_ROOT=$(CARMA_ROOT) library_aer

$(CARMA_BLD)/rad/librad.a :
	mkdir $(CARMA_BLD)/rad 2>/dev/null; \
	cd $(CARMA_BLD); \
	echo "Build the CARMA radiation library."; \
	echo "CARMA is NOT supported by the CCSM project, but by Ackerman and McKie at NASA Ames."; \
	$(MAKE) -j 1 -f $(CARMA_ROOT)/$(CARMA_MAKE) BOPT=$(CARMA_BOPT) CARMA_ROOT=$(CARMA_ROOT) library_rad

RM := rm
# Add user defined compiler and load flags if set, and replace FC if USER option set.
CFLAGS  += $(USER_CFLAGS)
FFLAGS  += $(USER_FFLAGS)
LDFLAGS  += $(USER_LDFLAGS)
ifneq ($(USER_FC),$(null))
FC := $(USER_FC)
endif

clean:
	$(RM) -rf carma
	$(RM) -rf esmf
	$(RM) -f Depends Srcfiles *.o *.mod *.stb *.f90 $(MODEL_EXEDIR)/$(EXENAME)

realclean: clean

include $(CURDIR)/Depends
