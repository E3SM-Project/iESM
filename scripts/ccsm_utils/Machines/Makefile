#===============================================================================
# Common Makefile: a framework for building all CESM components and more
#
# Command-line variables
#   MODEL=<model>  ~ a standard macro definition, often found in the included 
#                    MACFILE, used to trigger special compilation flags
# Supported compilers
# ibm, bgl, bgp, pgi, intel, pathscale, gnu, lahey
#===============================================================================

# Set up special characters
null  :=

# Load dependency search path.
dirs := . $(shell cat Filepath)
cpp_dirs := $(dirs)
# Add INCROOT to path for Depends and Include
MINCROOT := 
ifdef INCROOT
  cpp_dirs += $(INCROOT)
  MINCROOT := $(INCROOT)
endif

# Expand any tildes in directory names. Change spaces to colons.
VPATH := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
VPATH := $(subst $(space),:,$(VPATH))               

RM    := rm
CP    := cp

exec_se: $(EXEC_SE)  $(CURDIR)/Depends
complib: $(COMPLIB)  $(CURDIR)/Depends

# Determine whether to compile threaded or not
ifeq ($(strip $(BUILD_THREADED)),TRUE)
   compile_threaded = true
endif
ifeq ($(strip $(SMP)),TRUE)
   compile_threaded = true
endif


ifndef MOD_SUFFIX
   MOD_SUFFIX := mod
endif

#===============================================================================
# set CPP options (must use this before any flags or cflags settings)
#===============================================================================

CPPDEFS := $(USER_CPPDEFS) -D$(OS) 

# USE_ESMF_LIB is currently only defined in env_build.xml
ifeq ($(USE_ESMF_LIB), TRUE)
   CPPDEFS += -DUSE_ESMF_LIB
endif

# ESMF_INTERFACE is currently only defined in env_build.xml
ifeq ($(COMP_INTERFACE), ESMF)
   CPPDEFS += -DESMF_INTERFACE
else
   CPPDEFS += -DMCT_INTERFACE
endif

ifeq ($(strip $(MPILIB)),mpi-serial)
  CPPDEFS += -DNO_MPI2 -DNO_MPIMOD
else
  CPPDEFS += -DHAVE_MPI
endif
ifeq ($(compile_threaded), true)
  CPPDEFS += -DTHREADED_OMP
endif

LD = $(MPIFC)

include $(CASEROOT)/Macros

ifndef AR
   AR := ar
endif

ifdef NETCDF_PATH
  ifndef INC_NETCDF
    INC_NETCDF:=$(NETCDF_PATH)/include
  endif
  ifndef LIB_NETCDF
    LIB_NETCDF:=$(NETCDF_PATH)/lib
  endif
endif
ifdef PNETCDF_PATH
  ifndef $(INC_PNETCDF)
    INC_PNETCDF:=$(PNETCDF_PATH)/include
  endif
  ifndef LIB_PNETCDF
    LIB_PNETCDF:=$(PNETCDF_PATH)/lib
  endif
endif


# Set HAVE_SLASHPROC on LINUX systems which are not bluegene or Darwin (OSx)

ifeq ($(findstring -DLINUX,$(CPPDEFS)),-DLINUX)
  ifneq ($(findstring DBG,$(CPPDEFS)),DBG)
    ifneq ($(findstring Darwin,$(CPPDEFS)),Darwin)
      CPPDEFS += -DHAVE_SLASHPROC
    endif
  endif
endif

ifdef CPRE
  FPPDEFS := $(patsubst -D%,$(CPRE)%,$(CPPDEFS)) 
else
  FPPDEFS := $(CPPDEFS)
endif


#===============================================================================
# Set config args for pio and mct to blank and then enable serial 
#===============================================================================
ifndef CONFIG_ARGS
  CONFIG_ARGS := 
endif
ifeq ($(MPILIB),mpi-serial)
   CONFIG_ARGS+= --enable-mpiserial
endif
ifeq ($(MODEL),pio)
  CONFIG_ARGS+= --enable-timing
  ifeq ($DEBUG,TRUE)
     CONFIG_ARGS+= --enable-debug
  endif
endif

#===============================================================================
# User-specified INCLDIR
#===============================================================================

INCLDIR := -I. 
ifdef USER_INCLDIR
  INCLDIR += $(USER_INCLDIR)
endif

#===============================================================================
# MPI-serial library (part of MCT)
#===============================================================================

ifeq ($(strip $(MPILIB)), mpi-serial)
  CC      := $(SCC)
  FC      := $(SFC)
  MPIFC   := $(SFC)
  MPICC   := $(SCC)
  CONFIG_ARGS += --enable-mpiserial --disable-netcdf4 MCT_PATH=$(EXEROOT)/mct/mpi-serial
  INC_MPI := $(EXEROOT)/mct/mpi-serial
else
  CC := $(MPICC)
  FC := $(MPIFC)
  ifdef MPI_PATH
    INC_MPI := $(MPI_PATH)/include
    LIB_MPI := $(MPI_PATH)/lib
  endif
endif

#===============================================================================
# Set include paths (needed after override for any model specific builds below)
#===============================================================================

ifdef INC_NETCDF
  INCLDIR += -I$(INC_NETCDF)
endif
ifdef MOD_NETCDF
  INCLDIR += -I$(MOD_NETCDF)
endif
ifdef INC_MPI
  INCLDIR += -I$(INC_MPI)
endif 
ifdef INC_PNETCDF
  INCLDIR += -I$(INC_PNETCDF)
endif

ifeq ($(MODEL),driver)
  INCLDIR += -I$(EXEROOT)/atm/obj -I$(EXEROOT)/lnd/obj -I$(EXEROOT)/ice/obj -I$(EXEROOT)/ocn/obj -I$(EXEROOT)/glc/obj -I$(EXEROOT)/rof/obj -I$(EXEROOT)/iac/obj
endif

ifndef MCT_LIBDIR
  MCT_LIBDIR=$(EXEROOT)/mct
endif

ifndef PIO_LIBDIR
  PIO_LIBDIR=$(EXEROOT)/pio
endif

INCLDIR +=  -I$(MCT_LIBDIR)/mct -I$(MCT_LIBDIR)/mpeu -I$(PIO_LIBDIR) -I$(CODEROOT)/utils/pio -I$(EXEROOT)/lib/include

#
# Use the MCT dir for the cache for all configure calls because it is the first one 
#
CFLAGS+=$(CPPDEFS)

CONFIG_ARGS +=  CC="$(SCC)" FC="$(SFC)" MPICC="$(MPICC)" \
                MPIFC="$(MPIFC)" FCFLAGS="$(FFLAGS) $(FREEFLAGS) $(INCLDIR)" \
                CPPDEFS="$(CPPDEFS)" CFLAGS="$(CFLAGS) -I.. $(INCLDIR)" \
                NETCDF_PATH=$(NETCDF_PATH) LDFLAGS="$(SLIBS) $(LDFLAGS)"
ifdef PNETCDF_PATH
  CONFIG_ARGS+= PNETCDF_PATH=$(PNETCDF_PATH)
endif

ifdef PIO_CONFIG_OPTS
  CONFIG_ARGS+= $(PIO_CONFIG_OPTS)
endif 

FFLAGS+=$(FPPDEFS)

ifeq ($(findstring -cosp,$(CAM_CONFIG_OPTS)),-cosp)
# The following is for the COSP simulator code:
COSP_LIBDIR:=$(EXEROOT)/atm/obj/cosp
endif

ifeq ($(MODEL),cam)
   # These RRTMG files take an extraordinarily long time to compile with optimization.
   # Until mods are made to read the data from files, just remove optimization from
   # their compilation.
rrtmg_lw_k_g.o: rrtmg_lw_k_g.f90
	$(FC) -c $(FPPFLAGS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<
rrtmg_sw_k_g.o: rrtmg_sw_k_g.f90
	$(FC) -c $(FPPFLAGS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<


ifeq ($(strip $(COMPILER)),lahey)
binary_io.o: binary_io.F90
	$(FC) -c $(FPPDEFS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<
wrap_nf.o: wrap_nf.F90
	$(FC) -c $(FPPDEFS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<
wrap_mpi.o: wrap_mpi.F90
	$(FC) -c $(FPPDEFS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<
apex_subs.o: apex_subs.F90
	$(FC) -c $(FPPDEFS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<
endif
ifdef COSP_LIBDIR
INCLDIR+=-I$(COSP_LIBDIR) -I$(COSP_LIBDIR)/../
$(COSP_LIBDIR)/libcosp.a: abortutils.o
	$(MAKE) -C $(COSP_LIBDIR) F90='$(FC)' F90FLAGS='$(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS) $(FC_AUTO_R8)' \
	F90FLAGS_noauto='$(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS)' \
	F90FLAGS_fixed='$(INCLDIR) $(INCS) $(FIXEDFLAGS) $(FFLAGS) $(FC_AUTO_R8)'

cospsimulator_intr.o: $(COSP_LIBDIR)/libcosp.a
endif

endif

# The following rules are for the COSP simulator code:


# Set esmf.mk location with ESMF_LIBDIR having precedent over ESMFMKFILE
CCSM_ESMFMKFILE := undefined_CCSM_ESMFMKFILE
ifdef ESMFMKFILE
   CCSM_ESMFMKFILE := $(ESMFMKFILE)
endif
ifdef ESMF_LIBDIR
   CCSM_ESMFMKFILE := $(ESMF_LIBDIR)/esmf.mk
endif


# System libraries (netcdf, mpi, pnetcdf, esmf) 
ifndef SLIBS
   SLIBS := -L$(LIB_NETCDF) -lnetcdf 
endif
ifdef LIB_PNETCDF
   SLIBS += -L$(LIB_PNETCDF) -lpnetcdf
endif
ifdef LAPACK_LIBDIR
   SLIBS += -L$(LAPACK_LIBDIR) -llapack -lblas
endif
ifdef LIB_MPI
   ifndef MPI_LIB_NAME
      SLIBS += -L$(LIB_MPI) -lmpi
   else  
      SLIBS += -L$(LIB_MPI) -l$(MPI_LIB_NAME)
   endif
endif

# For compiling and linking with external ESMF.
# If linking to external ESMF library then include esmf.mk 
# ESMF_F90COMPILEPATHS
# ESMF_F90LINKPATHS
# ESMF_F90LINKRPATHS
# ESMF_F90ESMFLINKLIBS
ifeq ($(USE_ESMF_LIB), TRUE)
  include $(CCSM_ESMFMKFILE)
  FFLAGS += $(ESMF_F90COMPILEPATHS)
  SLIBS  += $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

# User libraries (mct, mpeu, pio)
ifndef ULIBS
   ifdef LIBROOT
     ULIBS += -L$(LIBROOT) -lcsm_share 
   endif
endif

# Component libraries - hard-wired for current ccsm components
ifndef CLIBS
   ifdef LIBROOT
     CLIBS  += -L$(LIBROOT) -latm -llnd -lice -locn -lglc -lrof -liac
   endif
   ifeq ($(OCN_SUBMODEL),moby)
     CLIBS  += -lmoby
   endif
endif

# Machine stuff to appear last on the link step
ifndef MLIBS
     MLIBS  :=
endif

#------------------------------------------------------------------------------
# Drive configure scripts for support libraries (mct and pio)
#------------------------------------------------------------------------------

$(MCT_LIBDIR)/Makefile.conf: 
	cd $(MCT_LIBDIR); \
	$(CONFIG_SHELL) $(CODEROOT)/utils/mct/configure $(CONFIG_ARGS) --srcdir $(CODEROOT)/utils/mct

$(MCT_LIBDIR)/mpeu/libmpeu.a: $(MCT_LIBDIR)/Makefile.conf
	$(MAKE) -C $(MCT_LIBDIR)/mpeu

$(MCT_LIBDIR)/mct/libmct.a: $(MCT_LIBDIR)/mpeu/libmpeu.a
	$(MAKE) -C $(MCT_LIBDIR)/mct

$(MCT_LIBDIR)/mpi-serial/libmpi-serial.a: $(MCT_LIBDIR)/Makefile.conf
	$(MAKE) -C $(MCT_LIBDIR)/mpi-serial


MCTLIBS = $(MCT_LIBDIR)/mct/libmct.a $(MCT_LIBDIR)/mpeu/libmpeu.a 
PIOLIB = $(PIO_LIBDIR)/libpio.a

ULIBS += -L$(MCT_LIBDIR)/mct -lmct -L$(MCT_LIBDIR)/mpeu -lmpeu -L$(PIO_LIBDIR) -lpio -lgptl
ifdef IAC_LIBS
  ULIBS += $(IAC_LIBS)
endif


$(PIO_LIBDIR)/Makefile.conf: 
	cd $(PIO_LIBDIR); \
	$(CONFIG_SHELL) $(CODEROOT)/utils/pio/configure $(CONFIG_ARGS) --srcdir $(CODEROOT)/utils/pio

$(PIO_LIBDIR)/libpio.a: $(PIO_LIBDIR)/Makefile.conf
	$(MAKE) -C $(PIO_LIBDIR)


#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------

touch_filepath: 
	touch $(CURDIR)/Filepath

# Get list of files and build dependency file for all .o files
#   using perl scripts mkSrcfiles and mkDepends

SOURCES := $(shell cat Srcfiles)
OBJS    := $(addsuffix .o, $(basename $(SOURCES)))
INCS    := $(foreach dir,$(cpp_dirs),-I$(dir)) 

CURDIR := $(shell pwd)

$(CURDIR)/Depends: $(CURDIR)/Srcfiles $(CURDIR)/Deppath
	$(CASETOOLS)/mkDepends Deppath Srcfiles > $@

$(CURDIR)/Deppath: $(CURDIR)/Filepath
	$(CP) -f $(CURDIR)/Filepath $@
	@echo "$(MINCROOT)" >> $@

$(CURDIR)/Srcfiles: $(CURDIR)/Filepath
	$(CASETOOLS)/mkSrcfiles 

$(CURDIR)/Filepath:
	@echo "$(VPATH)" > $@


#-------------------------------------------------------------------------------
# echo file names, paths, compile flags, etc. used during build
#-------------------------------------------------------------------------------

db_files:
	@echo " "
	@echo "* MACFILE := $(MACFILE)"
	@echo "* VPATH   := $(VPATH)"
	@echo "* INCS    := $(INCS)"
	@echo "* OBJS    := $(OBJS)"
db_flags:
	@echo " "
	@echo "* cc      := $(CC)  $(CFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F.o    := $(FC)  $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F90.o  := $(FC)  $(FFLAGS) $(FREEFLAGS) $(INCS) $(INCLDIR)"

#-------------------------------------------------------------------------------
# Rules used for the tests run by "configure -test"
#-------------------------------------------------------------------------------

test_fc: test_fc.o
	$(LD) -o $@ test_fc.o $(LDFLAGS)
test_nc: test_nc.o
	$(LD) -o $@ test_nc.o -L$(LIB_NETCDF) -lnetcdf $(LDFLAGS)
test_mpi: test_mpi.o
	$(LD) -o $@ test_mpi.o $(LDFLAGS)
test_esmf: test_esmf.o
	$(LD) -o $@ test_esmf.o $(LDFLAGS)

#-------------------------------------------------------------------------------
# build rules: 
#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .F .f90 .c .o

ifeq ($(ULIBDEP),$(null))
   ifneq ($(LIBROOT),$(null))
     ULIBDEP += $(LIBROOT)/libcsm_share.a
     ULIBDEP += $(LIBROOT)/libatm.a
     ULIBDEP += $(LIBROOT)/libice.a
     ULIBDEP += $(LIBROOT)/liblnd.a
     ULIBDEP += $(LIBROOT)/libocn.a
     ULIBDEP += $(LIBROOT)/libglc.a
     ULIBDEP += $(LIBROOT)/librof.a
     ULIBDEP += $(LIBROOT)/libiac.a
   endif
endif
ifdef COSP_LIBDIR
  CLIBS += -L$(COSP_LIBDIR) -lcosp
endif



ifeq ($(MPILIB),mpi-serial)
  MPISERIAL = $(MCT_LIBDIR)/mpi-serial/libmpi-serial.a
  MLIBS+= -L$(MCT_LIBDIR)/mpi-serial -lmpi-serial
endif


$(MCTLIBS) $(PIOLIB) : $(MPISERIAL)

$(OBJS): $(MCTLIBS) $(PIOLIB) $(SOURCES)

$(EXEC_SE): $(OBJS) $(ULIBDEP)
	$(LD) -o $(EXEC_SE) $(OBJS) $(CLIBS) $(ULIBS) $(SLIBS) $(MLIBS) $(LDFLAGS)

$(COMPLIB): $(OBJS)
	$(AR) -r $(COMPLIB) $(OBJS)

.c.o:
	$(CC) -c $(INCLDIR) $(INCS) $(CFLAGS)  $<
.F.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FIXEDFLAGS) $<
.f90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $<
.F90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $<

mostlyclean:
	$(RM) -f *.f *.f90 

clean: mostlyclean
	$(RM) -f *.d *.$(MOD_SUFFIX) $(OBJS) Srcfiles Filepath Deppath Depends

realclean: clean
	$(RM) -f $(EXEC_SE)

# the if-tests prevent DEPS files from being created when they're not needed
ifneq ($(MAKECMDGOALS), db_files)
ifneq ($(MAKECMDGOALS), db_flags)
ifneq ($(MAKECMDGOALS), mostlyclean)
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), realclean)
    -include $(CURDIR)/Depends $(CASEROOT)/Depends.$(COMPILER) $(CASEROOT)/Depends.$(MACH)
endif
endif
endif
endif
endif
