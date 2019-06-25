#! /bin/csh
#
#
#
#
#
#=======================================================================
#
# make-cam.csh
#
# Builds CAM for a variety of platforms. After checking the environment
# variables, the settings will be displayed and the user prompted to
# continue. If the build already exists the user will be prompted to
# build from scratch or to rebuild. If the build completes succesfully
# a namelist file already exists, the user will be prompted with whether
# it should be rebuilt. Providing a -q on the command line will prevent
# the prompts from being generated.
#
# The build assumes the following:
#   - the libraries are found in $CAM_BASE/Libraries
#   - the source files are in $CAM_BASE/Models/cam/$CAM_BRANCH/source
#   - the data files are in $CAM_BASE/Models/cam/$CAM_BRANCH/data
#   - the build goes to $CAM_BASE/Models/cam/$CAM_BRANCH/build
#   - the run goes to $CAM_BASE/Models/cam/$CAM_BRANCH/run/$CAM_RUN
#  
# Environment Variable
#   The following environmnt variables can be set before running the script
#   to override the default behavior:
#
#   CAM_BASE
#   (varies) The path to the cam directory.
#
#   CAM_BRANCH
#   (rad) The name of the subdirectory that is the branch containing the
#   source, data, and scripts.
#
#   CAM_CONFIGURE_OPTS
#   (varies) Optional flags to for the configure command. There are a couple of
#   special cases: '-pergro' and '-offline_test'. '-pergro' sets
#   CAM_CONFIGURE_OPTS and CAM_NAMELIST to the values required for the
#   perturbation growth test. '-offline_test' sets CAM_CONFIGURE_OPTS and
#   CAM_NAMELIST to the values required for CAM offline test.
#
#   CAM_DEBUG
#   (OFF) Controls whether a debug build is made [OFF, ON].
#
#   CAM_DYN
#   (fv) The model's dynamics package [eul, fv, sld].
#
#   CAM_BASE
#   (varies) The path to the directory containing the netcdf (and mpich)
#   libraries.
#
#   CAM_MODEL
#   (varies) The parallelism model used to run the code [CPU, HYBRID, MPI, OMP].
#
#   CAM_NAMELIST
#   (varies) Values to be added to the namelist
#   generated with the build-namelist utility.
#
#   CAM_OFFLINE
#   (OFF) Controls whether dynamics are generated internally or read in
#   from a meteorology file [OFF, NCEP ECMWF].
#
#   CAM_RES
#   (2x2.5) The model's resolution.
#
#-----------------------------------------------------------------------
# Usage: 
#   make-cam.csh [-q]
#
# NOTE: Since CARMA isn't thread safe, OMP and HYBRID can not be used.
#=======================================================================

# Many people have these aliased to add the -i which can cause problems with
# the scripts.
unalias cp
unalias mv
unalias rm

# Check the OS and MACHINE types.
setenv OS `uname -s`
setenv MACHINE_TYPE `uname -m`
setenv MACHINE `hostname`

   setenv CAM_F_COMPILER ifort
   setenv CAM_C_COMPILER icc
   setenv USER_FC $CAM_F_COMPILER
   setenv USER_CC $CAM_C_COMPILER
   setenv CAM_CPPDEFS '-D TITAN_FAO -D _DIML -D _TITAN'
## -D EQTORQ'

#HDF4/5 paths
setenv LIB_HDF4 $H4DIR/lib
setenv LIB_HDF5 $H5DIR/lib

#set MPI paths
setenv LIB_MPI $OMPIDIR/lib
setenv INC_MPI $OMPIDIR/include
#setenv INC_MPI /curc/tools/x_86_64/rh6/openmpi/1.6.4/intel/13.0.0/torque/4.1.4/include/
#setenv LIB_MPI /curc/tools/x_86_64/rh6/openmpi/1.6.4/intel/13.0.0/torque/4.1.4/lib
#setenv INC_MPI /curc/tools/free/redhat_5_x86_64/mpich2-1.4.1p1_ict-3.2.2.013_ib/include
#setenv INC_MPI /curc/tools/free/redhat_5_x86_64/mvapich2-1.6_ict-3.2.2.013/include
#setenv INC_MPI /curc/tools/free/redhat_5_x86_64/openmpi-1.4.3_ict-3.2.2.013_torque-2.5.8_ib/include
#setenv INC_MPI /home/larsonej/Libraries/mpich/include
echo "${0}:  Set INC_MPI to $INC_MPI"
echo "${0}:  Set LIB_MPI to $LIB_MPI"

#set mpirun = /home/larsonej/Libraries/mpich/bin/mpirun

#set NETCDF paths
setenv INC_NETCDF $NETCDF/include
setenv LIB_NETCDF $NETCDF/lib
#setenv INC_NETCDF /curc/tools/x_86_64/rh6/parallel-netcdf4/4.3/hdf5/1.8.11/hdf4/4.2.9/szip/2.1/zlib/1.2.78/jpeglib/8d/openmpi/1.6.4/intel/13.0.0/torque/4.1.4/ib/include/
#setenv LIB_NETCDF /curc/tools/x_86_64/rh6/parallel-netcdf4/4.3/hdf5/1.8.11/hdf4/4.2.9/szip/2.1/zlib/1.2.78/jpeglib/8d/openmpi/1.6.4/intel/13.0.0/torque/4.1.4/ib/lib
#setenv INC_NETCDF /curc/tools/x_86_64/rh5/netcdf/4.1.3/hdf5/1.8.9/openmpi/1.6/intel/12.1.4/torque/2.5.11/netcdf-4.1.3_hdf5-1.8.9_openmpi-1.6_intel-12.1.4_torque-2.5.11/include/
#setenv LIB_NETCDF /curc/tools/x_86_64/rh5/netcdf/4.1.3/hdf5/1.8.9/openmpi/1.6/intel/12.1.4/torque/2.5.11/netcdf-4.1.3_hdf5-1.8.9_openmpi-1.6_intel-12.1.4_torque-2.5.11/lib/
#setenv INC_NETCDF /curc/tools/x_86_64/rh6/netcdf4/4.3/hdf5/1.8.11/hdf4/4.2.9/szip/2.1/zlib/1.2.78/jpeglib/8d/intel/13.0.0/include/
#setenv LIB_NETCDF /curc/tools/x_86_64/rh6/netcdf4/4.3/hdf5/1.8.11/hdf4/4.2.9/szip/2.1/zlib/1.2.78/jpeglib/8d/intel/13.0.0/lib/
echo "${0}:  Set INC_NETCDF to $INC_NETCDF"
setenv MOD_NETCDF $INC_NETCDF
echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF"
echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF"
   
  
   #setenv CAM_BASE /projects/$LOGNAME/software/titancam
   setenv CAM_BASE /home/$LOGNAME/titancam 
   setenv CSMDATA     /projects/$LOGNAME/data    
   set rundir       = /lustre/janus_scratch/$LOGNAME/run
   setenv CAM_MODEL MPI
   setenv CAM_RUN default

switch ( $CAM_MODEL )
  case CPU:
    set use_smp 	= -nosmp
    set use_spmd 	= -nospmd  
    breaksw;
  case HYBRID:
    set use_smp 	= -smp
    set use_spmd 	= -spmd  
    breaksw;
  case MPI:
    set use_smp 	= -nosmp
    set use_spmd 	= -spmd  
    breaksw;
  case OMP:
    set use_smp 	= -smp
    set use_spmd 	= -nospmd  
    breaksw;
  default:
    echo "${0}:  ERROR: Unknown CAM_MODEL option, '$CAM_MODEL'"
    exit 1
    breaksw;
endsw
    
# Determine the src and data directories.
  setenv CAM_BRANCH rad_dfz

# Use finite volume at 10x15 as the default dynamics package.
  setenv CAM_DYN fv
  setenv CAM_RES 10x15

# Determine if we are debugging.
  setenv CAM_DEBUG ON

if ( $CAM_DEBUG == "ON" ) then
  setenv debugsw -debug
else
  setenv debugsw
endif

# Determine if we are using offline dynamics.
if ( ! $?CAM_OFFLINE ) then
  setenv CAM_OFFLINE OFF
endif

# ROOT OF CAM DISTRIBUTION - probably needs to be customized.
# Contains the source code for the CAM distribution.
# (the root directory contains the subdirectory "models")
if ( ! $?CAMROOT ) then
  set CAMROOT      = $CAM_BASE/$CAM_BRANCH/src
endif

# ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
# Contains the initial and boundary data for the CAM distribution.
# (the root directory contains the subdirectories "atm" and "lnd")
if ( ! $?CSMDATA ) then
  setenv CSMDATA     $CAM_BASE/titandata
endif

# $wrkdir is a working directory where the model will be built and run.
# $blddir is the directory where model will be compiled.
# $rundir is the directory where the model will be run.
# $cfgdir is the directory containing the CAM configuration scripts.
#set wrkdir       = $CAM_BASE
set wrkdir       = /projects/larsonej/titancam
set blddir       = $wrkdir/build_dfz-test
echo "${0}: Building in $blddir"
set cfgdir       = $CAMROOT/models/atm/cam/bld
#set rundir       = $wrkdir/run/CAM_RUN
set dyndir       = $wrkdir/run/dyn

if ( ! $?CAM_CONFIGURE_OPTS ) then
  setenv CAM_CONFIGURE_OPTS '-nadv 47'
  setenv CAM_CONFIGURE_OPTS "$CAM_CONFIGURE_OPTS -nnadv 9"
  setenv CAM_CONFIGURE_OPTS "$CAM_CONFIGURE_OPTS -nlev 61"
  setenv CAM_CONFIGURE_OPTS "$CAM_CONFIGURE_OPTS -pcols 16"
endif


# Default namelist settings:
# $runtype is the run type: initial, restart, or branch.
# $nelapse is the number of timesteps to integrate, or number of days if negative.
set runtype      = initial
set nelapse      = -10

if ( ! $?CAM_NAMELIST ) then
  setenv CAM_NAMELIST "&camexp nelapse=$nelapse rest_pfile='$wrkdir/run' mss_irt=0 carma_flag=.true. carma_do_coag=.true. carma_do_emission=.true. carma_do_vtran=.true. carma_do_wetdep=.true. / &clmexp rpntpath='$wrkdir/run' /"
endif

if ( $CAM_OFFLINE == "ECMWF" ) then
  set dyndir       = $wrkdir/dyn
  setenv CAM_CPPDEFS "-DOFFLINE_DYN"
  setenv CAM_NAMELIST "&camexp nelapse=$nelapse rest_pfile='$wrkdir/run' mss_irt=0 carma_flag=.true. carma_do_coag=.true. carma_do_emission=.true. carma_do_vtran=.true. carma_do_wetdep=.true. met_data_file='$dyndir/met0001.nc' ncdata='$dyndir/ic.nc' / &clmexp rpntpath='$wrkdir/run'/"
  setenv CAM_CONFIGURE_OPTS "$CAM_CONFIGURE_OPTS -nlev 61"
else if ( $CAM_OFFLINE == "NCEP" ) then
  set dyndir       = $wrkdir/dyn
  setenv CAM_CPPDEFS "-DOFFLINE_DYN"
  setenv CAM_NAMELIST "&camexp nelapse=$nelapse rest_pfile='$wrkdir/run' mss_irt=0 carma_flag=.true. carma_do_coag=.true. carma_do_emission=.true. carma_do_vtran=.true. carma_do_wetdep=.true. met_data_file='$dyndir/met0001.nc' ncdata='$dyndir/ic.nc' / &clmexp rpntpath='$wrkdir/run'/"
  setenv CAM_CONFIGURE_OPTS "$CAM_CONFIGURE_OPTS -nlev 61"
endif
echo "${0}:  Set offline to $CAM_OFFLINE"

echo "${0}:  Using configure options: " $CAM_CONFIGURE_OPTS

# Make sure that the user likes the setting of the environment variables
# before doing anything destructive.
if ( $#argv == 0 ) then 
  echo "${0}:  Continue? (Y/n)"
  set confirm = $<
  if ( $confirm == 'n' ) then
    echo "${0}:  Quiting ..."
    exit 0
  endif
endif

# Do our best to get sufficient stack memory
limit stacksize unlimited

# Ensure that run and build directories exist
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1


# Build CAM
cd $blddir                  || echo "cd $blddir failed" && exit 1

# Make sure that the user likes the setting of the environment variables
# before doing anything destructive.
set confirm = 'Y'

if ( -f $blddir/config_cache.xml ) then
  set confirm = 'n'
  if ( $#argv == 0 ) then 
    echo "${0}: Configure again? (complete rebuild) (Y/n)"
    set confirm = $<
  endif
endif

#set compilers (from richards make file)
set compilers = "-fc $CAM_F_COMPILER"
if ( $confirm != 'n' ) then
  if ( $?CAM_CPPDEFS ) then
    $cfgdir/configure $use_spmd $use_smp -dyn $CAM_DYN -res $CAM_RES $CAM_CONFIGURE_OPTS $debugsw -cppdefs "$CAM_CPPDEFS"  || echo "${0}:  configure failed" && exit 1
  else
    $cfgdir/configure $use_spmd $use_smp -dyn $CAM_DYN -res $CAM_RES $CAM_CONFIGURE_OPTS $debugsw -cppdefs "$CAM_CPPDEFS" || echo "${0}:  configure failed" && exit 1
  endif
endif 

echo "${0}:  building CAM in $blddir ..."
rm -f Depends
echo "${0}:  Started at `date`";
gmake >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
echo "${0}:  Completed at `date`";
