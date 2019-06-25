#! /bin/tcsh
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

setenv CAM_F_COMPILER mpif90
setenv CAM_C_COMPILER mpicc
setenv USER_FC $CAM_F_COMPILER    
setenv USER_CC $CAM_C_COMPILER
setenv CAM_BASE /home/$LOGNAME/
setenv CSMDATA /home/$LOGNAME/titandata


switch ( $OS )
  case Darwin:
    # Choose a compiler
    if ( ! $?CAM_F_COMPILER ) then
    
      # Default to PathScale on Opterons and Portland Group elsewhere.
      if ( $MACHINE_TYPE == "i386") then
        setenv CAM_F_COMPILER ifort
      else
        setenv CAM_F_COMPILER xlf90_r
      endif
    endif

    # Use gcc, since there was link trouble when using icc.
    if (! $?CAM_C_COMPILER ) then
      setenv CAM_C_COMPILER gcc
    endif

    echo "${0}:  Using $CAM_F_COMPILER compiler with CAM_C_COMPILER"
    
    # Determine the base directory.    
    if ( ! $?CAM_LIBRARY ) then
      setenv CAM_LIBRARY /Volumes/Data/Libraries
    endif
      
    # Chose the location for the libraries.
    setenv INC_MPI $CAM_LIBRARY/mpich/include
    echo "${0}:  Set INC_MPI to $INC_MPI"
    setenv LIB_MPI $CAM_LIBRARY/mpich/lib
    echo "${0}:  Set LIB_MPI to $LIB_MPI"
    set mpirun = $CAM_LIBRARY/mpich/bin/mpirun
    echo "${0}:  Set mpirun to $mpirun"
    
    setenv INC_NETCDF $CAM_LIBRARY/netcdf/include
    echo "${0}:  Set INC_NETCDF to $INC_NETCDF"
    setenv MOD_NETCDF $INC_NETCDF
    echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF"
    setenv LIB_NETCDF $CAM_LIBRARY/netcdf/lib
    echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF"
    
    # Determine the base directory.    
    if ( ! $?CAM_BASE ) then
      setenv CAM_BASE /Volumes/Data/Models/cam_carma/titancarma
      echo $CAM_BASE
    endif

    # Determine the parallelism model. The default is OMP.
    if ( ! $?CAM_MODEL ) then
      setenv CAM_MODEL MPI
    endif
      
    setenv USER_FC $CAM_F_COMPILER
    setenv USER_CC gcc

    breaksw;

  case Linux:
    # Choose a compiler
    if ( ! $?CAM_F_COMPILER ) then
    
      # Default to PathScale on Opterons and Portland Group elsewhere.
      if ( $MACHINE_TYPE == "x86_64") then
        setenv CAM_F_COMPILER pathf90
      else
        setenv CAM_F_COMPILER pgf90

        # LASP doesn't have a pgcc license, so switch to gcc
        if (! $?CAM_C_COMPILER ) then
          setenv CAM_C_COMPILER gcc
        endif
      endif
    endif
    
    # Chose the location for the libraries.
    if ( $CAM_F_COMPILER == "pathf90" ) then
      echo "${0}:  Using $CAM_F_COMPILER compiler with pathcc"
      
      if ( $MACHINE == "cynewulf.lasp.colorado.edu" ) then
        setenv INC_MPI /home/bardeen/libs/mpich-1.2.6/include
        echo "${0}:  Set INC_MPI to $INC_MPI"
        setenv LIB_MPI /home/bardeen/libs/mpich-1.2.6/lib
        echo "${0}:  Set LIB_MPI to $LIB_MPI"
        set mpirun = /home/bardeen/libs/mpich-1.2.6/bin/mpirun
        echo "${0}:  Set mpirun to $mpirun"
        setenv INC_NETCDF /home/bardeen/libs/netcdf-3.6.0-p1/pathscale/include
        echo "${0}:  Set INC_NETCDF to $INC_NETCDF"
        setenv MOD_NETCDF $INC_NETCDF
        echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF"
        setenv LIB_NETCDF /home/bardeen/libs/netcdf-3.6.0-p1/pathscale/lib
        echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF"
      else if ( $MACHINE == "cadfael.lasp.colorado.edu" ) then
        setenv INC_MPI /home/bardeen/libs/mpich-1.2.6/include
        echo "${0}:  Set INC_MPI to $INC_MPI"
        setenv LIB_MPI /home/bardeen/libs/mpich-1.2.6/lib
        echo "${0}:  Set LIB_MPI to $LIB_MPI"
        set mpirun = /home/bardeen/libs/mpich-1.2.6/bin/mpirun
        echo "${0}:  Set mpirun to $mpirun"
        setenv INC_NETCDF /contrib/2.6/netcdf/3.6.0-p1-pathscale-2.4-64/include
        echo "${0}:  Set INC_NETCDF to $INC_NETCDF"
        setenv MOD_NETCDF $INC_NETCDF
        echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF"
        setenv LIB_NETCDF /contrib/2.6/netcdf/3.6.0-p1-pathscale-2.4-64/lib
        echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF"
      else 
        echo "${0}: No library paths configured for this machine."
        exit 1
      endif
    
    else if ( $CAM_F_COMPILER == "pgf90" ) then

      if ( $?CAM_C_COMPILER ) then
        if ( $CAM_C_COMPILER == "gcc" ) then
          echo "${0}:  Using pgf90 compiler with gcc"
          setenv USER_CC gcc
          echo "${0}:  Set USER_CC to $USER_CC"
        else 
          echo "${0}:  Using pgf90 comiler with pgcc"
        endif
      else
        echo "${0}:  Using pgf90 compiler with pgcc"
      endif
    
      if ( $MACHINE == "cynewulf.lasp.colorado.edu" ) then
        setenv INC_MPI /usr/local/mpich/1.2.6/pgi/x86_64/include
        echo "${0}:  Set INC_MPI to $INC_MPI"
        setenv LIB_MPI /usr/local/mpich/1.2.6/pgi/x86_64/lib
        echo "${0}:  Set LIB_MPI to $LIB_MPI"
        set mpirun = /usr/local/mpich/1.2.6/pgi/x86_64/bin/mpirun
        echo "${0}:  Set mpirun to $mpirun"
        setenv INC_NETCDF /usr/local/netcdf/include
        echo "${0}:  Set INC_NETCDF to $INC_NETCDF"
        setenv MOD_NETCDF $INC_NETCDF
        echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF"
        setenv LIB_NETCDF /usr/local/netcdf/lib
        echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF"
      else 
        setenv INC_NETCDF /tools/netcdf/netcdf-3.6.1/include
        setenv MOD_NETCDF /tools/netcdf/netcdf-3.6.1/include
        setenv LIB_NETCDF /tools/netcdf/netcdf-3.6.1/lib
        echo "${0}:  Set INC_NETCDF to $INC_NETCDF"
        echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF"
        echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF"
      endif
    endif

    # Determine the base directory and parallelism models.    
    if (( $MACHINE == "cynewulf.lasp.colorado.edu" ) || \
        ( $MACHINE == "cadfael.lasp.colorado.edu" )) then
    
      if ( ! $?CAM_BASE ) then
        setenv CAM_BASE /home/$LOGNAME/cam
      endif

      if ( ! $?CAM_MODEL ) then
        setenv CAM_MODEL MPI
      endif
    else
      
      if ( ! $?CAM_BASE ) then
        setenv CAM_BASE /home/$LOGNAME/cam
      endif

      if ( ! $?CAM_MODEL ) then
        setenv CAM_MODEL CPU
      endif

    endif

    setenv USER_FC $CAM_F_COMPILER
    breaksw;

  default:
    echo "${0}:  ERROR - Unknown OS type $OS"
    exit 1
    breaksw;
endsw

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
    
echo "${0}:  Set model to $CAM_MODEL"

# Determine the src and data directories.
if ( ! $?CAM_BRANCH ) then
  setenv CAM_BRANCH titancam
  echo "Branch: '$CAM_BRANCH'"
endif
echo "${0}:  Set branch to $CAM_BRANCH"

# Use finite volume at 10x15 as the default dynamics package.
if ( ! $?CAM_DYN ) then
  setenv CAM_DYN fv
endif
echo "${0}:  Using " $CAM_DYN " dynamics"

if ( ! $?CAM_RES ) then
  setenv CAM_RES 10x15
endif
echo "${0}:  Using " $CAM_RES " resolution"

# Determine if we are debugging.
if ( ! $?CAM_DEBUG ) then
  setenv CAM_DEBUG OFF
endif

if ( $CAM_DEBUG == "ON" ) then
  setenv debugsw -debug
else
  setenv debugsw
endif
echo "${0}:  Set debug to $CAM_DEBUG"

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
set wrkdir       = $CAM_BASE/$CAM_BRANCH
set blddir       = $wrkdir/build
set cfgdir       = $CAMROOT/models/atm/cam/bld
set rundir       = $wrkdir/run/CAM_RUN
set dyndir       = $wrkdir/run/dyn

if ( ! $?CAM_CONFIGURE_OPTS ) then
  setenv CAM_CONFIGURE_OPTS '-nadv 47'
  setenv CAM_CONFIGURE_OPTS "$CAM_CONFIGURE_OPTS -nnadv 9"
  setenv CAM_CONFIGURE_OPTS "$CAM_CONFIGURE_OPTS -nlev 61"
  setenv CAM_CONFIGURE_OPTS "$CAM_CONFIGURE_OPTS -pcols 4"
endif


# Default namelist settings:
# $runtype is the run type: initial, restart, or branch.
# $nelapse is the number of timesteps to integrate, or number of days if negative.
set runtype      = initial
set nelapse      = -1

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

if ( $confirm != 'n' ) then
  if ( $?CAM_CPPDEFS ) then
    $cfgdir/configure $use_spmd $use_smp -dyn $CAM_DYN -res $CAM_RES $CAM_CONFIGURE_OPTS $debugsw -cppdefs "$CAM_CPPDEFS"  || echo "${0}:  configure failed" && exit 1
  else
    $cfgdir/configure $use_spmd $use_smp -dyn $CAM_DYN -res $CAM_RES $CAM_CONFIGURE_OPTS $debugsw || echo "${0}:  configure failed" && exit 1
  endif
endif 

echo "${0}:  building CAM in $blddir ..."
rm -f Depends
echo "${0}:  Started at `date`";
gmake >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
echo "${0}:  Completed at `date`";

## Create the namelist
cd $blddir                      || echo "${0}:  cd $blddir failed" && exit 1
set confirm = 'Y'

if ( -f $blddir/namelist ) then
  set confirm = 'n'
  if ( $#argv == 0 ) then 
    echo "${0}: Recreate namelist? (Y/n)"
    set confirm = $<
  endif
endif

if ( $confirm != 'n' ) then
  $cfgdir/build-namelist -s -case $CAM_BRANCH -runtype $runtype -o $blddir/namelist -namelist "$CAM_NAMELIST"  || echo "${0}:  build-namelist failed" && exit 1
endif
