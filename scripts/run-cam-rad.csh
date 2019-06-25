#'! /bin/tcsh 
#
#=======================================================================
#
#  run-mac.csh
#
#  Run CAM on a variety of machines.
#
# The build assumes the following:
#   - the data files are in $CAM_BASE/cam/$CAM_BRANCH/data
#   - the build is in $CAM_BASE/cam/$CAM_BRANCH/build
#   - the run goes to $CAM_BASE/cam/$CAM_BRANCH/run/$CAM_RUN
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
#   CAM_DEBUG
#   (OFF) Controls whether a debug build is made [OFF, ON].
#
#   CAM_MODEL
#   (varies) The parallelism model used to run the code [CPU, HYBRID, MPI, OMP].
#
#   CAM_THREADS
#   (2) The number of threads (nodes) to use when running in OMP (MPI) mode.
#
#   CAM_RUN
#   (varies) The subdirectory under the run directory from which the
#   model will be run.
#
#   CAM_RUNTYPE
#   (varies) The type of run [INITIAL, RESTART, BRANCH].
#
#
# NOTE: Since CARMA isn't thread safe, OMP and HYBRID can not be used.
#
#  The section below is the batch script for PBS machines.
#-----------------------------------------------------------------------
# Name of the queue
#PBS -q routeq
# Maximum number of processes
#PBS -l nodes=4:ppn=2:compute
# output file base name
#PBS -N cam
# Put standard error and standard out in same file
#PBS -j oe
# Pass on environemnt variables
#PBS -V
# Email notification
##PBS -M <your e-mail address>
##PBS -m abe
# End of options
#=======================================================================

# Many people have these aliased to add the -i which can cause problems with
# the scripts.
unalias cp
unalias mv
unalias rm
# rest_pfile		= '/Volumes/Data/Models/cam_carma/titancarma/run/apr1/cam.pointer'

# Check the OS and MACHINE types.
setenv OS `uname -s`
setenv MACHINE_TYPE `uname -m`
setenv MACHINE `hostname`
setenv NODENAME `echo $MACHINE | grep ^node`

switch ( $OS )
  case Darwin:

    # Determine the parallelism model. The default is OMP.
    if ( ! $?CAM_MODEL ) then
      setenv CAM_MODEL MPI
    endif
      
    # Determine the base directory.    
    if ( ! $?CAM_LIBRARY ) then
      setenv CAM_LIBRARY /Volumes/Data/Libraries
    endif
      
    set mpirun = $CAM_LIBRARY/mpich/bin/mpiexec
    echo "${0}:  Set mpirun to $mpirun"

    # Determine the base directory.    
    if ( ! $?CAM_BASE ) then
      setenv CAM_BASE /Volumes/Data/Models/cam_carma/titancarma
      echo $CAM_BASE
    endif
    breaksw;

  case Linux:
    # Choose a compiler
    if ( ! $?CAM_F_COMPILER ) then
    
      # Default to PathScale on Opterons and Portland Group elsewhere.
      if ( $MACHINE_TYPE == "x86_64") then
        setenv CAM_F_COMPILER pathf90
      else
        setenv CAM_F_COMPILER pgf90
      endif
    endif
    
    # Chose the location for the libraries.
    if ( $CAM_F_COMPILER == "pathf90" ) then
      echo "${0}:  Using $CAM_F_COMPILER compiler with pathcc"
      
      if ( $NODENAME != "" ) then
        set mpirun = /home/bardeen/libs/mpich-1.2.6/bin/mpirun
        echo "${0}:  Set mpirun to $mpirun"
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
    
      if ( $NODENAME != "" ) then
        set mpirun = /usr/local/mpich/1.2.6/pgi/x86_64/bin/mpirun
        echo "${0}:  Set mpirun to $mpirun"
      endif 
    endif

    # Determine the base directory and parallelism models.    
    if ( $NODENAME != "" ) then
    
      # We need to use PBS.
      if ( ! $?PBS_JOBID ) then
        echo "${0}: This batch script must be submitted via PBS."
      else 
        echo "${0}:  Running CAM on Linux using PBS."
        echo "${0}:  Set mpirun to $mpirun"
        echo "${0}:  Running on `hostname`"
        echo "${0}:  Working directory is $PBS_O_WORKDIR"
      endif

      if ( ! $?CAM_BASE ) then
        setenv CAM_BASE /home/$LOGNAME/cam
      endif

      if ( ! $?CAM_MODEL ) then
        setenv CAM_MODEL MPI
      endif
    else if ( $MACHINE == "cynewulf.lasp.colorado.edu" ) then
      echo "${0}: This batch script must be submitted via PBS."
      exit 1
    else if ( $MACHINE == "cadfael.lasp.colorado.edu" ) then
      echo "${0}: This batch script must be submitted via PBS."
      exit 1
    else
      if ( ! $?CAM_BASE ) then
        setenv CAM_BASE /home/$LOGNAME/cam
      endif

      if ( ! $?CAM_MODEL ) then
        setenv CAM_MODEL CPU
      endif
    endif
    breaksw;

  default:
    echo "${0}:  ERROR - Unknown OS type $OS"
    exit 1
    breaksw;
endsw

# Determine the branch name.
if ( ! $?CAM_BRANCH ) then
  setenv CAM_BRANCH rad
endif
echo "${0}:  Set branch to $CAM_BRANCH"

# Determine the run id.
if ( ! $?CAM_RUN ) then
  setenv CAM_RUN 090604_2
endif
echo "${0}:  Set run to $CAM_RUN"

# Determine if we are debugging.
if ( ! $?CAM_DEBUG ) then
  setenv CAM_DEBUG OFF
endif
echo "${0}:  Set debug to $CAM_DEBUG";

# Determine the number of threads.
if ( ! $?CAM_THREADS ) then
  setenv CAM_THREADS 8
endif
echo "${0}:  Set threads to $CAM_THREADS";

# Do our best to get sufficient stack memory
limit stacksize unlimited
limit datasize unlimited

# $blddir is the directory where model will be compiled.
# $rundir is the directory where the model will be run.
set blddir       = $CAM_BASE/$CAM_BRANCH/build
set rundir       = $CAM_BASE/$CAM_BRANCH/run/$CAM_RUN

# Ensure that run and build directories exist
mkdir -p $rundir                || echo "${0}:  cannot create $rundir" && exit 1

# Determine the run type.
if ( ! $?CAM_RUNTYPE ) then

  # If the run directory already exists and has the cam executable,
  # then consider it a restart. Otherwise it is a initial run.
    setenv CAM_RUNTYPE BRANCH
endif
echo "${0}:  Set runtype to $CAM_RUNTYPE"

switch ( $CAM_RUNTYPE )
  case INITIAL:
    set runtype 	= 0
    breaksw;
  case RESTART:
    set runtype 	= 1
    breaksw;
  case BRANCH:
    set runtype 	= 3
    breaksw;
  default:
    echo "${0}:  ERROR: Unknown CAM_RUNTYPE option, '$CAM_RUNTYPE'"
    exit 1
    breaksw;
endsw

echo "${0}:  Set model to $CAM_MODEL"

# Copy the CAM executable.
if ( $CAM_RUNTYPE == INITIAL ) then
  
  # Copy the executable from the build directory.
  if ( ! -x $blddir/cam ) then
    echo "${0}:  $blddir/cam not found, execute make-cam.csh first ..."
    exit 1
  else
    echo "${0}:  copying cam from $blddir ..." 
    cp $blddir/cam $rundir/cam
  endif

  # Copy the namelist file to the run directory.
  if ( ! -f namelist ) then
    echo "${0}:  namelist not found, execute make-cam.csh first ..." 
    exit 1
  else
    echo "${0}:  copying namelist from $blddir ..." 
    cp namelist $rundir/namelist
  endif
  
else

  # Use the executable from the run directory
  if ( ! -x $rundir/cam ) then
    echo "${0}:  $rundir/cam not found, do an INITIAL run first ..."
    exit 1
  else
    echo "${0}:  using cam from $rundir ..."   
  endif

  # Use the executable from the run directory
  if ( ! -f $rundir/namelist ) then
    echo "${0}:  $rundir/namelist not found, do an INITIAL run first ..."
    exit 1
  else
    echo "${0}:  using namelist from $rundir ..."  
  endif
endif

# Set the directory where the restart file goes to be the run directory.
# This is controled by the variables RESP_PFILE and RPNTPATH in the namelist.
#
# NOTE: Be careful if you edit this file. The characters in the brackets are
# a space and a tab. The tab Can NOT be converted to spaces and still have
# scripts work properly. Some editors try to do this automatically for you,
# so be careful.
sed "s,^[ 	]*rest_pfile[ 	]*=[ 	]*'[A-Za-z0-9/_ ]*', rest_pfile		= '$rundir/cam.rpointer'," < $rundir/namelist > $rundir/tmp.namelist
mv -f $rundir/tmp.namelist $rundir/namelist

sed "s,^[ 	]*rpntpath[ 	]*=[ 	]*'[A-Za-z0-9/_ ]*', rpntpath	 	= '$rundir/lnd.rpointer'," < $rundir/namelist > $rundir/tmp.namelist
mv -f $rundir/tmp.namelist $rundir/namelist

# Set the runtype to be the value in runtype.
sed "s,^[ 	]*nsrest[ 	]*=[ 	]*[0-9], nsrest		= $runtype," < $rundir/namelist > $rundir/tmp.namelist
mv -f $rundir/tmp.namelist $rundir/namelist

# Cleanup the any lines that don't start with an & or a space. These are
# probably errors from namelist where the text was too long for a line
# and was wrapped around.
sed -e :a -e '/,$/N; s/,\n[ ]*/, /; ta' < $rundir/namelist > $rundir/tmp.namelist
mv -f $rundir/tmp.namelist $rundir/namelist

# Run CAM
cd $rundir                      || echo "${0}:  cd $rundir failed" && exit 1

# Edit the namelist file so that the restart file goes
echo "${0}:  Running CAM in $rundir using $CAM_MODEL";
echo "${0}:  Started at `date`";

setenv outfile RUN`date "+%y%m%d%H%M%S"`.txt

switch ( $OS )
  case Darwin:
    if ( $CAM_MODEL == "CPU" ) then
      if ( $CAM_DEBUG == "ON" ) then
      
        # NOTE: Once in gdb do: run < namelist.
        gdb cam
      else
#        cam < namelist !>& $outfile || echo "CAM run failed" && exit 1
        cam < namelist || echo "CAM run failed" && exit 1
      endif
    else if ( $CAM_MODEL == "OMP" ) then
      if ( $MACHINE_TYPE == "i386" ) then 
        setenv OMP_NUM_THREADS $CAM_THREADS
        setenv KMP_STACKSIZE 64M
       else
        setenv XLSMPOPTS "stack=64000000:parthds=2"
      endif
#      cam < namelist >&! $outfile || echo "CAM run failed" && exit 1
      cam < namelist || echo "cam run failed" && exit 1 
    else if ( $CAM_MODEL == "MPI" ) then
#EJL
        $mpirun -n $CAM_THREADS cam < namelist !>& $outfile || echo "${0}:  CAM run failed" && exit 1
#        $mpirun -n $CAM_THREADS cam < namelist || echo "${0}: CAM run failed" && exit 1
#        $mpirun -n $CAM_THREADS cam < namelist || echo "${0}: CAM run failed" && exit 1 &
    else
      echo "${0}: Unsupported parallelism model $CAM_MODEL."
      exit 1
    endif
    breaksw;
  
  case Linux:
    if ( $?PBS_JOBID ) then
      echo "${0}: Running on nodes ..."
      echo `cat $PBS_NODEFILE`
      set nproc = `cat $PBS_NODEFILE | wc -w`

      if ( $CAM_MODEL == "MPI" ) then
        $mpirun -machinefile $PBS_NODEFILE -np $nproc cam < namelist || echo "${0}:  CAM run failed" && exit 1
      else if ($CAM_MODEL == "HYBRID" ) then
        $mpirun -machinefile $PBS_NODEFILE -np $nproc cam < namelist || echo "${0}:  CAM run failed" && exit 1
      else
        echo "${0}: Unsupported parallelism model $CAM_MODEL."
        exit 1
      endif
    else
      if ( $CAM_MODEL == "CPU" ) then
        if ( $CAM_DEBUG == "ON" ) then
          
          # NOTE: Once in gdb do: run < namelist.
          gdb cam
        else
          cam < namelist !>& $outfile || echo "${0}:  CAM run failed" && exit 1
        endif
      else if ( $CAM_MODEL == "OMP" ) then
        setenv OMP_NUM_THREADS $CAM_THREADS
        setenv MPSTKZ 256M
        cam < namelist >&! $outfile || echo "${0}:  CAM run failed" && exit 1
      else
        echo "${0}: Unsupported parallelism model $CAM_MODEL."
        exit 1
      endif
    endif
    breaksw;
endsw

echo "${0}:  Completed at `date`";

exit 0
