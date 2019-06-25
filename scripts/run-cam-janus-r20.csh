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
#PBS -q janus-long
# Maximum number of processes
#PBS -l nodes=4:ppn=12
#The following line makes it so I am not sharing nodes and crashing other people's code
##PBS -l naccesspolicy=singlejob
# Walltiem for the job, format D:HH:MM:SS
#PBS -l walltime=5:00:00:00
# output file base name
#PBS -N titancam
# Put standard error and standard out in same file
#PBS -j oe
# Pass on environemnt variables
#PBS -V
# Email notification
#PBS -M larsonej@colorado.edu
#PBS -m abe
# End of options
#=======================================================================

# Many people have these aliased to add the -i which can cause problems with
# the scripts.
unalias cp
unalias mv
unalias rm
# rest_pfile		= '/Volumes/Data/Models/cam_carma/titancarma/run/apr1/cam.pointer'


setenv runname p10m10r20-mon
setenv CAM_RUN $runname

# Check the OS and MACHINE types.
setenv OS `uname -s`
setenv MACHINE_TYPE `uname -m`
setenv MACHINE `hostname`
setenv NODENAME `echo $MACHINE | grep ^node`
set mpirun = /curc/tools/free/redhat_5_x86_64/openmpi-1.4.3_ict-3.2.2.013_torque-2.5.8_ib/bin/mpirun
set mpiboot = /curc/tools/free/redhat_5_x86_64/openmpi-1.4.3_ict-3.2.2.013_torque-2.5.8_ib/bin/mpiboot
set mpdallexit = /curc/tools/free/redhat_5_x86_64/openmpi-1.4.3_ict-3.2.2.013_torque-2.5.8_ib/bin/mpdallexit

#set mpirun = /curc/tools/free/redhat_5_x86_64/mpich2-1.4.1p1_ict-3.2.2.013_ib/bin/mpirun
#set MPDBOOT = /curc/tools/nonfree/redhat_5_x86_64/mpich2-1.4.1p1_ict-3.2.2.013_ib/bin/mpdboot
#set MPDALLEXIT = /curc/tools/nonfree/redhat_5_x86_64/mpich2-1.4.1p1_ict-3.2.2.013_ib/bin/mpdallexit

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

#      if ( ! $?CAM_BASE ) then
        setenv CAM_BASE /home/$LOGNAME/titancam
#      endif

      if ( ! $?CAM_MODEL ) then
        setenv CAM_MODEL MPI
      endif
        setenv CAM_BASE /home/$LOGNAME/titancam
      endif

#  setenv CAM_BRANCH rad

# Determine the run id.
#if ( ! $?CAM_RUN ) then
#  setenv CAM_RUN defaultrun
#endif
echo "${0}:  Set run to $CAM_RUN"

# Determine if we are debugging.
if ( ! $?CAM_DEBUG ) then
  setenv CAM_DEBUG OFF
endif
echo "${0}:  Set debug to $CAM_DEBUG";

# Determine the number of threads.
#if ( ! $?CAM_THREADS ) then
#  setenv CAM_THREADS 12
#endif
#echo "${0}:  Set threads to $CAM_THREADS";

# Do our best to get sufficient stack memory
limit stacksize unlimited
limit datasize unlimited

setenv CAM_BASE /home/$LOGNAME/titancam
echo $CAM_BASE

# $blddir is the directory where model will be compiled.
# $rundir is the directory where the model will be run.
set blddir       = $CAM_BASE/build_$runname
set rundir       = /lustre/janus_scratch/$LOGNAME/run/$CAM_RUN

# Ensure that run and build directories exist
mkdir -p $rundir                || echo "${0}:  cannot create $rundir" && exit 1

# Determine the run type.
if ( ! $?CAM_RUNTYPE ) then

  # If the run directory already exists and has the cam executable,
  # then consider it a restart. Otherwise it is a initial run.
    setenv CAM_RUNTYPE INITIAL
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
#  if ( ! -f namelist ) then
#    echo "${0}:  namelist not found, execute make-cam.csh first ..." 
#    exit 1
#  else
#    echo "${0}:  copying namelist from $blddir ..." 
    cp $CAM_BASE/lists/namelist_$runname $rundir/namelist
#  endif
  
else

  # Use the executable from the run directory
  if ( ! -x $rundir/cam ) then
    echo "${0}:  $rundir/cam not found, do an INITIAL run first ..."
    exit 1
  else
    echo "${0}:  using cam from $rundir ..."   
  endif

  # Use the executable from the run directory
  if ( ! -f $rundir/namelist_$runname ) then
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

    if ( $?PBS_JOBID ) then
      echo "${0}: Running on nodes ..."
      echo `cat $PBS_NODEFILE`
#      set nproc = `cat $PBS_NODEFILE | wc -w`
      set nproc = `wc -l < $PBS_NODEFILE`
      set NNODE = `uniq $PBS_NODEFILE | wc -l`

      if ( $CAM_MODEL == "MPI" ) then
        echo "${0}: NNODES, nproc, PBS_NODEFILE    $NNODE  $nproc  $PBS_NODEFILE "

        setenv MPD_CON_EXT ${PBS_JOBID}

#        $mpiboot -n $NNODE -f $PBS_NODEFILE -v --remcons
        $mpirun -np $nproc $rundir/cam < $rundir/namelist >&! $rundir/$outfile || echo "${0}:  CAM run failed" && exit 1
#        $mpirun -machinefile $PBS_NODEFILE -np $nproc cam < namelist || echo "${0}:  CAM run failed" && exit 1
        $mpdallexit
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

echo "${0}:  Completed at `date`";

exit 0
