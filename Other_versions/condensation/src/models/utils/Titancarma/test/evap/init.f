      subroutine init
c
c
c  @(#) init.f  McKie  Oct-1995
c  This routine performs all initializations at the beginning
c  of each run of the model.
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables
c
      include 'globals.h'
c
c
c  Define formats
c
    1 format('Initialization for ',a,' (Version ',a,')')
    2 format(a,':  ',i6)
    3 format(a,':  ',f12.2)
    4 format(a,':  ',a)
    5 format(/,'Model will run with the following values:')
    6 format(a,':  ',L7)
c
c
c---------------------------------------------------------------------------
c
c  Define run control values that can change from run to run.
c   (Could be input from a data file at this spot in the
c    code, but why not just explicitly define them here)
c
c
c  Define begin & end timestep indices for this run.
c   ibtime is 0 for a cold start new simulation.
c   ibtime is the ending timestep of previous run for a restart.
c
      ibtime = 0
      ietime = 60
c
c
c  Define names of input & output files for this run
c
      resifil = 'model_res.in'		! Restart input file
      iniifil = 'model_ini.in'          ! Cold start initialization file
c     prtofil = '/dev/tty'
      prtofil = 'model.p'               ! Output print file
      hisofil = 'model_his.out'         ! Output history file (binary)
      resofil = 'model_res.out'         ! Restart output file
c
c
c  Define # timestep cycles for print & history output
c
      nprint = 10          ! frequency of output to print file
      nhist = 5000         ! frequency of output to history file
      nrest = 5000         ! frequency of output to restart file
c
c
c  Define flags for various processes:
c
c
c  Define flag to control whether any coagulation is to be simulated
c
      do_coag = .false.
c
c
c  Define flag to control whether condensational growth is to be simulated
c  (evaporation and nucleation also)
c
      do_grow = .true.
c
c
c  Define flag to control if error trapping for debugging is to be done
c
      do_error = .true.
c
c
c  End of per run control values definition (usually no changes below here)
c
c---------------------------------------------------------------------------
c
c
c  Open output print file
c
      open(unit=LUNOPRT,file=prtofil,status='unknown')
c
c
c  Open output history file
c
      open(unit=LUNOHIS,file=hisofil,status='unknown',
     $     form='unformatted')
c
c
c  Open output restart file
c
      open(unit=LUNORES,file=resofil,status='unknown',
     $     form='unformatted')
c
c
c  Announce entry to this routine 
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter init'
c
c
c  Report model name & version tag
c
      write(LUNOPRT,1) PROGNAM, PROGTAG
c
c
c  Set up error trapping (for debugging) if it was requested
c
      if( do_error )then
       call setuperr
      endif
c
c
c  Initialize # history timepoints output in this run
c
      khist = 0
c
c
c  Do either:
c    A cold start to begin a new simulation, or
c    A restart from a previous simulation
c
      if( ibtime .eq. 0 )then
       call initnew
      else
       call initres
      endif
c
c
c  Report some initialization values
c
      write(LUNOPRT,5)
      write(LUNOPRT,2) 'ibtime', ibtime
      write(LUNOPRT,2) 'ietime', ietime
      write(LUNOPRT,2) 'nprint', nprint
      write(LUNOPRT,2) 'nhist', nhist
      write(LUNOPRT,2) 'nrest', nrest
      write(LUNOPRT,2) 'nx', nx
      write(LUNOPRT,2) 'ny', ny
      write(LUNOPRT,2) 'nz', nz
      write(LUNOPRT,3) 'time', time
      write(LUNOPRT,2) 'itime', itime
      write(LUNOPRT,3) 'dtime', dtime
      write(LUNOPRT,6) 'do_error', do_error
      write(LUNOPRT,6) 'do_coag', do_coag
      write(LUNOPRT,6) 'do_grow', do_grow
      call dblank(simtitle, ns)
      write(LUNOPRT,4) 'simtitle', simtitle(1:ns)
c
c
c  Write initial state to print file and history file
c
      call outprt
      call outhis
c
c
c  Return to caller with model initializations complete
c
      return
      end
