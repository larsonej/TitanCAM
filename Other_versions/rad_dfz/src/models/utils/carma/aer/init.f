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
      include 'globaer.h'
c
c
c  Declare local variables
c
      logical all_ok
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
    7 format('Error--(init) ',a,'=',i5,' not max of ',
     $  a,'=',i5,3x,a,'=',i5)
    8 format('Warning--(init): ',a,' because do_parcel = .true.')
    9 format(/,'End of model initialization')
c
c---------------------------------------------------------------------------
c
c
c  Define run control values that can change from run to run.
c   (Could be input from a data file at this spot in the
c    code, but why not just explicitly define them here)
c
c NOTE: A number of flags are now being set in via CAM namelist
c variables. The deault values for these are set in carma.F90 in
c the physics package. Do not change the values here. The values
c can be changed by setting variables in the CAM namelist file.
c
c NOTE: The variables that haven't been commented out in this file
c should probably not be changed as most of them required to have
c the current value to work with CAM.
c
c  Define begin & end timestep indices for this run.
c   <ibtime> is 0 for a cold start new simulation.
c   <ibtime> is the ending timestep of previous run for a restart.
c   <ietime> is the maximum ending timestep for current run
c  Note:  <time> .gt. <endtime> also ends current run, so
c         <ietime> and <endtime> control when current runs end.
c
      ibtime = 0
      ietime = 800000
c
c
c  Total simulation time for this run
c
      endtime = 1800.
c
c
c  Define timestep size [s].
c
c      dtime = 10.d0
c
c
c  Define flag to control if history output is done to netcdf format file
c    .true.  for netcdf history output
c    .false. for traditional Fortran binary output
c
      do_netcdf = .false.
c
c
c  Define names of input & output files for this run.
c
c      resifil = 'carma_res.in'		! Restart input file

c     prtofil = '/dev/tty'              ! Output print file
c      prtofil = 'carma.p'               ! Output print file

c      if( do_netcdf )then
c       hisofil = 'carma_his.cdf'        ! Output history file (netcdf)
c      else
c       hisofil = 'carma_his.bin'        ! Output history file (binary)
c      endif

c      resofil = 'carma_res.out'         ! Restart output file
c     resofil = '/dev/null'             ! Restart output file

c      stepofil = 'substep.out'          ! Timestepping output file
c     stepofil = '/dev/null'            ! Timestepping output file

c      radofil = 'carma_rad.out'         ! Radiation submodel print output
c
c
c  Define frequencies of print and history output:
c   use timestep period (nprint, nhist, nrest) when > 0, otherwise
c   use time period (pprint, phist, prest).
c
c
      nprint = -1         ! timestep period between outputs to print file
      nhist  = -1         ! timestep period between outputs to history file
      nrest  = -20        ! timestep period between outputs to restart file

      pprint = 86400.     ! time period between outputs to print file [s]
      phist  = 100000.    ! time period between outputs to history file [s]
      prest  = 100000.    ! time period between outputs to restart file [s]
c
c
c  Define frequency for radiation calcs (same convention as above)
c
      nrad   = -1         ! timestep period between radiation calcs
      prad   =  60.       ! time period between radiation calcs [s]
c
c
c  Define flags for various processes:
c
c
c  Define flags for overall timestepping output to
c   print, history, & restart files
c
c      do_print = .true.
      do_hist = .false.
      do_rest = .false.
c
c
c  Define flag to control whether radiative transfer is to be computed
c
c      do_rad = .false.
c
c
c  Define flag to control whether layers below the model domain are included in
c   the radiative transfer model
c
      do_below = .false.
c
c
c  Define flag to control whether the model is to be run as a parcel simulation
c   (multiple parcels may be simulated by allowing NZ > 1).
c
      do_parcel = .false.
c
c
c  Define flag to control whether any coagulation is to be simulated.
c
      do_coag = .true.
c
c
c  Define flag to control whether condensational growth is to be simulated
c   (evaporation and nucleation also).
c
      do_grow = .false.
c
c
c  Define flag to control whether temperature is changed by latent heating.
c   Note: Setting this flag to .false. will not prevent 
c   potential temperature <pt> or potential temperature concentration <ptc>
c   from changing due to transport.
c
      do_thermo = .false.
c
c
c  Define flag to control whether sedimentation is substepped
c   (in which case it is not done through vertical transport)
c
      do_implised = .false.
c
c
c  Define flag to control whether vertical transport occurs.
c
      do_vtran = .true.
c
c
c  Define flags to control vertical boundary conditions:
c    <itbnd_pc> = I_FIXED_CONC: use specified concentration at the boundary
c               = I_FLUX_SPEC:  use specified flux
c    <ibbnd_pc>: same as <itbnd_pc>, but for bottom boundary;
c   equivalent parameters are defined for gases and potential temperature
c   boundary conditions.
c
      itbnd_pc  = I_FLUX_SPEC
      ibbnd_pc  = I_FLUX_SPEC
      itbnd_gc  = I_FLUX_SPEC
      ibbnd_gc  = I_FLUX_SPEC
      itbnd_ptc = I_FIXED_CONC
      ibbnd_ptc = I_FIXED_CONC
c
c
c  Define flag to control whether horizontal transport occurs in the
c   east-west or north-south directions (always .false. for a parcel simulation).
c
      do_ew = .false.
      do_ns = .false.
c
c
c  Define flag to control which horizontal advection algorithm is
c   used:
c    <ihoradv> = I_PPM: Use piecewise polynomial method
c    <ihoradv> = I_GALERKIN: Use Galerkin method with Chapeau functions
c
      ihoradv = I_PPM
c
c
c  Define the minimum and maximum number of time substeps for fast
c   microphysics (nucleation and condensation), as well as the threshold
c   particle concentration [cm^-3] in a grid cell, below which the
c   minimum number of substeps is always used.
c
c      minsubsteps = 1
c      maxsubsteps = 1
c      conmax = 1.e-1
c
c
c  Define flag to control whether a variable time-step should be used.
c   Also define min and max time-steps, maximum tolerance in particle
c   concentration changes <dpctol>, maximum tolerance in gas concentration 
c   changes <dgstol>, and minimum relative concentration to consider <conmax>
c   (i.e., for bins with concentrations less than <conmax>*max(pc), we don't
c   worry about how large the changes were).
c
      do_varstep = .false.
c
c
c  Set up things that depend on variable timestepping
c
      if( do_varstep )then
        dtmin  = 2.e-3
        dtmax  = 5.e0
        dpctol = 0.8
        dgstol = 0.2
      else
        do_step = .true.
      endif
c
c
c  Define flag to control if error trapping for debugging is to be done
c   (May have no effect on some systems.  Mainly useful for sunos.)
c
c      do_error = .false.
c
c
c  End of per run control values definition (usually no changes below here)
c
c---------------------------------------------------------------------------
c
c
c  Open output print file
c
c  NOTE: Is now stdout, so doesn't need to be opened.
c
c      open(unit=LUNOPRT,file=prtofil,status='unknown')
c
c
c  Open output history file if traditional binary output is requested (non-netcdf)
c
c      if( do_hist )then
c       if( .not. do_netcdf )then
c        open(unit=LUNOHIS,file=hisofil,status='unknown',
c     $       form='unformatted')
c       endif
c      endif
c
c
c  Open output restart file
c
c      if( do_hist )then
c       open(unit=LUNORES,file=resofil,status='unknown',
c     $      form='unformatted')
c      endif
c
c
c  Open output file for timestep diagnostics
c
c     open(unit=LUNOSTEP,file=stepofil,status='unknown')
c
c
c  Open file for radiation submodel print output
c
c     open(unit=LUNORAD,file=radofil,status='unknown')
c
c
c  Announce entry to this routine 
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter init'
c
c
c  Report model name & version tag
c
c      if (do_print_setup) call prtsep
c      if (do_print_setup) write(LUNOPRT,1) PROGNAM, PROGTAG
c
c
c  Check critical symbolic constants for consistency
c 
      all_ok = .true.
      if( NXORNY .ne. max(NX,NY) )then
       write(LUNOPRT,7) 'NXORNY',NXORNY, 'NX',NX, 'NY',NY
       all_ok = .false.
      endif
      if( NXORNYP1 .ne. max(NXP1,NYP1) )then
       write(LUNOPRT,7) 'NXORNYP1',NXORNYP1, 'NXP1',NXP1, 'NYP1',NYP1
       all_ok = .false.
      endif
      if( .not. all_ok ) call endcarma
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
c      if( ibtime .eq. 0 )then
c       call initnew
c      else
c       call initres
c      endif
c
c
c  Ensure consistency of control flags
c 
      if( do_parcel )then

        if( NXY .ne. 1 )then
          write(LUNOPRT,'(/,a)')
     $      'do_parcel = .true. requires NXY = 1'
          call endcarma
        endif
          
        if( do_vtran )then
          do_vtran = .false.
          write(LUNOPRT,8) 'do_vtran set to .false.'
        endif

        if( do_ns )then
          do_ns = .false.
          write(LUNOPRT,8) 'do_ns set to .false.'
        endif

        if( do_ew )then
          do_ew = .false.
          write(LUNOPRT,8) 'do_ew set to .false.'
        endif

      endif
c
c
c  Check to make sure model includes at least 5 layers if <do_vtran>
c   is .true., <NX> is at least 5 if <do_ew> is .true., and <NY>
c   is at least 5 if <do_ns> is .true.
c
      if( do_vtran .and. NZ.lt.5 )then
        write(LUNOPRT,'(/,a)')
     $    'Cannot do vertical transport with NZ < 5'
        call endcarma
      endif

      if( do_ew .and. NX.lt.5 )then
        write(LUNOPRT,'(/,a)')
     $    'Cannot do east-west transport with NX < 5'
        call endcarma
      endif

      if( do_ns .and. NY.lt.5 )then
        write(LUNOPRT,'(/,a)')
     $    'Cannot do east-west transport with NY < 5'
        call endcarma
      endif
c
c
c  Report some initialization values
c
c      if (do_print_setup) then
c        write(LUNOPRT,5)
c
c        write(LUNOPRT,2) 'ibtime', ibtime
c        write(LUNOPRT,2) 'ietime', ietime
c        write(LUNOPRT,3) 'endtime', endtime
c
c        write(LUNOPRT,2) 'NX', NX
c        write(LUNOPRT,2) 'NY', NY
c        write(LUNOPRT,2) 'NZ', NZ
c
c        write(LUNOPRT,3) 'time', time
c        write(LUNOPRT,2) 'itime', itime
c        write(LUNOPRT,3) 'dtime', dtime
c
c        write(LUNOPRT,6) 'do_error', do_error
c        write(LUNOPRT,6) 'do_netcdf', do_netcdf
c        write(LUNOPRT,6) 'do_parcel', do_parcel
c        write(LUNOPRT,6) 'do_coag', do_coag
c        write(LUNOPRT,6) 'do_grow', do_grow
c        write(LUNOPRT,6) 'do_thermo', do_thermo
c        write(LUNOPRT,6) 'do_vtran', do_vtran
c        write(LUNOPRT,6) 'do_ew', do_ew
c        write(LUNOPRT,6) 'do_ns', do_ns
c
c        write(LUNOPRT,6) 'do_rad', do_rad
c        if( do_rad )then
c          if( nrad .gt. 0 )then
c            write(LUNOPRT,2) 'nrad', nrad
c          else
c            write(LUNOPRT,3) 'prad', prad
c          endif
c        endif
c
c        write(LUNOPRT,6) 'do_print', do_print
c        if( do_print )then
c          if( nprint .gt. 0 )then
c            write(LUNOPRT,2) 'nprint', nprint
c          else
c            write(LUNOPRT,3) 'pprint', pprint
c          endif
c        endif
c
c
c  Possibly write initial state to print file and history file
c
c      if( do_print )then
c        call outprt
c      endif
c
c      if( do_hist )then
c        call outhis
c      endif
c
c
c  Report end of initialization
c
c        write(LUNOPRT,9)
c        call prtsep
c      endif
c
c
c  Return to caller with model initializations complete
c
      return
      end
