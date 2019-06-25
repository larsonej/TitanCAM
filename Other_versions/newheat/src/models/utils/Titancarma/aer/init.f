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
    2 format(a,':  ',i10)
    3 format(a,':  ',f12.2)
    4 format(a,':  ',a)
    5 format(/,'Model will run with the following values:')
    6 format(a,':  ',L7)
    7 format('Error--(init) ',a,'=',i5,' not max of ',
     $  a,'=',i5,3x,a,'=',i5)
    8 format('Warning--(init): ',a,' because do_parcel = .true.')
    9 format(/,'End of model initialization')
   11 format(3(1pe13.6,3x))
   12 format(e6.2, e13.4)
c
c---------------------------------------------------------------------------
c
c
c  Define run control values that can change from run to run.
c   (Could be input from a data file at this spot in the
c    code, but why not just explicitly define them here)
c
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

c     ibtime = 1 !to call restart 
c     ibtime = 6967421  !ae187 - Huygens 

c   Restart time indexes
c
c     ietime = 32000000
c
c
c  Total simulation time for this run
c
c      endtime = 3.1536d7* 1000.5               ! 3.1536d7 = 1yr

      endtime = 1800.

c     endtime = 3.1536d7 * 25.0                ! 3.1536d7 = 1yr
c     endtime = 3.1536d7 * 2.0                ! 3.1536d7 = 1yr
c     endtime = 8.64d4 * 180.
c     endtime = 8.64d4 * 50.
c
c
c  Define timestep size [s].
c
c      dtime = 1.d0
c      dtime = 10.      ! 10 seconds
c      dtime = 8.64d4   ! one day

c      dtime = 8.64d4 * 6.6  ! keep largest tholin from falling 
                             ! more than 2 km in a timestep
c      dtime = 6.048d5  ! one week
c      dtime = 3.024d5  ! half week
c      dtime = 3.1536d7 ! one year
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
      ext = 'cd'
c      ext = 'dd'
c
c  Define flag to control temperature profile for model initialization
c    .true.  for Huygens HASI data
c    .false. for Lellouch & Hunten engineering model

      do_huygens_Tprofile = .true.

c     resifil = 'carma_res.jf.in'          ! Restart input file
       
      if( do_huygens_Tprofile ) then
        ! HASI data file 0-625 km; z, T, rho values
        trhoifil1 = '../data/HASI-temperature_density_625km.txt' 
      else
      !if( NZ .eq. 50 ) then              ! Temperature/Density input file
       trhoifil1 = '../data/Trho_adjust.mid' !  1 km - 99 km, +2 km 
      endif

c     tempofil = 'Files/pressure'         ! Output for graphing
c     tempofil = 'Files/diffusion'        ! Output for graphing
c     tempofil = 'Files/velocity1'        ! Output for graphing
c     tempofil = 'Files/coagkernel'       ! Output for graphing
c     tempofil = 'Files/N'                ! Output for graphing
      tempofil = 'Files/growth'//ext//'.p'  ! Output for graphing

      fluxifil = '../data/flux100'        ! Aerosol flux input file

      apcifil1  = '../data/ptclu_502.p'   ! Aerosol particles input file
      apcifil2  = '../data/pctop'         ! Aerosol particles input file

c      prtofil = 'Files/carma'//ext//'.p'  ! Output print file

c      if( do_netcdf )then
c       hisofil = 'tcarma_his'//ext//'.cdf' ! Output history file (netcdf)
c      else
c       hisofil = 'tcarma_his'//ext//'.bin' ! Output history file (binary)
c      endif

c      resofil = 'tcarma_res'//ext//'.out'  ! Restart output file

c      stepofil = 'substep'//ext//'.out'   ! Timestepping output file

c      radofil = 'tcarma_rad'//ext//'.out'  ! Radiation submodel print output
c
c
c  Define frequencies of print and history output:
c  use timestep period (nprint, nhist, nrest) when > 0, otherwise
c  use time period (pprint, phist, prest).
c
c
      nprint = -1         ! timestep period between outputs to print file
      nhist  = -1         ! timestep period between outputs to history file
c     nrest  = ietime     ! timestep period between outputs to restart file
      nrest  = -20

      pprint = 8.64d4     ! time period between outputs to print file [s]
c     pprint = YEAR       ! time period between outputs to print file [s]
c      pprint = -1.

c     write(*,*) 'Changing <phist> in {postep} after 1 year!!!'
c     phist  = YEAR       ! time period between outputs to history file [s]
c     phist  = 5.*YEAR   ! time period between outputs to history file [s]
c     phist  = YEAR/12.   ! time period between outputs to history file [s]
c     phist  = YEAR/5.    ! time period between outputs to history file [s]
c     phist  = 8.64d4     ! time period between outputs to history file [s]
c     phist  = 8.64d4/2.  ! time period between outputs to history file [s]
c     phist  = 8.64d4*10.  ! time period between outputs to history file [s]
c     phist  = 3.60d3     ! time period between outputs to history file [s]
c     phist  = 3.60d3*2.  ! time period between outputs to history file [s]

      phist  = 100000
      prest  = 100000.        ! time period between outputs to restart file [s]
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
      do_print = .true.
      do_hist = .false.
      do_rest = .false.
c
c
c  Define flag to control printing of yearly output
c
c      prt_year = .true.
c
c
c  Define flag to control whether radiative transfer is to be computed
c
      do_rad = .false.
c
c
c  Define flag to control whether layers below the model domain are included in
c  the radiative transfer model
c
      do_below = .false.
c
c
c  Define flag to control whether the model is to be run as a parcel simulation
c  (multiple parcels may be simulated by allowing NZ > 1).
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
c  (evaporation and nucleation also).

      do_grow = .false.
c
c
c  Define flag to control whether temperature is changed by latent heating.
c
c     do_thermo = .false.
c
c
c  Define flag to control whether vertical transport occurs.
c
c     do_vtran = .true.
c
c
c  Define flags to control vertical boundary conditions:
c    <itbnd_pc> = I_FIXED_CONC: use specified concentration at the boundary
c               = I_FLUX_SPEC:  use specified flux
c    <ibbnd_pc>: same as <itbnd_pc>, but for bottom boundary;
c    equivalent parameters are defined for gases and potential temperature
c    boundary conditions.
c
      itbnd_pc  = I_FLUX_SPEC
      ibbnd_pc  = I_FLUX_SPEC
      itbnd_gc  = I_FLUX_SPEC
      ibbnd_gc  = I_FLUX_SPEC
      itbnd_ptc = I_FIXED_CONC
      ibbnd_ptc = I_FIXED_CONC

c
c  Define flag for mass conservation testing: True sets <ftop> and <fbot>
c  to zero in <versol}
c
      test_mass_cons = .false.
c
c
c  Define flag to control whether horizontal transport occurs in the
c  east-west or north-south directions (always .false. for a parcel simulation).
c
      do_ew = .false.
      do_ns = .false.
c
c
c  Define flag to control which horizontal advection algorithm is
c  used:
c    <ihoradv> = I_PPM: Use piecewise polynomial method
c    <ihoradv> = I_GALERKIN: Use Galerkin method with Chapeau functions
c
      ihoradv = I_PPM
c
c
c  Define flag to control whether a variable time-step should be used.
c  Also define min and max time-steps, maximum tolerance in particle
c  concentration changes <dpctol>, maximum tolerance in gas concentration 
c  changes <dgstol>, and minimum relative concentration to consider <conmax>
c  (i.e., for bins with concentrations less than <conmax>*max(pc), we don't
c  worry about how large the changes were).
c
      do_varstep = .false.
c
c
c  Set up things that depend on variable timestepping
c
c      maxsubsteps = 1 !1000 !7200
c      minsubsteps = 1 !2
c      conmax = 1.e-1
c      if( do_varstep )then
c        dtmin  = 1.d-2
c        dtmax  = 8.64d4   ! one day
c        dpctol = 2.
c        dgstol = 0.2
c        conmax = 1.e-3
c      else
        do_step = .true.
c      endif
c
c
c  Define flag to control if error trapping for debugging is to be done
c   (May have no effect on some systems.  Mainly useful for sunos.)
c
c      do_error = .false.
c
c
c  Print some model setup parameters to the screen if they differ from a 
c  default model run
c
      if( do_coag .eqv. .false.) write(*,*) 'do_coag = .false. !!!!!'
      if( do_grow .eqv. .false.) write(*,*) 'do_grow = .false. !!!!!'
      if( do_thermo .eqv. .true.) write(*,*) 'do_thermo = .true. !!!!!'
      if( do_vtran .eqv. .false.) write(*,*) 'do_vtran = .false. !!!!!'
      if( do_varstep .eqv. .false.) write(*,*) 'do_varstep = .false. !!'
      if( test_mass_cons .eqv. .true. ) 
     $  write(*,*) 'Testing mass conservation'
c
c
c  End of per run control values definition (usually no changes below here)
c
c---------------------------------------------------------------------------
c
c
c  Open output print file
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
cc    open(unit=LUNOSTEP,file=stepofil,status='unknown')
c
c
c  Open file for radiation submodel print output
c
c      open(unit=LUNORAD,file=radofil,status='unknown')
c
c
c  Open output file for graphing
c
       open(unit=LUNOTEMP,file=tempofil,status='unknown')
c
c
c  Open input temperatures and air density file
c
       open(unit=LUNITRHO1,file=trhoifil1,status='old')
cc     if( do_huygens_Tprofile ) 
cc   $  open(unit=LUNITRHO2,file=trhoifil2,status='old')
c
c
c  Open input flux file
c
       open(unit=LUNIFLX,file=fluxifil,status='old')
c
c
c  Open input aerosol particle file
c
       open(unit=LUNIAER1,file=apcifil1,status='old')
       open(unit=LUNIAER2,file=apcifil2,status='old')
c
c
c  Announce entry to this routine 
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter init'
c
c
c  Report model name & version tag
c
c      call prtsep
c      write(LUNOPRT,1) PROGNAM, PROGTAG
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
c  Set up execution of write statements for debugging
c
      call testprt
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
c      call initres
       do_netcdf = .true.    !when restarting from carma_2.0 files
       ct(1,1) = 0.986
c      do_thermo = .true.    !when restarting with T perturbation
c      do k=1,NZ
c        read (LUNITAEM,12) tinit(k), rjunk
c        tinit(k) = t3(k)
c
c        if( zl3(k) .lt. 50.d5 ) then        !for 10x tholin case
c          do j=1,NBIN
c            pc3(k,j,1) = 10. * pc3(k,j,1)
c          enddo
c        endif
c
c      enddo
       call restart
      endif

      do igas=1,NGAS
       do k=1,NZ
        gcb(k,igas) = gc3(k,igas)
       enddo
      enddo
      tbegin = time
      tinc = YEAR         !extrapolate over time tinc
      trecov = 2.*YEAR    !time to next extrapolation
      tlex = 3000.*YEAR    !time of last extrapolation  

c
c  Set up some parameters for parcel test
c
      if(do_parcel) then
        simtitle = 'PARCEL TEST'
        ibtime = 0
        ietime = 8.64d3 * 200.
        itime = ibtime
        time = 0.
        endtime = YEAR
        dtime = 10.
        maxsubsteps = 1
        nprint = 43200
        pprint = -1
        nhist = 43200
        phist = -1
        do_vtran = .false.
      
        write(*,*) '{INIT}: Parameters changed for parcel test.'

      endif
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
c  is .true., <NX> is at least 5 if <do_ew> is .true., and <NY>
c  is at least 5 if <do_ns> is .true.
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
     $    'Cannot do north-south transport with NY < 5'
        call endcarma
      endif
c
c
c  Report some initialization values
c
c      write(LUNOPRT,5)
c
c      write(LUNOPRT,2) 'ibtime', ibtime
c      write(LUNOPRT,2) 'ietime', ietime
c      write(LUNOPRT,3) 'endtime', endtime
c
c      write(LUNOPRT,2) 'NX', NX
c      write(LUNOPRT,2) 'NY', NY
c      write(LUNOPRT,2) 'NZ', NZ
c
c      write(LUNOPRT,3) 'time', time
c      write(LUNOPRT,2) 'itime', itime
c      write(LUNOPRT,3) 'dtime', dtime
c
c      write(LUNOPRT,6) 'do_error', do_error
c      write(LUNOPRT,6) 'do_netcdf', do_netcdf
c      write(LUNOPRT,6) 'do_parcel', do_parcel
c      write(LUNOPRT,6) 'do_coag', do_coag
c      write(LUNOPRT,6) 'do_grow', do_grow
c      write(LUNOPRT,6) 'do_thermo', do_thermo
c      write(*,6) 'do_thermo', do_thermo
c      write(LUNOPRT,6) 'do_vtran', do_vtran
c      write(LUNOPRT,6) 'do_ew', do_ew
c      write(LUNOPRT,6) 'do_ns', do_ns
c
c      write(LUNOPRT,6) 'do_rad', do_rad
c      if( do_rad )then
c        if( nrad .gt. 0 )then
c          write(LUNOPRT,2) 'nrad', nrad
c        else
c          write(LUNOPRT,3) 'prad', prad
c        endif
c      endif
c
c      write(LUNOPRT,6) 'do_print', do_print
c      if( do_print )then
c        if( nprint .gt. 0 )then
c          write(LUNOPRT,2) 'nprint', nprint
c        else
c          write(LUNOPRT,3) 'pprint', pprint
c        endif
c      endif
c
c      write(LUNOPRT,6) 'do_hist', do_hist
c      if( do_hist )then
c        if( nhist .gt. 0 )then
c          write(LUNOPRT,2) 'nhist', nhist
c        else
c          write(LUNOPRT,3) 'phist', phist
c        endif
c      endif
c
c      write(LUNOPRT,6) 'do_rest', do_rest
c      if( do_rest) then
c        if( nrest .gt. 0 )then
c          write(LUNOPRT,2) 'nrest', nrest
c        else
c          write(LUNOPRT,3) 'prest', prest
c        endif
c      endif
c
c      call dblank(simtitle, ns)
c      write(LUNOPRT,4) 'simtitle', simtitle(1:ns)
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
c      write(LUNOPRT,9)
c      call prtsep
c
c
c  If parcel test, generate some output
c
c     if(do_parcel) 
c    $  call conserv
c
c
c  Return to caller with model initializations complete
c
      return
      end
