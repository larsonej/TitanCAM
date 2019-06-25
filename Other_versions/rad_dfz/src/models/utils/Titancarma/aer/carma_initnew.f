      subroutine carma_initnew
c
c
c  @(#) initnew.f  McKie  Oct-1995
c  This routine performs a cold start initialization for a new
c  model simulation.
c
c  Note:  Array dimensions such as NX, NY, NZ, ... are
c  defined in the globaer.h file, and also serve as the
c  upper limits to index values for looping through arrays.
c
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables.
c
      include 'globaer.h'
c
c
c  Define formats.
c
    1 format(/,'Doing a cold start initialization for a new simulation')
c
c---------------------------------------------------------------------------
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter initnew'
c
c
c  Open output print file
c
      ext = 'xx'
      prtofil = 'Files/carma'//ext//'.p'  ! Output print file
      open(unit=LUNOPRT,file=prtofil,status='unknown')
c
c
c  Report cold start initialization being done.
c
      write(LUNOPRT,1)
c
c
c  Define simulation title.
c
      simtitle = 'ETHANE / METHANE clouds
     $            GCM bin model Test'
c
c
      do_vtran = .true.
      do_coag = .true.
      do_grow = .true.
c
c  Initialize timestep index.
c
      itime = 0 
c
c
c  Initialize simulation time [s].
c
      time = 0.
c
c
c  Initialize atmospheric structure.
c
cc    call carma_initatm
      call carma_initmp
c
c
c  Define mapping arrays and time-independent parameters for
c  aerosol and cloud microphysics.
c
      call setupaer
c
c
c  Initialize particle concentrations.
c
      call initaer
c
c
c  Initialize gas concentrations.
c
      call initgas
c
c
c  Return to caller with cold start initialization complete.
c
      return
      end
