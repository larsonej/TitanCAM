      subroutine endcarma
c
c
c  @(#) endrun.f  Bardeen  Aug-2006
c  This routine prepares for abnormal termination of a model run. It relies on the
c  parent model to actual terminate the run, but cleans up CARMA and signals that
c  the run should be terminated.
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
      use abortutils
      include 'globaer.h'
c
c
c  Define formats
c
    1 format(a,' terminating abnormally')
    2 format('Last timestep index: ',i6)
    3 format('Last time: ',f12.2)
    4 format(a,' output was to file ',a)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter endcarma'
c
c
c  Report last timestep & time
c
c  NOTE: Even with the flush, these messages didn't make it out to carma.p,
c  so I have change this to go to the CAM log (stdout).
      call prtsep
      write(*,1)
      write(*,2) itime
      write(*,3) time
c
c
c  Close files for time step diagnostic output
c
c     call flush(LUNOSTEP)
c     close(unit=LUNOSTEP)
c
c
c  Close output print file
c
c  NOTE: The flush and close didn't work with MPI, so now LUNOPRT goes to
c  stdout.
c
      call prtsep
c      call flush(LUNOPRT)
c      close(unit=LUNOPRT)
c
c
c  Terminate the program's execution.
c  
c  NOTE: This is an f90 routine.
c
c
c      stop 1
      call endrun("CARMA aborting the run ...")
      
      end
