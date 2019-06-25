      subroutine step
c
c
c  @(#) step.f  McKie  Oct-1995
c  This routine performs all calculations necessary to
c  take one timestep.
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
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter step'
c
c
c  Do pre-timestep processing
c
      call prestep
c
c
c  Update model state at new time
c
      call newstate
c
c
c  Modify time-step if necessary
c
      if( do_varstep ) call varstep
c
c
c  Do post-timestep processing
c
      if( do_step ) call postep 
c
c
c  Return to caller with one timestep taken
c
      return
      end
