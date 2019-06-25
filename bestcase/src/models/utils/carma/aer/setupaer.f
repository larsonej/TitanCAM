       subroutine setupaer
c
c
c  @(#) setupaer.f  Ackerman  Jan-1996
c  This routine calls all the other setup routines to calculate other
c  time-independent parameters for aerosol and cloud microphysics. The
c  defintion of the arrays that was formerly done in this routine is now
c  done in defineaer.f.
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
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer'
!      write(*,*) 'enter setupaer'
c
c
c  Evaluate fall velocities.
c
      call setupvf

c
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for condensational growth and evaporation.
c
      if ( do_grow ) then
        call setupgrow
        call setupgkern

c
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for nucleation.
c
        call setupnuc
      endif
      
c
c-------------------------------------------------------------------------------
c
c
c==Set options for coagulation kernel:
c
      if( do_coag )then

c
c
c  Evaluate derived coagulation mapping arrays and kernels.
c
        call setupckern
        call setupcoag

      endif
c
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined.
c
      return
      end