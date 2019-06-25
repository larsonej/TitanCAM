      subroutine microslow
c
c
c  @(#) microslow.f  McKie  Sep-1997
c  This routine drives the potentially slower microphysics calculations.
c
c  Originally part of microphy.  Now in this separate routine to allow
c  time splitting of coagulation at a different timestep size from
c  other microphysical calcs.
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
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter microslow'
c
c
c  Do coagulation if requested
c
      if( do_coag )then
c
c  ! from E.Barth
c  Determine particle charging and set up new coagulation kernels
c
!     open(unit=99,file='ckern1.txt',status='unknown')
!	  write(99,*) ckernel(24,:,:,1,1)
!	  close(unit=99)
  
!       call newckern
   
!        open(unit=99,file='ckern2.txt',status='unknown')
!	  write(99,*) ckernel(24,:,:,1,1)
!	  close(unit=99)
c
c  Calculate (implicit) particle loss rates for coagulation.
c
        call coagl
c
c
c  Calculate particle production terms and solve for particle 
c  concentrations at end of time step.
c
!        do ielem = 1,NELEM
!          do ibin = 1,NBIN
!            call coagp(ibin,ielem) !EJL 6-8-10
!            call csolve(ibin,ielem)
!          enddo
!        enddo
c
c
c  End of conditional coagulation calcs
c
      endif
c
c
c  Return to caller with new particle concentrations.
c
      return
      end
