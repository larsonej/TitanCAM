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
c  Define formats
c
    1 format('Cloud: ',1pe16.9,'  Ethane: ',1pe15.8,
     $       '  Methane: ',1pe15.8,'  Core: ',1pe15.8,i4) 
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
c
c  Determine particle charging and set up new coagulation kernels
c
       call newckern
c
c
c  Calculate (implicit) particle loss rates for coagulation.
c
        call coagl
c
c
c  Calculate particle production terms and solve for particle 
c  concentrations at end of time step.
c
        do ielem = 1,NELEM
          do ibin = 1,NBIN
            call coagp(ibin,ielem)
            call csolve(ibin,ielem)
          enddo
        enddo
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