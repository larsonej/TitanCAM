      subroutine zeromicro
c
c
c  @(#) zeromicro.f  Ackerman  Oct-1997
c  This routine zeroes the fast microphysics sinks and sources, 
c  at one spatial point per call.
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
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter zeromicro'
c
c
c  Set production terms and loss rates due to nucleation, growth,
c  and evaporation to zero.  Also set index of smallest bin nuceleated
c  during time step equal to <NBIN> first time through spatial loop.
c
      rlheat = 0.

      do igas = 1,NGAS
        gasprod(igas) = 0.
        do igrp = 1,NGROUP
          gprod_grow(igrp,igas) = ZERO
          gprod_drop(igrp,igas) = ZERO
          gprod_mono(igrp,igas) = ZERO
          gprod_poly(igrp,igas) = ZERO
        enddo
      enddo

      do ielem = 1,NELEM
        do i = 1,NBIN
          rnucpe(i,ielem) = 0.
          growpe(i,ielem) = 0.
          evappe(i,ielem) = 0.
        enddo
      enddo

      do igroup = 1,NGROUP

        if( ixyz .eq. 1 )then
          inucstep(igroup) = NBIN
        endif

        do i = 1,NBIN
          growlg(i,igroup) = 0.
          evaplg(i,igroup) = 0.
          do igto = 1,NGROUP
            rnuclg(i,igroup,igto) = 0.
          enddo
        enddo
      enddo
c
c
c  Return to caller with fast microphysics sinks and sources zeroed.
c
      return
      end
