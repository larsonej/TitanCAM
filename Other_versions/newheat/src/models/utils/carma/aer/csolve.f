       subroutine csolve(ibin,ielem)
c
c
c  @(#) csolve.f  McKie  Sep-1997
c  This routine calculates new particle concentrations from coagulation
c  microphysical processes.
c
c   The basic form from which the solution is derived is
c   ( new_value - old_value ) / dtime = source_term - loss_rate*new_value
c
c  This routine derived from psolve.f code, in which particle concentrations
c  due to coagulation were formerly included, before the relatively slow
c  coagulation calcs were separated from the other microphysical processes
c  so that time splitting could be applied to these fast & slow calcs.
c
c  Argument list input:
c    ielem, ibin
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter csolve'
c
c
c  Define current group & particle number concentration element indices
c
      igroup = igelem(ielem)         ! particle group
      iepart = ienconc(igroup)       ! particle number concentration element
c
c
c  Visit each spatial point
c
      do ixyz = 1,NXYZ
c
c
c  Metric scaling factor
c
        xyzmet = xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)
c
c
c  Compute total production rate due to coagulation
c
        ppd = coagpe(ixyz,ibin,ielem) / xyzmet
c EJL 6-9-10 comment out ppd, pls, pc3
c
c  Compute total loss rate due to coagulation
c
        pls = coaglg(ixyz,ibin,igroup) / xyzmet
c
c
c  Update net particle number concentration during current timestep
c  due to production and loss rates for coagulation
c
        pc3(ixyz,ibin,ielem) = ( pc3(ixyz,ibin,ielem) + dtime*ppd ) /
     $                         ( ONE + pls*dtime )

      enddo
c
c
c  Return to caller with new particle number concentrations.
c
      return
      end
