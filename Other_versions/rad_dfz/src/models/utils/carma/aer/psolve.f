       subroutine psolve(ibin,ielem)
c
c
c  @(#) psolve.f  Jensen  Oct-1995
c  This routine calculates new particle concentrations.
c
c   The basic form from which the solution is derived is
c   ( new_value - old_value ) / dtime = source_term - loss_rate*new_value
c
c  Modified  Sep-1997  (McKie)
c    New particle concentrations due to coagulation processes
c    were moved to the csolve routine.  Csolve is called to
c    update particle concentrations due to coagulation.
c    This new psolve now updates particle concentrations due
c    to the faster calcs of the non-coag microphysical processes.
c
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter psolve'
c
c
c  Define current group & particle number concentration element indices
c
      igroup = igelem(ielem)
      iepart = ienconc(igroup)
c
c
c  Compute total production rate
c
        ppd = rnucpe(ibin,ielem) 
     $      + growpe(ibin,ielem) 
     $      + evappe(ibin,ielem) 
c
c
c  Sum up nucleation loss rates
c
        rnuclgtot = 0.
        do igto = 1,NGROUP
          rnuclgtot = rnuclgtot + rnuclg(ibin,igroup,igto)
        enddo
c
c
c  Compute total loss rate
c
        pls = rnuclgtot
     $      + growlg(ibin,igroup) 
     $      + evaplg(ibin,igroup) 

!EJL - CARMA rainout from E.Barth's titan code 9-28-10
!calculated as 1/time   .17 = 30km
        if(zl3(ixyz) .gt. 0.17) pls = pls + 6.337618d-10    !6.338e-10 is 50 year rainout time
c
c
c  Update net particle number concentration during current timestep
c  due to production and loss rates
c
        pc3(ixyz,ibin,ielem) = ( pc3(ixyz,ibin,ielem) + dtime*ppd ) /
     $                         ( ONE + pls*dtime )
c
c
c  Prevent particle concentrations from dropping below SMALL_PC
c
         call smallconc(ibin,ielem)
c
c
c  Return to caller with new particle number concentrations.
c
      return
      end
