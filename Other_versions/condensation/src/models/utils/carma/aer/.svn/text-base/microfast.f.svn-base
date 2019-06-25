      subroutine microfast
c
c
c  @(#) microfast.f  Jensen  Oct-1995
c  This routine drives the microphysics calculations.
c
c  Modified Sep-1997  (McKie)
c    Name changed from microphy to microfast.
c    Remove coag{l,p} calls moved to separate routine microslow.
c    Both micro{fast,slow} called from newstate.
c    Allows time-splitting slower coagulation from other microphysics.
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter microfast'
c
c
c  Set production and loss rates to zero.
c
      call zeromicro
c
c
c  Calculate (implicit) particle loss rates for nucleation, growth,
c  evaporation, melting, etc.
c
      if( do_grow )then
        call supersat
        call growevapl
        call actdropl
        call freezaerl
        call freezdropl
        call melticel
        call coremelt
        call freezmixedl
        call meltmixedl
      endif
c
c
c  Calculate particle production terms and solve for particle 
c  concentrations at end of time step.
c
      do ielem = 1,NELEM
        do ibin = 1,NBIN

          if( do_grow )then
            call growp(ibin,ielem)
            call upgxfer(ibin,ielem)
          endif

          call psolve(ibin,ielem)

        enddo
      enddo
c
c
c  Calculate particle production terms for evaporation;
c  gas loss rates and production terms due to particle nucleation;
c  growth, and evaporation;
c  apply evaporation production terms to particle concentrations;
c  and solve for gas concentrations at end of time step.
c
      if( do_grow )then
        call evapp
        call downgxfer
        call gasexchange
        call downgevapply
        call gsolve
      endif
c
c
c  Update temperature if thermal processes requested
c
      if( do_thermo )then
        call tsolve
      endif
c
c
c  Update saturation ratios (note: the simulated microphysical
c  processes do not depend on these updated values)
c
      if( do_grow .or. do_thermo )then
        call supersat
      endif
c
c
c  Update particle densities
c
      call rhopart
c
c
c  Return to caller with new particle and gas concentrations.
c
      return
      end
