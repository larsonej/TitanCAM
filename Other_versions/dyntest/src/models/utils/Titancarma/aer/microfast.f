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
c  Coefficients for 3rd order polynomial fit to freezing temperature
c  data for methane+nitrogen system, from Omar et al. 1962, Table I
      dimension tcoeff_ch4(4), xcoeff_n2(6)

      data tcoeff_ch4/ 52.4593, 54.4942, -39.0671, 22.9201/

      data xcoeff_n2/ -0.054949131, 0.0020078157, -4.1577243e-06,
     $                3.9162634e-09, -1.7868039e-12, 3.1424387e-16/
c
c  Define formats
c
    1 format(a,i6,' cloud: ',1pe13.6,'  vapor: ',1pe13.6)
    2 format(4(i3,1pe13.6,a3,'|'))
    3 format(i3,1pe13.6)
    5 format(a,i4,2x,'gc/vol:',f15.10,2x,'gc/cld:',f15.10,2x,
     $       'cmf:',f15.10,2x,'cloud(pc):',1pe13.6,i4)
    6 format(a,i4,2x,'primary: ',1pe15.8,'growcore: ',1pe15.8)
    7 format(i4,2x,'Primary: ',1pe11.4,'  Growcore: ',1pe11.4,
     $       2(2x,1pe11.4),'  add mass frac: ',1pe11.4)
    8 format(2(i4),'  Nucl: ',1pe13.6,'  Grow: ',1pe13.6,
     $       '  Evap: ',1pe13.6,'  ss: ',0p,f10.5,i10)
    9 format(a,'  C2H6(loss/prod): ',1pe13.6,1pe13.6,
     $       '  Core(loss/prod): ',1pe13.6,1pe13.6,i10)
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter microfast'
c
c
c  Do growth microphysical processes if requested
c
      if( do_grow )then

c
c  Update maximum particle concentration for each spatial grid box
c  (in units of cm^-3)
c
       do igroup = 1,NGROUP 
        iepart = ienconc(igroup)
        pconmax(ixyz,igroup) = SMALL_PC
        do ibin = 1,NBIN
          pconmax(ixyz,igroup) = max( pconmax(ixyz,igroup),
     $                                pc3(ixyz,ibin,iepart) )
        enddo
        pconmax(ixyz,igroup) = pconmax(ixyz,igroup) /
     $        ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )
       enddo
c
c
c  Calculate (implicit) particle loss rates for nucleation, growth,
c  evaporation, melting, etc.
c

c  Calculate new freezing point temperatures based on methane content of droplet
c
        Tfreez(1) = 80.6
c       igas = 1
c       yP = p3(ixyz)/1000. !pressure in mbar
c       xfrac = xcoeff_n2(1) + xcoeff_n2(2)*yP +
c    $           xcoeff_n2(3)*yP**2 + xcoeff_n2(4)*yP**3 +
c    $           xcoeff_n2(5)*yP**4 + xcoeff_n2(6)*yP**5
c       xfrac = ONE - xfrac
c       Tfreez(igas) = tcoeff_ch4(1) + tcoeff_ch4(2)*xfrac + 
c    $                     tcoeff_ch4(3)*xfrac**2 + 
c    $                        tcoeff_ch4(4)*xfrac**3

cc      Tfreez(igas) = tcoeff_ch4(1) + tcoeff_ch4(2)*xfrac(igas,ig) + 
cc   $                     tcoeff_ch4(3)*xfrac(igas,ig)**2 + 
cc   $                        tcoeff_ch4(4)*xfrac(igas,ig)**3

        call zeromicro
        call supersat
        call growevapl
        call vapordepl  ! Nucleation from Pruppacher & Klett vapor deposition eqns
        call actdropl   ! Droplet activation (constant rate)
        call freezaerl
        call freezdropl ! Droplet freezing (constant rate)
        call melticel   ! Melt ice crystal (constant rate)
        call coremelt   ! Adjust ice core part of mixed phase particle 

        call freezmixedl ! Freeze droplet part of mixed phase particle (const rate)
        call meltmixedl  ! Melt ice core part of mixed phase particle (const rate)
c
c
c  Calculate particle production terms and solve for particle 
c  concentrations at end of time step.
c
        do ielem = 1,NELEM
          do ibin = 1,NBIN
            call growp(ibin,ielem)
            call upgxfer(ibin,ielem)
            call psolve(ibin,ielem)
          enddo
        enddo
 
c
c  Calculate particle production terms for evaporation.
c
        call evapp
        call downgxfer
c
c
c  Calculate gas loss rates and production terms due to particle nucleation,
c  growth, and evaporation.
c
        call gasexchange
c
c
c  Apply evaporation production terms to particle concentrations.
c
        call downgevapply
c
c
c  Solve for gas concentrations at end of time step.
c
        call gsolve
c
c
c  End of conditional microphysical growth processing
c
      else
c  Add tholin production term here if not calling {psolve}
cc     write(*,*) itime,' Haze Production:',ixyz,pc3(ixyz,1,1),
cc   $             pc3(ixyz,1,1) + dtime*rprod(ixyz)
       pc3(ixyz,1,1) = pc3(ixyz,1,1) + dtime*rprod(ixyz)

      endif !do_grow
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
      call supersat
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
