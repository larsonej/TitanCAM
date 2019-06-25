      subroutine vapordepl
c
c
c  @(#) vapordepl.f  Barth  
c  This routine evaluates particle loss rates due to nucleation <rnuclg>:
c  vapor deposition only.
c  
c  The loss rates for all particle elements in a particle group are equal.
c
c  To avoid nucleation into an evaporating bin, this subroutine must
c  be called after growevapl, which evaluates evaporation loss rates <evaplg>.
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c  General structure taken from {actdropl} but modified to include nucleation
c  equations from Pruppacher and Klett
c
c  Argument list input:
c
c  Argument list output:
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c
      logical evapfrom_nucto, print_message

      dimension fm(4)
c
c
c  Define formats
c
    1  format(f4.1,i3,1p3e12.4)
    2  format('nucleation', 3x,i3,'km',3x,f8.4,i11,f12.2,3x,1p4e12.4)
    3  format(f4.1,3x,3(i3,3x),2(1p2e12.4))
    4  format(f12.2,3x,f4.1,3x,2(1pe12.4))
    5  format(6(1pe12.4))
    6  format(a,2x,3(i5),2x,'rnuclg: ',1pe13.6,' tholin: ',1pe13.6)
    7  format(a,f12.4,2x,i9,2x,i5,2x,'rnuclg: ',2(1pe11.4,3x))
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter vapordepl'
cc    return  !!! Test - turn off all nucleation
c     if( itime.gt.10000 ) return
c
c
      print_message = .false.
c
c
c  Loop over particle groups.
c
      do igroup = 1,NGROUP
       do igs = 1,NGAS

        igas = inucgas(igs,igroup)                ! condensing gas
        iepart = ienconc( igroup )            ! particle number density element
 
        if( igas .ne. 0 )then
c       if( mod(itime,2) ) then      !Test - only one gas nucl in a timestep
c        igasturn = 1
c       else
c        igasturn = 2
c       endif
c       if( igas .eq. igasturn ) then
c       if( igas.eq.1) then
c
c
c  Calculate nucleation loss rates.  Do not allow nucleation into
c  an evaporating bin.
c
c  <ienucto> is index of target nucleation element;
c  <ignucto> is index of target nucleation group.
c
        do inuc = 1,nnuc2elem(iepart)

         ienucto = inuc2elem(inuc,iepart)
         if( ienucto .ne. 0 )then
           ignucto = igelem( ienucto )
         else
           ignucto = 0
         endif
c
c
c  Only compute nucleation rate for vapor deposition
c
         if( inucproc(iepart,ienucto) .eq. I_VAPORDEP ) then
c
c
c  Loop over particle bins.  Loop from largest to smallest for 
c  evaluation of index of smallest bin nucleated during time step <inucstep>.
c
c    Since tholin only cover 35 bins, don't want to nucleate from larger
c    bins; but want to consider all bins for clouds
          if( iepart.eq.1 ) then
           lbin = NCCNBIN 
          else 
           lbin = NBIN
          endif

          do ibin =lbin,1,-1
c
c
c  <inucto> is index of target nucleation bin.
c
           if( ignucto .ne. 0 )then
             inucto = inuc2bin(ibin,igroup,ignucto)
           else
             inucto = 0
           endif
c
c
c  Bypass calculation if few particles are present 
c
cc         if( pconmax(ixyz,igroup) .gt. FEW_PC )then
           if( pc3(ixyz,ibin,iepart) .gt. FEW_PC )then
c
c
c  Set <evapfrom_nucto> to .true. when target droplets are evaporating
c  
            if( inucto .ne. 0 )then
              evapfrom_nucto = evaplg(inucto,ignucto) .gt. 0.
            else
              evapfrom_nucto = .false.
            endif
c
c  Set ice or liquid versions of <supsat>, <pvap>, and <surfct> to be used
c
            if( t3(ixyz) .le. Tfreez(igas) ) then
cc           return !test - no ice cloud nucleation

             ! No nucleation if target group is single-phase droplet
             ! (too cold)
             if((.not. is_grp_mixed_phase(ignucto)) .and.
     $            (.not. is_grp_ice(ignucto)))  return  
             ! Otherwise, set parameters for ice clouds
              ss = supsati3(ixyz,igas)
              pvap = pvapi3(ixyz,igas)
              surfct = surfctia(ixyz,igas)
            else
             ! No nucleation if target group is single-phase ice cloud
             ! (too warm)
             if( (.not. is_grp_mixed_phase(ignucto)) .and. 
     $               is_grp_ice(ignucto) )  return  
             ! Otherwise, set parameters for droplet clouds
              ss = supsatl3(ixyz,igas)
              pvap = pvapl3(ixyz,igas)
              surfct = surfctwa(ixyz,igas)
            endif
            
            if( ( ss .gt. scrit(iz,ibin,igroup) ).and.
     $           .not. evapfrom_nucto )then

             ielem = ienconc(igroup)
             ag=2.*gwtmol(igas)*surfct/
     $             (RGAS*t3(ixyz)*rhoelem(ielem)*log(ss+1.))       !RHO_I


c
c  Nucleation for germ growth by direct vapor deposition (PK 9-8)
c
c
c   ..for spherical surface:
c

      fx=r(ibin,igroup)/ag

        IF(fx.gt.3000.) fx = 3000.
        Phi=(1.-2.*ct(igas,igroup)*fx+fx*fx)**(1./2.)
        fm(1)=1.+((1.-ct(igas,igroup)*fx)/Phi)**3
        fma=fx*fx*fx*2.
        fmb=fx*fx*fx*(-3.*(fx-ct(igas,igroup)) /
     $              Phi+((fx-ct(igas,igroup))/Phi)**3)
        fm(2)=fma*1.d0+fmb*1.d0
        fm(3)=3.*ct(igas,igroup)*fx*fx*((fx-ct(igas,igroup))/Phi-1.)
        fm(4)=fm(1)+fm(2)+fm(3)
        fm(4)=(fm(4)/2.)

           if(itime.eq.14695 .and. ixyz.eq.3 .and.
     $         inuc2bin(ibin,igroup,ignucto).gt.34 )
     $         write(*,*) itime,' vapordepl:',ag,ss,t3(ixyz),fm(4),
     $             dFg,gk**2,dFg/(3*PI*BK*t3(ixyz)*gk**2)
c
c   ..for flat surface:
c
c       fm(4)=(2.+ct(igas,igroup))*(1.-ct(igas,igroup))**2/4.
c
        dFg=4.*PI*ag**2*surfct/3.

        rMw = gwtmol(igas)/AVG

          C1 = 1.d15
          gk=4.*PI*ag**3/(3.*vol(ibin,igroup))
          con1=pvap*(ss + 1.)
     $          / sqrt(2.*PI *rMw *BK *t3(ixyz))
          Zs=sqrt(dFg/(3*PI*BK*t3(ixyz)*gk**2))

        dFg = dFg * fm(4)

c  ..And here is the prefactor <Pfact>:
c
        Pfact = 4.*(PI**2)*(r(ibin,igroup)**2)*(ag**2)*con1*Zs*C1
c
c  ..And finally here is the Nucleation rate <rnuclg>:
c
      rnuclg(ibin,igroup,ignucto)=Pfact*EXP(-dFg/(BK*t3(ixyz)))

c  Test - Turn off methane-tholin nucl
cc    if(ignucto.eq.4) rnuclg(ibin,igroup,ignucto)=ZERO
cc    if(ignucto.eq.2) rnuclg(ibin,igroup,ignucto)=ZERO
cc    if(ignucto.eq.3) rnuclg(ibin,igroup,ignucto)=ZERO

c     if(igroup.eq.1)
c    $  write(*,*) ibin,rnuclg(ibin,igroup,ignucto)
c     if(ixyz.eq.iwa)
c    $ write(*,*) 'Nucl:',igas,rnuclg(ibin,igroup,ignucto),itime


cc      if(rnuclg(ibin,igroup,ignucto).gt.SMALL_PC .and. do_vtran) then
cc        dtmax = 3.6d3 
cc        if( pc3(ixyz,35,2) .gt. 1.d-25 ) dtmax = 1.2d3
cc      endif
        
c       if(print_message) then
c       if(rnuclg(ibin,igroup,ignucto) .ge. SMALL_PC) then
c               write(*,2) ixyz*2-2,r(ibin,1)*1.d4,itime,
c    $                 time/60.**2/24.,
c    $                 dtime,pc3(ixyz,ibin,1),
c    $                 rnuclg(ibin,igroup,ignucto),
c    $                 ss
c       endif
c       endif
c        ..........
c        write(LUNOTEMP,3) zl3(ixyz)/1.d5,ibin,igroup,itime,
c    $            rnuclg(ibin,igroup,ignucto),ss,
c    $            pc3(ixyz,ibin,1)
c
c       enddo     !supsati
c        ..........
c
c            if(rnuclg(ibin,igroup,ignucto) .ne. 0. .and.
             rnp = rnuclg(ibin,igroup,ignucto)
     $                  *pc3(ixyz,ibin,iepart)*dtime 

c     if( rnp .gt. SMALL_PC )
c    $write(*,*) itime,' Nucl:',ixyz,ibin,ss,rnuclg(ibin,igroup,ignucto)

             if(rnp .gt. SMALL_PC .and.
     $          ibin .lt. inucstep(igroup) ) then
              inucstep(igroup) = ibin
             endif

            endif  ! supsat

           endif   ! pconmax(ixyz,igroup) .gt. FEW_PC
          enddo      ! ibin = 1,NBIN
         endif       ! inucproc(iepart,ienucto) .eq. I_VAPORDEP
        enddo       ! inuc = 1,nnuc2elem(iepart)
       endif        ! (igas = inucgas(igroup)) .ne. 0 
       enddo        ! igs = 1,NGAS
      enddo         ! igroup = 1,NGROUP
c     if(ixyz.eq.iwa)
c    $ write(*,*) '-------------------------------------------------'
c
c
c  Return to caller with particle loss rates due to nucleation evaluated.
c
      return
      end
