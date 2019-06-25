      subroutine growevapl
c
c
c  @(#) growevapl_multi.f  Barth  Jan-2003
c  from growevapl.f  Ackerman  Dec-1995
c  This routine evaluate particle loss rates due to condensational
c  growth and evaporation for all condensing gases.  Handles growth
c  of multiple gases onto a particle
c
c  The loss rates for each group are <growlg> and <evaplg>.
c
c  Units are [s^-1].
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
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
      dimension otherm(NELEM), ieother(NELEM),
     $  dela(NBIN), delma(NBIN), aju(NBIN),
     $  ar(NBIN), al(NBIN), a6(NBIN),dHv(NGAS),
     $  alpha(NBIN,NGAS),beta(NBIN,NGAS),gamma(NBIN,NGAS),
     $  dmdt(NBIN),gromat(NGAS,NGAS+1),xx(NGAS+1),
     $  fracmol(NELEM),elemmol(NELEM),
     $  dmdt_grow(NBIN),dmdt_evap(NBIN),rnormfrac(NBIN,NELEM)
c
c
c  Define formats
c
    1 format(2i5,3x,f7.3,3x,f6.4,3(3x,1pe13.6))
    2 format(a,f5.2,1pe13.6,i6,i3,1pe13.6)
    3 format(a,1pe11.4,' dmdt: ',1pe11.4,' pc(i+1): ',1pe11.4,
     $       ' al: ',1pe11.4,' dela: ',1pe11.4,' a6: ',1pe11.4,
     $       ' x: ',1pe11.4,' bin ',i3)

    4 format(a,1pe11.4,' dmdt: ',1pe11.4,' com2: ',1pe11.4,
     $       ' bin ',i3)
    5 format(a,1pe11.4,' dmdt_c2h6: ',1pe11.4,' dmdt_ch4: ',
     $       1pe11.4,' bin',i3,' dmdt_gro: ',1pe11.4,
     $       '  dmdte_gro: ',1pe11.4,' frac: ',0p,f10.8)
    6 format(a,i4,2x,'growlg: ',1pe11.4,i4,2x,'evaplg: ',
     $       1pe11.4,2x,'previous: (g) ',1pe11.4,' (e) ',
     $       1pe11.4)
    7 format(a,'  alpha: ',1pe11.4,'  beta: ',1pe11.4,
     $       '  gamma: ',1pe11.4,'  dHv: ',1pe11.4)
    8 format(a,i4,'  dmdt: ',1pe11.4,'  dmdt_pr: ',1pe11.4,
     $       '  dmdt_gc: ',1pe11.4,'  dmdt_gro: ',1pe11.4,
     $       '  pc: ',1pe11.4)
    9 format(a,2(i6),3(3x,1pe13.6))
   10 format(a,2(i6),5(3x,1pe13.6))
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter growevapl'
c
c
c  Calculate particle fluxes due to growth and evaporation for all
c  particle groups.
c
c  Treat condensation of gas <igas> to/from particle group <igroup>.
c
      do igroup = 1,NGROUP

       iepart = ienconc(igroup)     !number concentration element 
       imaingas = igrowgas(iepart)  !primary condensing gas
       
       small = rmass(1,igroup) * SMALL_PC

       if( imaingas .ne. 0 ) then
c
c
c  Bypass calculation if few particles are present 
c
        if( pconmax(ixyz,igroup) .gt. FEW_PC )then
c
c
c        NO SOLUTE STUFF IN HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
         argsol = ZERO 
c
          do ibin = 1,NBIN-1
c
c  <binmass> is mass concentration of cloud particle.
c  <condmass> is mass concentration of primary condensed gas <igas> in particle.
c  <volmass> is mass concentration of all condensed gases <igas> in particle.
c  <suminvol> is total mass concentration of involatile cores in particle group.
c  <sumcores> is total mass concentration of all cores in particle group.
c
           binmass = pc3(ixyz,ibin,iepart) * rmass(ibin,igroup)
           fracmol(iepart) = ONE

           suminvol = ZERO
           sumcores = ZERO

c  First, loop over elements in group to get mass sums
           do ielemg = 2,nelemg(igroup) 
            iesub = iepart - 1 + ielemg  

            if( itype(iesub) .eq. I_COREMASS )then

              suminvol = suminvol + pc3(ixyz,ibin,iesub)
              sumcores = sumcores + pc3(ixyz,ibin,iesub)
              elemmol(iesub) = ZERO

            elseif( itype(iesub) .eq. I_GROWCORE ) then

              isubgas = igrowgas(iesub)
              sumcores = sumcores + pc3(ixyz,ibin,iesub)
              elemmol(iesub) = pc3(ixyz,ibin,iesub) / gwtmol(isubgas)

            elseif( itype(iesub) .eq. I_CORE2MOM ) then
              elemmol(iesub) = ZERO

            endif
           enddo !elements in group

           volmass = binmass - suminvol

c  Next, loop over elements in group to get mass fraction(s) -> m'/m
           do ielemg = 2,nelemg(igroup)
            iesub = iepart - 1 + ielemg

            if( itype(iesub) .eq. I_GROWCORE ) then
             if(volmass .lt. small) then
              rnormfrac(ibin,iesub) = ZERO
             else
              rnormfrac(ibin,iesub) = pc3(ixyz,ibin,iesub)/ volmass
             endif

c            if( do_write(1) ) then
c             if(rnormfrac(ibin,iesub) .gt. ONE) then
c              write(*,2) 'gcmf ',zl3(ixyz)/1.d5,time/24./60.**2,
c    $                  itime,ibin,rnormfrac(ibin,iesub)
c             endif
c            endif

             rnormfrac(ibin,iesub) = 
     $         max( ZERO, min(ONE, rnormfrac(ibin,iesub)) )

            endif
           enddo

           condmass = binmass - sumcores

           elemmol(iepart) = condmass / gwtmol(imaingas)

           totmol = ZERO
           do ielemg = 1,nelemg(igroup)
             iesub = iepart - 1 + ielemg
             totmol = totmol + elemmol(iesub)
           enddo

           do ielemg = 1,nelemg(igroup)
             iesub = iepart - 1 + ielemg
             fracmol(iesub) = elemmol(iesub) / totmol
           enddo
c
c   Now, loop over elements in group to set up dmdt terms (gro kernels)
           do ielemg = 1,nelemg(igroup)
c                     ------------------>need to include numconc elem
c                                        here since we are setting up
c                                        for dmdt calculations
            ielem = iepart - 1 + ielemg
            igas = igrowgas(ielem)

            if(igas .ne. 0) then
c
c  <akas> is combined kelvin (curvature) and solute factors.
c
c  Ignore solute factor for ice particles.
c
             if( t3(ixyz) .le. Tfreez(igas) ) then !is_grp_ice(igroup) ) then
              expon = akelvini(iz,igas) / rup(ibin,igroup)
             else
              expon = akelvin(iz,igas) / rup(ibin,igroup) - argsol 
             endif
             expon = max(-POWMAX, expon)
             akas = exp( expon )
c 
c 
c  Trick for removing haze droplets from droplet bins:
c  allows haze droplets to exist under supersaturated conditions;
c  when below supersaturation, haze droplets will evaporate.
c
             if( (.not. is_grp_ice(igroup)) .and. (akas .lt. 1.) .and.
     $        (supsatl3(ixyz,igas) .lt. ZERO) )  akas = 1.
c
c
c  <dmdt> is growth rate in mass space [g/s].
c
cc Test - set akas to one
cc           akas = ONE
cc End Test
             g0 =  gro(iz,ibin,ielem) 
             g1 = gro1(iz,ibin,ielem) 
             g2 = gro2(iz,ielem)
             g3 = gro3(iz,ielem)
 
            !if( is_grp_ice(igroup) )then
             if( t3(ixyz) .le. Tfreez(igas) )then
              ss = supsati3(ixyz,igas)
              pvap = pvapi3(ixyz,igas)
c             if(itime.gt.17855.and.ixyz.eq.31.and.igroup.eq.2)
c    $         write(*,9) 'growevapl:',itime,ibin,ss,akas,ss-akas
             else
              ss = supsatl3(ixyz,igas)
              pvap = pvapl3(ixyz,igas)
             endif                         !!!! Take <fracmol> out when not 
                                           !!!! mixing gases
              alpha(ibin,igas) = g0 * pvap / (ONE + g0*g1*pvap)
              beta(ibin,igas) = ss + ONE - akas !*fracmol(ielem)
              gamma(ibin,igas) = g1 * g2 * akas !*fracmol(ielem)
            
c            if( do_write(2) ) then
c            if(ixyz.eq.iwa) then
c             write(*,1) ixyz,ibin,time/60.**2/24.,
c    $        akas,alpha(ibin,igas),beta(ibin,igas),gamma(ibin,igas)
c             write(*,1) ixyz,ibin,time/60.**2/365.25,
c    $         akas,g0,g1,g3
c            endif
c            endif

             dHv(igas) = ONE/g2   !Latent heat

            endif !igas .ne. 0
           enddo !ielemg=1,nelemg(igroup)

c  dmdt's for individutal gases (ignoring dependence on each other)
           do i = 1,NGAS
             xx(i) = alpha(ibin,i) * beta(ibin,i)
           enddo !igas
c
c
c  Growth rate for particle in <ibin> is sum of growth rates for each gas
c
c  Calculate dm/dt <dmdt_gro> and dm'/dt <dmdte_gro> for use later to get 
c  mass fraction rate of change.  Also use dm/dt <dmdt> for the cloud  particle 
c  for calculation of growth/evap loss terms
c
     !!!!!! Test !!!!!!!
     !!!    xx(1) = -2.d-13
     !!!    xx(2) =  8.d-13
     !!!!!!!!!!!!!!!!!!!

            dmdt_gro(ibin,igroup) = ZERO             
            dmdt_grow(ibin) = ZERO
            dmdt_evap(ibin) = ZERO

            do i=1,NGAS
             dmdt_gro(ibin,igroup) = dmdt_gro(ibin,igroup) + xx(i)
            enddo

            do ielemg = 2,nelemg(igroup)
             iesub = ielemg - 1 + iepart

             dmdte_gro(ibin,iesub) = ZERO
             if( itype(iesub) .eq. I_GROWCORE ) then

              igas = igrowgas(iesub)       ! growcore gas

ccc Test - Hold one component constant
c!         !(1) Hold growcore constant if it evaps
c!            if( xx(igas) .lt. 0.) then
c!              dmdt_gro(ibin,igroup) = xx(imaingas)
c!              xx(igas) = 0.
c!            endif
c!         !(2) Hold primary constant if it evaps
c!            if( xx(imaingas) .lt. 0.) then
c!              dmdt_gro(ibin,igroup) = xx(igas)
c!              xx(imaingas) = 0.
c!            endif
ccc End Test

              if( xx(igas) .lt. ZERO .and. 
     $            rnormfrac(ibin,iesub) .lt. FIX_COREF ) then

c         Don't let growcore evaporate newly formed particles,
c         reset dmdt to be that of primary when there is no growcore
c         present yet

                dmdt_gro(ibin,igroup) = dmdt_gro(ibin,igroup) 
     $                                     - xx(igas)
              else 
               dmdte_gro(ibin,iesub) = xx(igas)
              endif

              if( xx(imaingas) .lt. ZERO .and. 
     $            rnormfrac(ibin,iesub) .gt. ONE - FIX_COREF ) then

c         Don't let primary evaporate totally growcore particles,
c         reset dmdt to be that of growcore when there is no primary
c         present ---fix for multiple growcores????????????

                dmdt_gro(ibin,igroup) = xx(igas)
              endif

c         Total dmdt's for growth and evap
              if( dmdte_gro(ibin,iesub) .gt. ZERO ) then
               dmdt_grow(ibin) = dmdt_grow(ibin) 
     $                            + dmdte_gro(ibin,iesub)
              else
               dmdt_evap(ibin) = dmdt_evap(ibin) 
     $                            + dmdte_gro(ibin,iesub)
              endif

             endif !GROWCORE element
            enddo  !sub elements in group
 
c         Include primary's growth rate in appropriate dmdt term
            if( xx(imaingas) .gt. ZERO ) then
             dmdt_grow(ibin) = dmdt_grow(ibin) + xx(imaingas)
            else
             dmdt_evap(ibin) = dmdt_evap(ibin) + xx(imaingas)
            endif

c          if( ixyz.eq.iwa ) then
c           if( itime.eq.1132 .or. itime.eq.1362 .or. 
c    $             itime.eq.2042 .or. itime.eq.2587) then
c            if(ibin.eq.1) 
c    $        write(*,*) '----------------',itime,'--------------'
c            write(*,5) 'dm/dt: ',dmdt(ibin),xx(1),xx(2),ibin,
c    $            dmdt_gro(ibin,igroup),dmdte_gro(ibin,5),
c    $            rnormfrac(ibin,5)
c           endif
c            write(*,7) 'C2H6: ',alpha(ibin,1),beta(ibin,1),
c    $            gamma(ibin,1),dHv(1)
c            write(*,7) 'CH4:  ',alpha(ibin,2),beta(ibin,2),
c    $            gamma(ibin,2),dHv(2)
c          endif

          enddo     ! ibin = 1,NBIN-1
c
c
c  Now calculate condensation/evaporation production and loss rates.
c  Use Piecewise Polynomial Method [Colela and Woodard, J. Comp. Phys.,
c  54, 174-201, 1984]
c
c  First, use cubic fits to estimate concentration values at bin
c  boundaries
c
         do ibin = 2,NBIN-1

          dpc = pc3(ixyz,ibin,iepart) / dm(ibin,igroup)
          dpc1 = pc3(ixyz,ibin+1,iepart) / dm(ibin+1,igroup)
          dpcm1 = pc3(ixyz,ibin-1,iepart) / dm(ibin-1,igroup)
          ratt1 = dm(ibin,igroup) / 
     $      ( dm(ibin-1,igroup) + dm(ibin,igroup) + dm(ibin+1,igroup) )
          ratt2 = ( 2.*dm(ibin-1,igroup) + dm(ibin,igroup) ) /
     $            ( dm(ibin+1,igroup) + dm(ibin,igroup) )
          ratt3 = ( 2.*dm(ibin+1,igroup) + dm(ibin,igroup) ) /
     $            ( dm(ibin-1,igroup) + dm(ibin,igroup) )
          dela(ibin) = ratt1 *
     $             ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )
          delma(ibin) = ZERO

          if( (dpc1-dpc)*(dpc-dpcm1) .gt. ZERO )
     $       delma(ibin) = min( abs(dela(ibin)), 2.*abs(dpc-dpc1),
     $            2.*abs(dpc-dpcm1) ) * abs(dela(ibin))/dela(ibin)

         enddo     ! ibin = 2,NBIN-2

         do ibin = 2,NBIN-2

          dpc = pc3(ixyz,ibin,iepart) / dm(ibin,igroup)
          dpc1 = pc3(ixyz,ibin+1,iepart) / dm(ibin+1,igroup)
          dpcm1 = pc3(ixyz,ibin-1,iepart) / dm(ibin-1,igroup)
c
c
c  <ratx> and <den1> should be calculated in a setup routine
c  (since they are time-independent).
c
          rat1 = dm(ibin,igroup) /
     $            ( dm(ibin,igroup) + dm(ibin+1,igroup) )
          rat2 = 2. * dm(ibin+1,igroup) * dm(ibin,igroup) /
     $           ( dm(ibin,igroup) + dm(ibin+1,igroup) )
          rat3 = ( dm(ibin-1,igroup) + dm(ibin,igroup) ) /
     $           ( 2.*dm(ibin,igroup) + dm(ibin+1,igroup) )
          rat4 = ( dm(ibin+2,igroup) + dm(ibin+1,igroup) ) /
     $           ( 2.*dm(ibin+1,igroup) + dm(ibin,igroup) )
          den1 = dm(ibin-1,igroup) + dm(ibin,igroup) +
     $           dm(ibin+1,igroup) + dm(ibin+2,igroup)
c
c
c  <aju(ibin)> is the estimate for concentration (dn/dm) at bin
c  boundary <ibin>+1/2.
c
          aju(ibin) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $             ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $             dm(ibin,igroup)*rat3*delma(ibin+1) +
     $             dm(ibin+1,igroup)*rat4*delma(ibin) )

         enddo     ! ibin = 2,NBIN-2
c
c
c  Now construct polynomial functions in each bin
c
         do ibin = 3,NBIN-2

          al(ibin) = aju(ibin-1)
          ar(ibin) = aju(ibin)

         enddo
c
c
c  Use linear functions in first two and last two bins
c
         if( NBIN .gt. 1 ) then
         ar(2) = aju(2)
         al(2) = pc3(ixyz,1,iepart)/dm(1,igroup) +
     $           (rmassup(1,igroup)-rmass(1,igroup)) /
     $           (rmass(2,igroup)-rmass(1,igroup)) *
     $           (pc3(ixyz,2,iepart)/dm(2,igroup)-
     $           pc3(ixyz,1,iepart)/dm(1,igroup))
         ar(1) = al(2)
         al(1) = pc3(ixyz,1,iepart)/dm(1,igroup) +
     $           (rmassup(1,igroup)/rmrat(igroup)-rmass(1,igroup)) /
     $           (rmass(2,igroup)-rmass(1,igroup)) *
     $           (pc3(ixyz,2,iepart)/dm(2,igroup)-
     $           pc3(ixyz,1,iepart)/dm(1,igroup))

         al(NBIN-1) = aju(NBIN-2)
         ar(NBIN-1) = pc3(ixyz,NBIN-1,iepart)/dm(NBIN-1,igroup) +
     $           (rmassup(NBIN-1,igroup)-rmass(NBIN-1,igroup))
     $           / (rmass(NBIN,igroup)-rmass(NBIN-1,igroup)) *
     $           (pc3(ixyz,NBIN,iepart)/dm(NBIN,igroup)-
     $           pc3(ixyz,NBIN-1,iepart)/dm(NBIN-1,igroup))
         al(NBIN) = ar(NBIN-1)
         ar(NBIN) = pc3(ixyz,NBIN-1,iepart)/dm(NBIN-1,igroup) +
     $           (rmassup(NBIN,igroup)-rmass(NBIN-1,igroup))
     $           / (rmass(NBIN,igroup)-rmass(NBIN-1,igroup)) *
     $           (pc3(ixyz,NBIN,iepart)/dm(NBIN,igroup)-
     $           pc3(ixyz,NBIN-1,iepart)/dm(NBIN-1,igroup))
         endif
c
c
c  Next, ensure that polynomial functions do not deviate beyond the
c  range [<al(ibin)>,<ar(ibin)>]
c
         do ibin = 1,NBIN

          dpc = pc3(ixyz,ibin,iepart) / dm(ibin,igroup)

          if( (ar(ibin)-dpc)*(dpc-al(ibin)) .le. ZERO ) then
            al(ibin) = dpc
            ar(ibin) = dpc
          endif

          test1 = (ar(ibin)-al(ibin))*(dpc - 0.5*(al(ibin)+ar(ibin)))
          test2 = 1./6.*(ar(ibin)-al(ibin))**2

          if( test1 .gt. test2 ) then
             al(ibin) = 3.*dpc - 2.*ar(ibin)
          elseif( test1 .lt. -test2 ) then
             ar(ibin) = 3.*dpc - 2.*al(ibin)
          endif

         enddo
c
c
c  Lastly, calculate fluxes across each bin boundary.
c
c  Use upwind advection when courant number > 1.
c  
c
         do ibin = 1,NBIN

          dpc = pc3(ixyz,ibin,iepart) / dm(ibin,igroup)
          dela(ibin) = ar(ibin) - al(ibin)
          a6(ibin) = 6. * ( dpc - 0.5*(ar(ibin)+al(ibin)) )

         enddo

         do ibin = 1,NBIN-1

          com2  = ( dm(ibin,igroup) + dm(ibin+1,igroup) ) / 2.
c
c   Calculate growth loss rate     
c
          if( dmdt_grow(ibin) .gt. ZERO )then

           x = dmdt_grow(ibin)*dtime/dm(ibin,igroup) 
           if( x .lt. 1. )then
            growlg(ibin,igroup) = dmdt_grow(ibin)/pc3(ixyz,ibin,iepart)
     $                 * ( ar(ibin) - 0.5*dela(ibin)*x +
     $                 (x/2. - x**2/3.)*a6(ibin) )
           else
            growlg(ibin,igroup) = dmdt_grow(ibin) / com2      ! upwind advection


           endif
          endif
c
c   Calculate evaporation loss rate     
c
          if( dmdt_evap(ibin) .lt. ZERO )then

           x = -dmdt_evap(ibin)*dtime/dm(ibin+1,igroup) 
           if( x .lt. 1. )then
            evaplg(ibin+1,igroup) = -dmdt_evap(ibin)/
     $                 pc3(ixyz,ibin+1,iepart)
     $                 * ( al(ibin+1) + 0.5*dela(ibin+1)*x +
     $                 (x/2. - (x**2)/3.)*a6(ibin+1) )

c          if(ixyz.eq.iwa .and. ibin.ge.iwb-2 .and. ibin.le.NBIN)
c    $      write(*,3) 'Evaplg:',evaplg(ibin+1,igroup),dmdt(ibin),
c    $                pc3(ixyz,ibin+1,iepart),al(ibin+1),dela(ibin+1),
c    $                a6(ibin+1),x,ibin

           else
            evaplg(ibin+1,igroup) = -dmdt_evap(ibin) / com2   ! upwind advection

c          if(ixyz.eq.iwa .and. ibin.ge.iwb-2 .and. ibin.le.NBIN)
c    $      write(*,4) 'Evaplg:',evaplg(ibin+1,igroup),dmdt(ibin),
c    $                com2,ibin

c           if(do_write(11))
c    $       write(*,*) 'Courant Num > 1 (evap)',itime,
c    $               time/60.**2/24.,x,ixyz,ibin,ntsubsteps

           endif


c
c
c  Boundary conditions: for evaporation out of first bin (with cores), 
c  use evaporation rate from second bin.
c  
           if( ibin .eq. 1 .and. ncore(igroup) .gt. 0) then
             evaplg(1,igroup) = -dmdt_evap(1) / dm(1,igroup)
           endif

          endif
c
c
c  Limit growth rates to reasonable values
c
          growlg(ibin,igroup) = min( growlg(ibin,igroup), 1.e10*ONE )
          evaplg(ibin+1,igroup) = min( evaplg(ibin+1,igroup),
     $                                 1.e10*ONE )

cc Test --- no evap / growth
cc
c        if(igroup .ne. 4) then
c         growlg(ibin,igroup) = 0.d0
c         if(ibin .eq. 1 ) evaplg(1,igroup) = 0.d0
c         evaplg(ibin+1,igroup) = 0.d0
c        endif
cc
cc End Test - no evap / growth

         enddo    ! ibin = 1,NBIN-1

        endif     ! (pconmax .gt. FEW_PC)
       endif      ! (imaingas = igrowgas(iepart)) .ne. 0 
      enddo       ! igroup = 1,NGROUP
c
c
c  Return to caller with particle loss rates for growth and evaporation
c  evaluated.
c
      return
      end