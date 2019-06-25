      subroutine evapp
c
c
c  @(#) evapp.f  Ackerman  Dec-1995
c  This routine calculates particle source terms due to evaporation <evappe>.
c
c  Distinct evaporation of cores has not been treated.
c
c  When the core second moment is used, it is assumed that the probability
c  distribution function of core mass is a log-normal in mass space muliplied
c  by mass raised to the -3/2 power (which implies that the average core
c  mass from the probability distribution function is the same as the
c  average core mass).
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c  Argument list input:
c    
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
      dimension prob(NBIN), evcore(NELEM)

      logical evap_mass, evap_mom, evap_mono, off_cn_grid,
     $        core_too_small, core_too_big
c
c
c  <SIG_MONO> is criterion for a monodisperse core mass distribution
c 
      parameter( SIG_MONO = ALMOST_ZERO )
c
c
c  <CMPDF_EXP> is exponent of core mass in the probability distribution function
c  
      parameter( CMPDF_EXP = -1.5 )
c
c
c   Define formats
c
    1 format(/,'warning in evapp: tot evap. core mass > rmass(NBIN)',
     $       /,'  itime, ixyz, i, number created = ',3i8,1pe11.3)
    2 format(/,'warning in evapp: tot evap. core mass < rmass(1)',
     $       /,'  itime, ixyz, i, number depleted= ',3i8,1pe11.3)
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter evapp'
c
c
c-------------------------------------------------------------------------------
c
c
c  Set evaporation production rates to zero to avoid double-application
c  of rates calculated in growp.f
c
      do ielem = 1,NELEM
        do i = 1,NBIN
          evappe(i,ielem) = 0.
        enddo
      enddo
c
c
c  Calculate source terms due to evaporation
c
c
c
c  <ig> is source group (from which evaporation is being treated).
c
      do ig = 1,NGROUP
 
       iepart = ienconc(ig)           ! source number concentration element
       iecor1 = icorelem(1,ig)        ! element of first core mass in group

       if( iecor1 .ne. 0 .and. igrowgas(iepart) .ne. 0 ) then
       if( itype(iecor1) .ne. I_VOLCORE ) then

        ieto = ievp2elem(iecor1)      ! target number concentration element
                                      ! for total evaporation
        igto = igelem(ieto)           ! target group for total evaporation
 
        do i = 1,NBIN                ! source bin for evaporation
c
c
c  First calculate temporary evaporation source for droplets <evdrop> in
c  bin <i-1> assuming no total evaporation.
c
          evdrop = pc3(ixyz,i,iepart)*evaplg(i,ig)

          if( evdrop .gt. 0. )then

           totevap(i,ig) = .false.
c
c
c  Calculate total and average particle core mass <coretot,coreavg> and
c  core mass fraction <cmf> by summing over all cores in group.
c
           coretot = pc3(ixyz,i,iecor1)
           do ic = 2,ncore(ig)
             iecore = icorelem(ic,ig)
             coretot = coretot + pc3(ixyz,i,iecore)
           enddo
           coreavg = coretot / pc3(ixyz,i,iepart) 
           coreavg = min( rmass(i,ig), coreavg )
           cmf(i,ig) = coreavg / rmass(i,ig)
c
c  Calculate <evcore>, the amount of the source term by number <evdrop>
c  associated with the secondary core <ic>.
c
           do ic = 2,ncore(ig)
             iecore = icorelem(ic,ig)
             evcore(ic) = evdrop * pc3(ixyz,i,iecore) / coretot  
           enddo
c
c
c  Want total evaporation when droplets will be created with core
c  mass fracion > 1 OR evaporating droplets are in bin 1.
c   
           evap_mass = ( rmrat(ig)*cmf(i,ig) .ge. ONE )
     $            .or. ( i .eq. 1 )
c
c
c  Test whether or not average core mass (which is totally evaporating)
c  falls within target CN grid.
c
           core_too_small = coreavg .lt. rmass(1,igto)
           core_too_big   = coreavg .ge. rmass(NBIN,igto)
           off_cn_grid    = core_too_small .or. core_too_big
c
c
c  If evaporating core falls within target CN grid and core second moment is
c  used, calculate average core second moment <coremom>, second moment
c  fraction <smf>, square of the logarithm of the geometric standard deviation
c  of the assumed core mass distribution <coresig>, and index of the target
c  CN bin of the average core mass <iavg>.
c
           if( sec_mom(ig) .and. .not. off_cn_grid )then

             coremom = pc3(ixyz,i,imomelem(ig)) / pc3(ixyz,i,iepart)
  
             iavg = log( coreavg / rmassmin(igto) ) /
     $              log( rmrat(igto) ) + 2
             iavg = min( iavg, NBIN )
       
             smf = coremom / rmass(i,ig)**2
             coresig = log( smf / cmf(i,ig)**2 )
             coresig = max( SIG_MONO, coresig )
c
c
c  Want monodisperse total evaporation when [ droplets will be created with
c  core second moment fraction > 1 OR core mass fraction > 1 ] AND
c  [ [ core mass distribution is monodisperse ] OR
c    [ core is smaller than the smallest CN recently nucleated ] ]
c
             evap_mom  = rmrat(ig)**2 * smf .ge. ONE

             evap_mono = ( evap_mass .or. evap_mom ) .and.
     $                 ( ( coresig .le. SIG_MONO ) .or.
     $                   ( iavg .le. inucmin(igto) ) )

           else
c
c
c  Either there is no core second moment or the average core mass does not
c  fall within target CN grid: monodisperse core distribution is implied.
c
             evap_mom  = .false. 
             evap_mono = evap_mass

           endif  ! sec_mom(ig)
c 
c
c  Treat total evaporation from a monodisperse core mass distribution.
c  
           if( evap_mono )then

            totevap(i,ig) = .true.

            if( core_too_big )then
c
c  Put all of core mass into largest CN bin and print warning
c  (conserves mass, number increases).
c
             factor = coreavg / rmass(NBIN,igto)

             created = ( evdrop * ( factor - ONE ) * dtime )
     $          / ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )

c            if( created .gt. FEW_PC )then
c              write(LUNOPRT,1) itime, ixyz, i, created
c            endif
c
c  First the CN number concentration element
c 
             evappe(NBIN,ieto) = evappe(NBIN,ieto) +
     $                                 factor * evdrop

             if(evappe(NBIN,ieto) .lt. 0.)
     $         print *,'ev1',evappe(NBIN,ieto)
c
c  Now the CN cores
c
             do ic = 2,ncore(ig)
               iecore = icorelem(ic,ig)
               ie2cn  = ievp2elem(iecore)

               evappe(NBIN,ie2cn) = evappe(NBIN,ie2cn) +
     $            factor * evcore(ic) * rmass(NBIN,igto)

             if(evappe(NBIN,ie2cn) .lt. 0.)
     $         print *,'ev2',evappe(NBIN,ie2cn)

             enddo

            elseif( core_too_small )then
c
c  Put all of core mass into smallest CN bin and print warning
c  (conserves mass, number decreases).
c
             factor = coreavg / rmass(1,igto)

             depleted = ( evdrop * ( ONE - factor ) * dtime )
     $          / ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )

c            if( depleted .gt. FEW_PC )then
c              write(LUNOPRT,2) itime, ixyz, i, depleted
c            endif
c
c  First the CN number concentration element
c 
             evappe(1,ieto) = evappe(1,ieto) + factor * evdrop

             if(evappe(1,ieto) .lt. 0.)
     $           print *,'ev3',evappe(1,ieto)
c
c  Now the CN cores
c
             do ic = 2,ncore(ig)
               iecore = icorelem(ic,ig)
               ie2cn  = ievp2elem(iecore)

               evappe(1,ie2cn) = evappe(1,ie2cn) +
     $            factor * evcore(ic) * rmass(1,igto)

             if(evappe(1,ie2cn) .lt. 0.)
     $         print *,'ev4',evappe(1,ie2cn)

             enddo
 
            else
c
c  Partition core mass between two CN bins, conserving total core mass
c  and number.  The number will be subdivided into bins <ito> and <ito>-1.
c
             if( sec_mom(ig) )then
               ito = iavg
             else
               ito = log( coreavg / rmassmin(igto) ) /
     $               log( rmrat(igto) ) + 2
               ito = min( ito, NBIN )
             endif

             fracmass = ( rmass(ito,igto) - coreavg ) /
     $                    diffmass(ito,igto,ito-1,igto)
c
c  Partition the source to the CN number concentration element.
c 
             evappe(ito-1,ieto) = evappe(ito-1,ieto) +
     $              evdrop * fracmass

             if(evappe(ito-1,ieto) .lt. 0.)
     $         print *,'ev5',evappe(ito-1,ieto),ito-1,ieto

             evappe(ito,ieto) = evappe(ito,ieto) +
     $              evdrop * ( ONE - fracmass )

             if(evappe(ito,ieto) .lt. 0.) then
c                print *,'ev6',evappe(ito,ieto)
                 evappe(ito,ieto) = 0.
             endif

c
c  Partition mass from secondary cores to CN so that their core mass fraction
c  as target CN preserves the fraction of core mass in the source droplet.
c
             do ic = 2,ncore(ig)
               iecore = icorelem(ic,ig)
               ie2cn  = ievp2elem(iecore)
c
c  <evcore> is the fraction of the source term by number <evdrop> that is
c  from the secondary core <ic>.
c
c              evcore = evdrop * pc3(ixyz,i,iecore) / coretot  

               evappe(ito-1,ie2cn) = evappe(ito-1,ie2cn) +
     $                rmass(ito-1,igto) * evcore(ic) * fracmass

             if(evappe(ito-1,ie2cn) .lt. 0.)
     $         print *,'ev7',evappe(ito-1,ie2cn)


               evappe(ito,ie2cn) = evappe(ito  ,ie2cn) +
     $                rmass(ito,igto) * evcore(ic) * ( ONE - fracmass )

             if(evappe(ito,ie2cn) .lt. 0.)
     $         print *,'ev8',evappe(ito,ie2cn)

             enddo
 
            endif  ! core_too_big etc.

           elseif( evap_mass .or. evap_mom )then
c
c
c  Treat total evaporation from a polydisperse core mass distribution:
c  assume a log-normal CN size distribution and conserve number and mass as
c  described by Turco (NASA Technical Paper 1362).  Don't put anything
c  in CN bins with index less than <inucmin>.
c  
            totevap(i,ig) = .true.
c
c
c  Calculate number <rn_norms,rn_norml> and mass <rm_norms,rm_norml>
c  normalization factors for cores smaller and larger than <rmass(m,igto)>.
c
            rn_norms = 0.
            rn_norml = 0.
            rm_norms = 0.
            rm_norml = 0.
            kount_s = 0
            kount_l = 0

            do ito = inucmin(igto), NBIN

              rmassto = rmass(ito,igto)           
              dmto    = dm(ito,igto)
c
c  <prob> is probability that core mass is in CN bin <ito>.
c
              expon = -log( rmassto/coreavg )**2 / ( 2.*coresig )
              expon = max(-POWMAX, expon)
              prob(ito) = rmassto**CMPDF_EXP * exp( expon )

              if( ito .le. iavg )then
                rn_norms = rn_norms + prob(ito)*dmto
                rm_norms = rm_norms + prob(ito)*dmto*rmassto
                kount_s = kount_s + 1
              else
                rn_norml = rn_norml + prob(ito)*dmto
                rm_norml = rm_norml + prob(ito)*dmto*rmassto
                kount_l = kount_l + 1
              endif

            enddo
c
c
c  Calculate mass weighting factors <weights,weightl> for small and
c  large cores.
c
            if( kount_s .eq. 0 )then
              rm_norml = rm_norml/rn_norml
              weightl = 1.
            elseif( kount_l .eq. 0 )then
              rm_norms = rm_norms/rn_norms
              weightl = 0.
            else
              rm_norms = rm_norms/rn_norms
              rm_norml = rm_norml/rn_norml
              weightl = (coreavg - rm_norms) / (rm_norml - rm_norms)
              if( (ONE + weightl) .eq. ONE )then
                weightl = 0.
              else if( (ONE - weightl) .le. 0. )then
                weightl = ONE
              endif
            endif
            weights = ONE - weightl
c
c
c Renormalize probability distribution function and evaluate the CN
c evaporation source term <evappe>.
c
            do ito = inucmin(igto), NBIN

              if( ito .le. iavg )then
                prob(ito) = prob(ito)*weights/rn_norms
              else
                prob(ito) = prob(ito)*weightl/rn_norml
              endif
c
c  First the CN number concentration element
c 
              evappe(ito,ieto) = evappe(ito,ieto) +
     $              evdrop * prob(ito) * dm(ito,igto)

             if(evappe(ito,ieto) .lt. 0.)
     $         print *,'ev9',evappe(ito,ieto)

c
c  Now the CN core elements
c
              do ic = 2,ncore(ig)
                iecore = icorelem(ic,ig)
                ie2cn  = ievp2elem(iecore)

                evappe(ito,ie2cn) = evappe(ito,ie2cn) +
     $             rmass(ito,igto) * evcore(ic) * 
     $             prob(ito) * dm(ito,igto)

             if(evappe(ito,ie2cn) .lt. 0.)
     $         print *,'ev10',evappe(ito,ie2cn)

              enddo

            enddo
c
c
c  Evaluate evaporation source term <evappe> for droplet number concentration,
c  core mass, and core second moment.
c
           else

            evappe(i-1,iepart) = evdrop
 
             if(evappe(i-1,iepart) .lt. 0.)
     $         print *,'ev11',evappe(i-1,iepart),i-1,iepart

            do iesub = 2,nelemg(ig)
              ielem = iepart + iesub - 1

              evappe(i-1,ielem) = evappe(i-1,ielem) +
     $           pc3(ixyz,i,ielem) * evaplg(i,ig)

             if(evappe(i-1,ielem) .lt. 0.)
     $         print *,'ev12',evappe(i-1,ielem),i-1,ielem

            enddo  ! iesub=1,nelemg
           endif   ! (1-rmrat*cmf) gt 0 else i gt 1

          endif    ! evdrop gt 0
        enddo      ! i=1,NBIN
      
       endif
       elseif( igrowgas(iepart) .ne. 0 )then     ! ieto eq 0
c
c 
c  Treat evaporation of droplets that have no cores (for which there
c  is no evaporation out of bin 1).
c
        do i = 2,NBIN
          if( evaplg(i,ig) .gt. 0. )then

            evappe(i-1,iepart) = evappe(i-1,iepart) +
     $         pc3(ixyz,i,iepart) * evaplg(i,ig)

             if(evappe(i-1,iepart) .lt. 0.)
     $         print *,'ev13',evappe(i-1,iepart)

          endif    ! evaplg gt 0
        enddo      ! i=2,NBIN
       endif       ! igrowgas ne 0
      enddo        ! ig=1,NGROUP
c
c
c  Return to caller with evaporation production terms evaluated.
c
      return
      end
