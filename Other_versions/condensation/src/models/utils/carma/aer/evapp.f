      subroutine evapp
c
c
c  @(#) evapp.f  Ackerman  Aug-2001
c  This routine calculates particle source terms due to evaporation <evappe>.
c
c  Distinct evaporation of cores has not been treated.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c  Local declarations
c
      logical evap_total
      character*(4) totevp_comp
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
c  Define criterion for monodisperse core mass distributions
c
      sig_mono = sqrt( ALMOST_ZERO )
c
c
c  Loop over source groups (from which evaporation is being treated)
c
      do ig = 1, NGROUP
c
c  source number concentration element
c 
        ip = ienconc(ig)
c
c  No evaporation unless particles are volatile
c  
        if( itype(ip) .eq. I_VOLATILE )then
c
c  Element of first core mass in group
c
          ic1 = icorelem(1,ig)
          if( itype(ic1) .ne. I_COREMASS ) then
           do ic = 2,ncore(ig)
            if( itype(icorelem(ic,ig)) .eq. I_COREMASS )
     $        ic1 = icorelem(ic,ig)
           enddo
          endif
c
c  Loop over source bins and calculate temporary evaporation source
c  for droplets in next smaller bin assuming no total evaporation <evdrop>
c
          do ibin = 1, NBIN
            totevap(ibin,ig) = .false.
            evdrop = pc3(ixyz,ibin,ip)*evaplg(ibin,ig)
c
c  Check for evaporation of a sufficient number of droplets
c  
            if( evdrop .gt. 0. .and.
     $          pc3(ixyz,ibin,ip) .gt. SMALL_PC )then
c
c  No cores: transfer droplets within group
c
              if( ic1 .eq. 0 )then
                call evap_drops(ibin,ig,ip)
              else
c
c  First core is not involatile (therefore none are)
c  -- this is a hack until enforced/checked in setupbins() --
c  transfer droplets within group
c
                if( itype(ic1) .ne. I_COREMASS .and.
     $              itype(ic1) .ne. I_CCN )then
                  call evap_drops(ibin,ig,ip)
                else
c
c  Have cores: calculate <evcore> the amount of the source term 
c  by number <evdrop> associated with total evaporation of secondary cores
c EJL 2-14-13 Added code below to replace commented out code
c
                  ieto = ievp2elem(ic1)

                  if( itype(ieto) .eq. I_COREMASS .or.
     $                itype(ieto) .eq. I_GROWCORE  )  then
                    do ic = 2, ncore(ig)
                     iecore = icorelem(ic,ig)
                     ie2cn = ievp2elem(iecore)
                     if( itype(ie2cn) .eq. I_INVOLATILE .or.
     $                   itype(ie2cn) .eq. I_VOLATILE  ) then
                       ieto = ievp2elem(iecore)
                     endif
                    enddo
                  endif

                  igto = igelem(ieto)

                  coretot = ZERO
                  do ic = 1, ncore(ig)
                    iecore = icorelem(ic,ig)
                    if( ievp2elem(iecore) .ne. 0 ) 
     $               coretot = coretot + pc3(ixyz,ibin,iecore)
                  enddo

                  do ic = 1, ncore(ig)
                    iecore = icorelem(ic,ig)
                    if( ievp2elem(iecore) .ne. ieto ) 
     $                evcore(ic) = evdrop*pc3(ixyz,ibin,iecore)/coretot  
                  enddo 

!                  coretot = pc3(ixyz,ibin,ic1)
!                  do ic = 2, ncore(ig)
!                    iecore = icorelem(ic,ig)
!                    if( itype(iecore) .eq. I_COREMASS .or.
!     $                  itype(iecore) .eq. I_CCN )then
!                      coretot = coretot + pc3(ixyz,ibin,iecore)
!                    endif
!                  enddo
!                  do ic = 2, ncore(ig)
!                    iecore = icorelem(ic,ig)
!                    if( itype(iecore) .eq. I_COREMASS .or.
!     $                  itype(iecore) .eq. I_CCN )then
!                      evcore(ic) = evdrop*pc3(ixyz,ibin,iecore)/coretot  
!                    endif
!                  enddo 
c
c  Calculate average particle core mass and fraction
c
                  coreavg = coretot / pc3(ixyz,ibin,ip) 
                  coreavg = min( rmass(ibin,ig), coreavg )
                  cmf(ibin,ig) = coreavg / rmass(ibin,ig)
c                 cmf(ibin,ig) = max( 0., min( ONE, cmf(ibin,ig) ) )
c
c  Get target number concentration element and group for total evaporation
c  and evaluate logical flags regarding position on CN bin and index of
c  target CN bin
c

                  too_small = coreavg .lt. rmass(1,igto)
                  too_big   = coreavg .gt. rmass(NBIN,igto)

                  if( .not. (too_small .or. too_big) )then
                    iavg = log( coreavg / rmassmin(igto) ) /
     $                     log( rmrat(igto) ) + 2
                    iavg = min( iavg, NBIN )
                  endif
c
c
c  Only consider size of evaporating cores relative to nuc_small
c  when treating core second moment for this particle group
c
                  if( if_sec_mom(ig) )then
                    nuc_small = coreavg .lt. rmass(inucmin(igto),igto)
                  else
                    nuc_small = .false.
                  endif
c
c  Want total evaporation when 
c   cores smaller than smallest nucleated 
c   OR evaporating droplets are in bin 1
c   OR droplets will be created with core mass fraction > 1
c 
                  evap_total = nuc_small .or. ibin .eq. 1 .or.
     $                rmrat(ig)*cmf(ibin,ig) .gt. ONE
c
c  No core second moment: evaporate to monodisperse CN cores or within group.
c
                  if( .not. if_sec_mom(ig) )then
 
                    if( evap_total )then
                      call evap_mono(ibin,ig,iavg,ieto,igto)
                    else
                      call evap_drops(ibin,ig,ip)
                    endif
c
c  Have core second moments: evaporate to mono- or polydisperse CN cores
c  or within group.  First calculate average core second moment <coremom>, 
c  second moment fraction <smf>, and square of the logarithm of the geometric
c  standard deviation of the assumed core mass distribution <coresig>.
c
                  else
 
                    coremom = pc3(ixyz,ibin,imomelem(ig)) / 
     $                        pc3(ixyz,ibin,ip)
                    smf = coremom / rmass(ibin,ig)**2
                    coresig = log( smf / cmf(ibin,ig)**2 )
c
c  Want total evaporation for above reasons 
c   OR droplets will be created with core moment fraction > 1
c 
                    evap_total = evap_total .or.
     $                  rmrat(ig)**2*smf .gt. ONE 

                    if( evap_total )then
c
c  Want monodisperse total evaporation when
c   cores smaller than smallest nucleated 
c   OR evaporating core distribution is narrow
c  Otherwise want polydisperse total evaporation 
c 
                      if( nuc_small .or. coresig .le. sig_mono .or.
     $                    ievp2elem(imomelem(ig)) .ne. 0      )then
                        call evap_mono(ibin,ig,iavg,ieto,igto)
                      else
                        call evap_poly(ibin,ig,iavg,ieto,igto)
                      endif
c
c  Droplet evaporation within group
c
                    else
                      call evap_drops(ibin,ig,ip)
                    endif

                  endif      ! if_sec_mom(ig)
                endif        ! itype(ic1)
              endif          ! ic1=0 
            endif            ! evaplg > 0
          enddo              ! ibin=1,NBIN
        endif                ! volatile particles
      enddo                  ! ig=1,NGROUP
c
c
c  Return to caller with evaporation production terms evaluated.
c
      return
      end
