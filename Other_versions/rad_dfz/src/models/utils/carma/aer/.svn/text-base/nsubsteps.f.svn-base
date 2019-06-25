      subroutine nsubsteps
c
c
c  @(#) nsubsteps.f  Jensen  Apr-2000
c  This routine calculates the number of sub-timesteps <ntsubsteps>
c  for the current model spatial point.
c
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c  Argument list input:
c    None.
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
c  Local declarations
c
      dimension ibin_small(NGROUP)
c
c---------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter nsubsteps'
c
c
c  Set default values
c
      ntsubsteps = minsubsteps
c
c
c  If the supersaturation crosses 0 or significant nucleation will occur,
c  then use <maxsubsteps>
c
      do ig = 1, NGROUP

        iepart = ienconc(ig)     ! element of particle number concentration
        igas = igrowgas(iepart)  ! condensing gas

        if( itype(iepart) .eq. I_VOLATILE ) then

         if( is_grp_ice(ig) )then
          ss = supsati3(ixyz,igas)
          ssold = supsatiold(ixyz,igas)
         else
          ss = supsatl3(ixyz,igas)
          ssold = supsatlold(ixyz,igas)
         endif
 
         if( ( ssold .lt. 0. .and. ss .ge. 0. ) .or.
     $       ( ssold .ge. 0. .and. ss .lt. 0. ) ) then
           ntsubsteps = maxsubsteps
         endif

        endif
      enddo
c
c
c  Find the bin number of the smallest particle bin that
c  contains a significant number of particles.
c  Also check for significant activation of water droplets.
c
      if( ntsubsteps .lt. maxsubsteps )then

        do ig = 1, NGROUP

          if( pconmax(ixyz,ig) .gt. conmax ) then

            ibin_small(ig) = NBIN

            iepart = ienconc(ig)     ! element of particle number concentration
            igas = inucgas(iepart)   ! condensing gas

            if( itype(iepart) .eq. I_INVOLATILE ) then

              ss = max( supsatl3(ixyz,igas), supsatlold(ixyz,igas) )

              do inuc = 1,nnuc2elem(iepart)
                ienucto = inuc2elem(inuc,iepart)

                if( inucproc(iepart,ienucto) .eq. I_DROPACT ) then
                  do ibin = 1, NBIN
                    if( pc3(ixyz,ibin,iepart) .gt. conmax*10 .and.
     $                 ss .gt. scrit(iz,ibin,ig) )then
                      ntsubsteps = maxsubsteps
                    endif
                  enddo
                endif

              enddo
  
            elseif( itype(iepart) .eq. I_VOLATILE ) then

              do ibin = NBIN-1, 1, -1
                if( pc3(ixyz,ibin,iepart) .gt. conmax/10 )then
                  ibin_small(ig) = ibin
                endif
              enddo

            endif
          endif
        enddo
      endif
c
c
c  Calculate the growth rate of a particle with the mode radius for
c  each volatile group.  The maximum time-step to use is then the
c  mass growth rate divided by the mass bin width / 2.
c
      if( ntsubsteps .lt. maxsubsteps )then

        dt_adv = dtime_save
        do ig = 1, NGROUP
        
          iepart = ienconc(ig)     ! element of particle number concentration
          igas = igrowgas(iepart)  ! condensing gas

          if( pconmax(ixyz,ig) .gt. conmax ) then

            if( itype(iepart) .eq. I_VOLATILE ) then

              if( is_grp_ice(ig) )then
                ss = supsati3(ixyz,igas)
                pvap = pvapi3(ixyz,igas)
              else
                ss = supsatl3(ixyz,igas)
                pvap = pvapl3(ixyz,igas)
              endif

              g0 = gro(iz,ibin_small(ig),ig)
              g1 = gro1(iz,ibin_small(ig),ig)
              dmdt = abs( pvap * ss * g0 / ( 1. + g0*g1*pvap ) )
              dt_adv = min( dt_adv, 0.5*dm(ibin_small(ig),ig)/dmdt )
  
            endif
          endif
        enddo
  
        ntsubsteps = min( maxsubsteps, nint( dtime_save/dt_adv ) )
        ntsubsteps = max( minsubsteps, ntsubsteps )

      endif
c
c
c  If the ice supersaturation is large enough for homogeneous freezing
c  of sulfate aerosols, then use maximum number of substeps
c
c      if( supsati3(ixyz,1) .gt. 0.4 .and.
c     $    t3(ixyz) .lt. 233.16 ) ntsubsteps = maxsubsteps
c
c     ntsubsteps = minsubsteps    ! hack
c     if(itime.eq.1)print*,'hack in nsubsteps().'
c
c  Return to caller with number of sub-timesteps evaluated.
c
      return
      end
