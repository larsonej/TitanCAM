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
      dimension rmass_mode(NGROUP), ibin_mode(NGROUP),
     $          ibin_small(NGROUP)

      parameter( BAL = 6.1121e3 )           
      parameter( BBL = 18.729 )
      parameter( BCL = 257.87 )
      parameter( BDL = 227.3 )
c
c
c  Define formats
c
    1 format(/,'Timestep     dtime    time     dpcmax   ixyzmax',
     $         '  ielemax  ibinmax    dsmax   ixyzsmax'/)
    2 format(i6,3x,1pe8.2,3x,1pe8.2,3x,1pe9.2,3x,
     $       i3,5x,i4,5x,i4,5x,1pe8.2,4x,i4)
c
c---------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter nsubsteps'
c
c
c  Find the mode mass (volume weighted) of each hydrometeor group
c  and find the mass bin closest to this mass
c
      do ig = 1, NGROUP

        if( pconmax(ixyz,ig) .gt. conmax ) then

          rmass_mode(ig) = 0.
          ibin_mode(ig) = 0
          iepart = ienconc(ig)     ! element of particle number concentration

          if( itype(iepart) .eq. I_VOLATILE ) then

            voltot = 0.
            rntot = 0.
            do ibin = 1, NBIN
              voltot = voltot + pc3(ixyz,ibin,iepart)*vol(ibin,ig)
              rntot = rntot + pc3(ixyz,ibin,iepart)
            enddo
            rmass_mode(ig) = voltot/rntot * rhoelem(iepart)

            ibin_mode(ig) = log( rmass_mode(ig) / rmassmin(ig) ) /
     $                      log( rmrat(ig) ) + 2
            ibin_mode(ig) = max( 1, min( ibin_mode(ig), NBIN-1 ) )

          endif

        endif

      enddo
c
c
c  Find the bin number of the smallest particle bin that
c  contains a significant number of particles
c
      do ig = 1, NGROUP

        if( pconmax(ixyz,ig) .gt. conmax ) then

          ibin_small(ig) = NBIN-1
          iepart = ienconc(ig)     ! element of particle number concentration

          if( itype(iepart) .eq. I_VOLATILE ) then

            do ibin = NBIN-1,1,-1
              if( pc3(ixyz,ibin,iepart) .gt. conmax/10. ) then
                ibin_small(ig) = ibin
              endif
            enddo

          endif

        endif

      enddo
c
c
c  Calculate the growth rate of a particle with the mode radius for
c  each volatile group.  The maximum time-step to use is then the
c  mass growth rate divided by the mass bin width / 2.
c
      dt_adv = 1.e10
      do ig = 1, NGROUP
        
        if( pconmax(ixyz,ig) .gt. conmax ) then

          iepart = ienconc(ig)     ! element of particle number concentration
          igas = igrowgas(iepart)      ! condensing gas

          if( itype(iepart) .eq. I_VOLATILE ) then

            if( t3(iz) .le. Tfreez(igas) ) then !is_grp_ice(ig) )then
              ss = supsati3(ixyz,igas)
              pvap = pvapi3(ixyz,igas)
            else
              ss = supsatl3(ixyz,igas)
              pvap = pvapl3(ixyz,igas)
            endif

            g0 = gro(iz,ibin_small(ig)+1,ig)
            g1 = gro1(iz,ibin_small(ig)+1,ig)
            dmdt = abs( pvap * ss * g0 / ( 1. + g0*g1*pvap ) )
            dt_adv = min( dt_adv, 0.5*dm(ibin_small(ig)+1,ig)/dmdt )

          endif

        endif

      enddo
      ntsubsteps = min( maxsubsteps, nint(dtime_save/dt_adv) )
      ntsubsteps = max( minsubsteps, ntsubsteps )
c
c
c  If the supersaturation crosses 0, then use <maxsubsteps>
c
      do ig = 1, NGROUP

       iepart = ienconc(ig)     ! element of particle number concentration
       igas = igrowgas(iepart)      ! condensing gas

       if( itype(iepart) .eq. I_VOLATILE ) then

        if( t3(iz) .le. Tfreez(igas) ) then !is_grp_ice(ig) )then
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
c  If the ice supersaturation is large enough for homogeneous freezing
c  of sulfate aerosols, then use maximum number of substeps
c
c      if( supsati3(ixyz,1) .gt. 0.4 .and.
c     $    t3(ixyz) .lt. 233.16 ) ntsubsteps = maxsubsteps
c
c
c  Return to caller with number of sub-timesteps evaluated.
c
      return
      end
