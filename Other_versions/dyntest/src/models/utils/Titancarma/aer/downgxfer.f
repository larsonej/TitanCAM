      subroutine downgxfer
c
c
c  @(#) downgxfer.f  Ackerman  Dec-1995
c  This routine calculates particle source terms <rnucpe> due to particle
c  element transfer processes for which the source element number is larger
c  than the target element number.
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
c
c
c   Define formats
c
    1 format(/,'warning in downgxfer: tot evap. core mass > rmass(NBIN)',
     $       /,'  itime, ixyz, i, number created = ',3i8,1pe11.3)
    2 format(/,'warning in downgxfer: tot evap. core mass < rmass(1)',
     $       /,'  itime, ixyz, i, number depleted= ',3i8,1pe11.3)
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter downgxfer'
c
c
c-------------------------------------------------------------------------------
c
c
c  Atmospheric density in cgs units
c
       xyzmet = xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)
       rhoa_cgs = rhoa3(ixyz) /  xyzmet
c
c
c  Calculate nucleation source terms for which the source element
c  number is greater than the target element number
c
c
c
c  Set nucleation production rates to zero to avoid double-application
c  of rates calculated in upgxfer.f
c
      do ielem = 1,NELEM
        do i = 1,NBIN
          rnucpe(i,ielem) = 0.
        enddo
      enddo
c
c
c  Loop over particle elements and bins
c
      do ielem = 1, NELEM
      do ibin = 1, NBIN
c
c
c  Define group & particle # concentration indices for current element
c
       igroup = igelem(ielem)      ! target particle group
       iepart = ienconc(igroup)    ! target particle number concentration element
c
c
c  First calculate production terms due to nucleation <rnucpe>.
c
c
c  Loop over elements that nucleate to element <ielem>.
c
       do jefrom = 1,nnucelem(ielem)

        iefrom = inucelem(jefrom,ielem)    ! source particle element
c
c
c  Only calculate production rates here if <ielem> is less than
c  <iefrom>.  Otherwise, production is calculated in upgxfer.f
c
        if( ielem .lt. iefrom ) then

         igfrom = igelem(iefrom)            ! source particle group
c
c  Loop over bins that nucleate to bin <ibin>.
c
         do jfrom = 1,nnucbin(igfrom,ibin,igroup)

          ifrom = inucbin(jfrom,igfrom,ibin,igroup)    ! bin of source
c
c
c  Bypass calculation if few source particles are present
c
          if( pconmax(ixyz,igfrom) .gt. FEW_PC )then

           if( rnuclg(ifrom,igfrom,igroup) .gt. ZERO )then
c
c
c  First calculate mass associated with the source element <elemass>
c  (this is <rmass> for all source elements except particle number
c  concentration in a multicomponent particle group).
c
            if( ncore(igfrom) .eq. 0 .or.
     $          itype(iefrom) .gt. I_VOLATILE )then
             elemass = rmass(ifrom,igfrom)
            else
              totmass  = pc3(ixyz,ifrom,iefrom) * rmass(ifrom,igfrom)
              rmasscore = pc3(ixyz,ifrom,icorelem(1,igfrom))
              do ic = 2,ncore(igfrom)
                iecore = icorelem(ic,igfrom)
                rmasscore = rmasscore + pc3(ixyz,ifrom,iecore)
              enddo
              fracmass = 1. - rmasscore/totmass
              fracmass = max( ZERO, min( ONE, fracmass ) )
              elemass  = fracmass * rmass(ifrom,igfrom)
            endif

            rnucprod = rnuclg(ifrom,igfrom,igroup) *
     $         pc3(ixyz,ifrom,iefrom) * elemass**ipownuc(iefrom,ielem)

            rnucpe(ibin,ielem) = rnucpe(ibin,ielem) + rnucprod
c
c
c  Calculate latent heat associated with nucleation to <ibin,ielem>
c  from <ifrom,iefrom>
c
c-from 2.0- if( if_nuc_lh(iefrom,ielem) ) then
              rlheat = rlheat + rnucprod * rlh_nuc(iefrom,ielem) /
     $                 ( CP * rhoa_cgs ) * elemass
c           endif
           endif  ! (rnuclg > 0.)
          endif   ! (pconmax > FEW_PC)
         enddo    ! (jfrom = 1,nnucbin)
        endif     ! (ielem < iefrom)
       enddo      ! (jefrom = 1,nnucelem)
      enddo       ! (ibin = 1, NBIN)
      enddo       ! (ielem = 1, NELEM)
c
c
c  Return to caller with down-grid production terms evaluated.
c
      return
      end
