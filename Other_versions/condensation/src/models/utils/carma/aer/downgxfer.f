      subroutine downgxfer
c
c
c  @(#) downgxfer.f  Ackerman  Dec-1995
c  This routine calculates particle source terms <rnucpe> due to particle
c  element transfer processes for which the source element number is larger
c  than the target element number.
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
c
c  <ipow> is the power to which the source particle mass must be taken
c  to match the type of the target element.  This ugliness could be
c  handled much more slickly in setupnuc()
c
         if( itype(iefrom) .eq. I_INVOLATILE .or.
     $       itype(iefrom) .eq. I_VOLATILE )then
           ipow_from = 0
         elseif ( itype(iefrom) .eq. I_COREMASS .or.
     $       itype(iefrom) .eq. I_VOLCORE )then
           ipow_from = 1
         else
           ipow_from = 2
         endif

         if( itype(ielem) .eq. I_INVOLATILE .or.
     $       itype(ielem) .eq. I_VOLATILE )then
           ipow_to = 0
         elseif ( itype(ielem) .eq. I_COREMASS .or.
     $       itype(ielem) .eq. I_VOLCORE )then
           ipow_to = 1
         else
           ipow_to = 2
         endif

         ipow = ipow_to - ipow_from
c
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

           if( rnuclg(ifrom,igfrom,igroup) .gt. 0. )then
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
              elemass  = fracmass * rmass(ifrom,igfrom)
            endif

            rnucprod = rnuclg(ifrom,igfrom,igroup) *
     $        pc3(ixyz,ifrom,iefrom) * elemass**ipownuc(iefrom,ielem)

            rnucpe(ibin,ielem) = rnucpe(ibin,ielem) + rnucprod
c
c
c  Calculate latent heat associated with nucleation to <ibin,ielem>
c  from <ifrom,iefrom>
c
            rlheat = rlheat + rnucprod * rlh_nuc(iefrom,ielem) /
     $               ( CP * rhoa_cgs ) * elemass

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
