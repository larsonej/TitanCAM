      subroutine upgxfer(ibin,ielem)
c
c
c  @(#) upgxfer.f  Ackerman  Dec-1995
c  This routine calculates particle source terms <rnucpe> due to element transfer
c  processes for which the target element number is greater than the source element
c  number.  (Otherwise, the source terms are calculated in downgxfer.f.)
c  The calculation is done for one particle size bin at one spatial grid point per
c  call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c
c  Argument list input:
c    ibin, ielem
c
c  Argument list output:
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Define formats
c
    1 format(a,'Target(elem/group/#conc)',3(i5),
     $       '  Source(index/elem/group)',3(i5),
     $       '  ipow(from/to)',2(i5))
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter upgxfer'
c
c
c  Calculate the atmospheric density in cgs units
c
      xyzmet = xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)
      rhoa_cgs = rhoa3(ixyz) /  xyzmet
c
c
c  Define group & particle # concentration indices for current element
c
      igroup = igelem(ielem)      ! target particle group 
      iepart = ienconc(igroup)	  ! target particle number concentration element
c
c
c  Calculate production terms due to nucleation <rnucpe>.
c
c
c  Loop over elements that nucleate to element <ielem>.
c
      do jefrom = 1,nnucelem(ielem)

       iefrom = inucelem(jefrom,ielem)    ! source particle element
c
c
c  Only calculate production rates here if <ielem> is greater than
c  <iefrom>.  Otherwise, production is calculated in downgxfer.f
c
       if( ielem .gt. iefrom ) then

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
     $         itype(iefrom) .gt. I_VOLATILE )then
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
     $       pc3(ixyz,ifrom,iefrom) * elemass**ipownuc(iefrom,ielem)

           rnucpe(ibin,ielem) = rnucpe(ibin,ielem) + rnucprod
c
c  Production rate for VOLCORE element of a mixed-phase cloud group
c  from melting/freezing of a single-phase cloud group
c
           if( itype(ielem) .eq. I_VOLATILE  .and. 
     $          is_grp_mixed_phase(igroup)           ) then

               igas = igrowgas(ielem)
               i2e = inuc2elem(1,iefrom)
c
c  For nucleation of mixed phase particles from ice crystals,
c  we need to make sure the ice core mass fraction of the
c  nucleated particle is nearly 1.
c
            if( inucproc(iefrom,ielem) .eq. I_ICEMELT ) then
             rnucpe(ibin,ielem+1) = rnucpe(ibin,ielem+1) +
     $                        rnucprod*elemass * 0.99
c
c  For nucleation of mixed phase particles due to droplet freezing,
c  we need to make sure the ice core mass fraction of the
c  nucleated particle is very small.
c
            elseif( inucproc(iefrom,ielem) .eq. I_DROPFREEZE ) then
             rnucpe(ibin,ielem+1) = rnucpe(ibin,ielem+1) +
     $                       rnucprod*elemass * 0.0001
c
c  For nucleation (vapor deposition) of mixed phase particles at
c  cold temperatures, all mass goes to ice core
c
            elseif( inucproc(iefrom,i2e) .eq. I_VAPORDEP ) then
              if( t3(ixyz) .le. Tfreez(igas) ) then
               rnucpe(ibin,ielem+1) = rnucpe(ibin,ielem+1) +
     $               rnucprod*diffmass(ibin,igroup,ifrom,igfrom)
     $                *ALMOST_ONE
              else
               rnucpe(ibin,ielem+1) = rnucpe(ibin,ielem+1) +
     $               rnucprod*diffmass(ibin,igroup,ifrom,igfrom)
     $                *ALMOST_ZERO
              endif
cc             write(*,*) itime,' ice core nucprod:',ielem+1,
cc   $            ibin,rnucpe(ibin,ielem+1),rnucprod,
cc   $            diffmass(ibin,igroup,ifrom,igfrom)

            endif ! inucproc
           endif ! mixed_phase group
c
c  Calculate latent heat associated with nucleation to <ibin,ielem>
c  from <ifrom,iefrom>
c
c-from 2.0-if( if_nuc_lh(iefrom,ielem) ) then
             rlheat = rlheat + rnucprod * rlh_nuc(iefrom,ielem) /
     $                ( CP * rhoa_cgs ) * elemass


c          endif

          endif  ! (rnuclg > 0.)
         endif   ! (pconmax > FEW_PC)
        enddo    ! (jfrom = 1,nnucbin)
       endif     ! (ielem > iefrom)
      enddo      ! (jefrom = 1,nnucelem)
c
c
c  Return to caller with nucleation production terms evaluated.
c
      return
      end
