      subroutine meltmixedl
c
c
c  @(#) meltmixedl.f  Jensen  Jan-2000
c  This routine evaluates particle loss rates due to nucleation <rnuclg>:
c  total melting of mixed particle volatile core.
c  
c  The loss rates for all particle elements in a particle group are equal.
c
c  Calculations are done at one spatial point per call.
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
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter meltmixedl'
c
c
c  Loop over particle groups.
c
      do igroup = 1,NGROUP

       if( is_grp_mixed(igroup) ) then

        iepart = ienconc( igroup )            ! particle number density element
c
c
c  Calculate nucleation loss rates.
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
c  Only compute nucleation rate for mixed freezing
c
         if( inucproc(iepart,ienucto) .eq. I_MIXEDMELT ) then
c
c
c  Loop over particle bins.  Loop from largest to smallest for 
c  evaluation of index of smallest bin nucleated during time step <inucstep>.
c
          do ibin =NBIN,1,-1
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
           if( pconmax(ixyz,igroup) .gt. FEW_PC )then
c
c
c  Set <rnuclg> to 1.e3 if T > 0C and the total core mass fraction = 0
c  (meaning the mixed particle is entirely melted)
c
            if( t3(ixyz) .gt. T0 ) then

             rmass_core = 0.
             do jcore = 1, ncore(igroup)
               iecore = icorelem(jcore,igroup)
               if( itype(iecore) .eq. I_VOLCORE )
     $           rmass_core = rmass_core + pc3(ixyz,ibin,iecore)
             enddo

             core_mass_frac = rmass_core /
     $             ( pc3(ixyz,ibin,iepart) *rmass(ibin,igroup) )

             if( core_mass_frac .lt. 1.e-4 ) then

              rnuclg(ibin,igroup,ignucto) = 1.e3

              inucstep(igroup) = ibin

             endif

            endif

           endif   ! pconmax(ixyz,igroup) .gt. FEW_PC
          enddo      ! ibin = 1,NBIN
         endif       ! inucproc(iepart,ienucto) .eq. I_DROPFREEZE
        enddo       ! inuc = 1,nnuc2elem(iepart)
       endif        ! is_grp_mixed(igroup)
      enddo         ! igroup = 1,NGROUP
c
c
c  Return to caller with particle loss rates due to nucleation evaluated.
c
      return
      end
