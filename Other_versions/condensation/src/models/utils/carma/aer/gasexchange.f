      subroutine gasexchange
c
c
c  @(#) gexchange.f  Ackerman  Dec-1995
c  This routine calculates the total production of gases due to nucleation,
c  growth, and evaporation <gasprod> [g/x_units/y_units/z_units/s].
c  It also calculates the latent heating rate from a condensing gas
c  <rlheat> [deg_K/s]
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c  
c
c  Declare local variables
c
      dimension gprod_nuc(NGROUP,NGAS) !, gprod_grow(NGROUP,NGAS)
      logical if_do_nuc
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter gasexchange'
c
c
c-------------------------------------------------------------------------------
c
c
c  Initialize local variables for keeping track of gas changes due
c  to nucleation and growth in each particle group.
c
      do igroup = 1,NGROUP
        do igas = 1,NGAS
          gprod_nuc(igroup,igas) = 0.
          gprod_grow(igroup,igas) = 0.
        enddo
      enddo
c
c  First calculate gas loss and latent heat gain rates due to nucleation.
c
      do igroup = 1,NGROUP
       do igs =1,NGAS
 
        igas = inucgas(igs,igroup)      ! condensing gas
        ielem = ienconc(igroup)     ! element of particle number concentration

        if( igas .ne. 0 .and. nnuc2elem(ielem) .gt. 0 )then

         do ienuc2 = 1,NELEM

          ig2 = igelem( ienuc2 )    ! target particle group

          if( if_nuc(ielem,ienuc2) ) then
 
           do i = 1,NBIN
            
            i2 = inuc2bin(i,igroup,ig2)            ! target bin

            gprod_nuc(igroup,igas) = gprod_nuc(igroup,igas) -
     $          pc3(ixyz,i,ielem) * rnuclg(i,igroup,ig2) *
     $          diffmass(i2,ig2,i,igroup)

           enddo    
c
c  Latent heating rate from condensing gas: <rlh> is latent heat of evaporation 
c  ( + fusion, for ice deposition ) [erg/g]
c
           if( inucproc(ielem,ienuc2) .eq. I_DROPACT )then
             rlh = rlhe(iz,igas)
           elseif( inucproc(ielem,ienuc2) .eq. I_AERFREEZE )then
             rlh = rlhe(iz,igas) + rlhm(iz,igas)
           endif
           if( is_grp_ice(ig2) ) then
             rlh = rlhe(iz,igas) + rlhm(iz,igas)
           else
             rlh = rlhe(iz,igas)
           endif

           rlheat = rlheat - rlh * gprod_nuc(igroup,igas) /
     $              ( CP * rhoa3(ixyz) )

          endif
         enddo     ! ienuc2 = 1,NELEM
        endif      ! (igas = inucgas(ielem) .ne. 0
       enddo !igs
c
c
c  Next calculate gas lost/gained due to and heat gained/lost from 
c  growth/evaporation.
c
        igas = igrowgas(ielem)     ! condensing gas

        if( igas .ne. 0 )then

          do i = 1,NBIN-1
c
c
c  Calculate <gasgain>, mass concentration of gas gained due to evaporation
c  from each droplet in bin <i+1>.  First check for total evaporation.
c
            if( totevap(i+1,igroup) )then
              gasgain = ( 1. - cmf(i+1,igroup) )*rmass(i+1,igroup)
            else
              gasgain = diffmass(i+1,igroup,i,igroup)
            endif

            gprod_grow(igroup,igas) = gprod_grow(igroup,igas)
     $        + evaplg(i+1,igroup) * pc3(ixyz,i+1,ielem) * 
     $          gasgain
     $        - growlg(i  ,igroup) * pc3(ixyz,i  ,ielem) *
     $          diffmass(i+1,igroup,i,igroup)

          enddo    
c
c
c  Add evaporation out of smallest bin (always total evaporation).
c
          gprod_grow(igroup,igas) = gprod_grow(igroup,igas) +
     $      evaplg(1,igroup) * pc3(ixyz,1,ielem) *
     $      ( 1. - cmf(1,igroup) ) * rmass(1,igroup)
c
c  Now account for the growcore gas 

          do ielemg = 2,nelemg(igroup)
            ieaux = ielem - 1 + ielemg
            iauxgas = igrowgas(ieaux)
            if(iauxgas .ne. 0) then

              gprod_grow(igroup,igas) = gprod_grow(igroup,igas)
     $                                   - gprod_grow(igroup,iauxgas)

            endif
          enddo !other elements in group
c
c
c  Latent heating rate from condensing gas: <rlh> is latent heat of evaporation 
c  ( + fusion, for ice deposition ) [erg/g]
c
          if( is_grp_ice(igroup) )then
            rlh = rlhe(iz,igas) + rlhm(iz,igas)
          else
            rlh = rlhe(iz,igas)
          endif

          rlheat = rlheat - rlh * gprod_grow(igroup,igas) /
     $             ( CP * rhoa3(ixyz) )
c
c EJL 3-14-13 Added below
          do ielemg = 2,nelemg(igroup)
            ieaux = ielem - 1 + ielemg
            iauxgas = igrowgas(ieaux)
            if(iauxgas .ne. 0) then

            if( is_grp_ice(igroup) )then
              rlh = rlhe(iz,iauxgas) + rlhm(iz,iauxgas)
            else
              rlh = rlhe(iz,iauxgas)
            endif

            rlheat = rlheat - rlh * gprod_grow(igroup,iauxgas) /
     $               ( CP * rhoa3(ixyz) )
 
            endif
          enddo !other elements in group

        endif      ! (igas = igrowgas(ielem)) .ne. 0

      enddo        ! igroup=1,NGROUP
c
c
c  Sum up gas production from nucleation and growth terms.
c
      do igroup = 1,NGROUP
        do igas = 1,NGAS
          gasprod(igas) = gasprod(igas) +
     $       gprod_nuc(igroup,igas) + gprod_grow(igroup,igas)
        enddo
      enddo
c
c
c  Return to caller with <gasprod> evaluated.
c
      return
      end
