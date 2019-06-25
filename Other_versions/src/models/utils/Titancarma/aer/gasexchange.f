      subroutine gasexchange
c
c
c  @(#) gexchange.f  Ackerman  Dec-1995
c  This routine calculates the total production of gases due to nucleation,
c  growth, and evaporation <gasprod> [g/x_units/y_units/z_units/s].
c  It also calculates the latent heating rate from a condensing gas
c  <rlheat> [deg_K/s]
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c
      dimension gprod_nuc(NGROUP,NGAS)
      logical if_do_nuc
c
c
c  Define Formats
c
    1 format(a,2(i3),'  gasprod: ',1pe11.4,'  nucl: ',1pe11.4,
     $       '  grow: ',1pe11.4,'  ex: ',1pe11.4)
    2 format(a,i3,' evaplg: ',1pe13.6,'  growlg: ',1pe13.6)
    3 format(3(i3),3x,'evaplg:',3(1pe13.6,2x),'growlg:',1pe13.6)
    4 format(a,2(i3),'  Both:',1pe11.4,'  CH4:',1pe11.4,
     $       '  drop:',1pe11.4,'  mono:',1pe11.4,'  poly:',1pe11.4)
    5 format(a,i4,2x,'gc/vol:',f15.10,2x,'gc/cld:',f15.10,2x,
     $       'cmf:',f15.10,2x,'cloud(pc):',1pe13.6)
    6 format(a,i4,2x,1pe11.4,' evap: ',1pe11.4,' pc(i+1): ',1pe11.4,
     $       ' gain: ',1pe11.4,' grow: ',1pe11.4,' pc(i): ',1pe11.4,
     $       ' diff: ',1pe11.4,' bin ',i3)
    7 format(a,i4,2x,1pe11.4,' evap: ',1pe11.4,' pc(1): ',1pe11.4,
     $       ' cmf: ',1pe11.4,' mass: ',1pe11.4)
    8 format(a,2(i4),2x,'vapor loss: ',1pe15.8,'  cond gain: ',1pe15.8)
    9 format(a,2(i4),2x,'primary cond mass(tot): ',1pe15.8,2x,
     $       'growcore cond mass: ',1pe15.8,2x,'invol core mass: ',
     $       1pe15.8,2x,'diff: ',1pe15.8)
   10 format(a,'EvapTerm: ',1pe13.6,'  GrowTerm: ',1pe13.6,
     $       '  Primary: ',1pe13.6,'  Growcore: ',1pe13.6,0p,
     $       '  ss(pr): ',f9.6,'  ss(gc): ',f9.6,'  Max(pc): ',1pe13.6)
   11 format(a,2(i4),2(1pe13.6,3x))
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
c  Metric scale factor
c
      xyzmet = xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)
      rhoa_cgs = rhoa3(ixyz) /  xyzmet
c
c
c  Initialize local variable for keeping track of gas changes due
c  to nucleation in each particle group.
c
      do igroup = 1,NGROUP
        do igas = 1,NGAS
          gprod_nuc(igroup,igas) = 0.
        enddo
      enddo

c
c
c  First calculate gas loss and latent heat gain rates due to nucleation.
c
      do igroup = 1,NGROUP
 
ccccccccc
ccccccccc
ccccccccc      No nucl with growcores as souce elem
ccccccccc      Will need to fix gprod_nuc and rlheat
ccccccccc
ccccccccc
       do igs = 1,NGAS
 
        igas = inucgas(igs,igroup)      ! condensing gas

        ielem = ienconc(igroup)     ! element of particle number concentration

        if( igas .ne. 0 .and. nnuc2elem(ielem) .gt. 0 )then

         do ienuc2 = 1,NELEM

          ig2 = igelem( ienuc2 )    ! target particle group

ccc Test
c         if( (ienuc2.eq.6 .and. igas.eq.1) .or.
c    $        (ienuc2.eq.3  .and. igas.eq.2)     ) then
c              if_do_nuc = .false.
c         else 
c              if_do_nuc = if_nuc(ielem,ienuc2)
c         endif

          if( if_nuc(ielem,ienuc2) ) then      !uncomment one
cc        if( if_do_nuc ) then                 !of these
ccc End Test
 
           do i = 1,NBIN
            
            i2 = inuc2bin(i,igroup,ig2)            ! target bin

            gprod_nuc(igroup,igas) = gprod_nuc(igroup,igas) -
     $          pc3(ixyz,i,ielem) * rnuclg(i,igroup,ig2) *
     $          diffmass(i2,ig2,i,igroup)

           enddo    
          
cc         if( gprod_nuc(igroup,igas) .ne. 0. ) 
cc   $       write(*,*) ixyz,igroup,igas,gprod_nuc(igroup,igas)
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

           rlheat = rlheat - rlh * gprod_nuc(igroup,igas) / xyzmet /
     $              ( CP * rhoa_cgs )

          endif
         enddo     ! ienuc2 = 1,NELEM
        endif      ! (igas = inucgas(ielem) .ne. 0
       enddo       ! igs = 1,NGAS
c
c
c  Next calculate gas lost/gained due to and heat gained/lost from 
c  growth/evaporation.  Do this only for the main cloud gas as 
c  production term for growcore gases was calculated elsewhere.
c
         igas = igrowgas(ielem)     ! condensing gas (main gas)

         if( igas .ne. 0 )then

          do i = 1,NBIN-1
c
c
c  Calculate <gasgain>, mass concentration of gas gained due to evaporation
c  from each droplet in bin <i+1>.  First check for total evaporation.
c
            if( totevap(i+1,igroup) )then
              gasgain = ( ONE - cmf(i+1,igroup) )*rmass(i+1,igroup)

c             if( pc3(ixyz,i+1,5).ge.pc3(ixyz,i+1,ielem)*gasgain )
c             if( volpart.ne.pc3(ixyz,i+1,ielem)*gasgain )
c    $         write(*,*) ixyz,i+1,pc3(ixyz,i+1,5)/volpart

            else
              gasgain = diffmass(i+1,igroup,i,igroup)
            endif

c      <gasgain> is total mass of gas gained through evap
c      will separate into amounts for each gas later

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

c
c
c         ------------some write statements------------------
c         if( do_write(4) ) then
c          if(igroup .eq. 2) then
c           if(ixyz .ge. iwa-1 .and. ixyz .le. iwa+1) then
c           if(ixyz .eq. iwa) then
c            write(*,1) 'gasexchange',ixyz,igas,gasprod(igas),
c    $               gprod_nuc(igroup,igas),gprod_grow(igroup,igas)
c            ielem=2
c            if(igas.eq.1) then
c             do i=1,NBIN-1
c              write(*,3) ixyz,igas,i,evaplg(i+1,igroup),
c    $                    pc3(ixyz,i+1,ielem),gasgain,
c    $                    growlg(i  ,igroup)
c             enddo
c            endif!alt
c           endif !alt
c          endif  !group
c         endif   !write
c         ---------------------------------------------------

        enddo

      enddo
c
c
c  Return to caller with <gasprod> evaluated.
c
      return
      end
