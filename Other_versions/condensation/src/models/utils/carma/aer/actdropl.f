      subroutine actdropl
c
c
c  @(#) actdropl.f  Ackerman  Dec-1995
c  This routine evaluates particle loss rates due to nucleation <rnuclg>:
c  droplet activation only.
c  
c  The loss rates for all particle elements in a particle group are equal.
c
c  To avoid nucleation into an evaporating bin, this subroutine must
c  be called after growp, which evaluates evaporation loss rates <evaplg>.
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
      logical evapfrom_nucto
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter actdropl'
c
c
c  This calculation is only necessary for temperatures greater
c  than -40C.
c
c
c
c  Loop over particle groups.
c
       do igroup = 1,NGROUP
        do igs = 1,NGAS
c
c
c  Bypass calculation if few particles are present
c
        if( pconmax(ixyz,igroup) .gt. FEW_PC )then

          igas = inucgas(igs,igroup)                ! condensing gas
          iepart = ienconc( igroup )            ! particle number density element
 
          if( igas .ne. 0 )then
c
c
c  Calculate nucleation loss rates.  Do not allow nucleation into
c  an evaporating bin.
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
c  Only compute nucleation rate for droplet activation
c
            if( inucproc(iepart,ienucto) .eq. I_DROPACT ) then
c
c
c  Loop over particle bins.  Loop from largest to smallest for 
c  evaluation of index of smallest bin nucleated during time step <inucstep>.
c
             do ibin = NBIN, 1, -1
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
c  Set <evapfrom_nucto> to .true. when target droplets are evaporating
c  
             if( inucto .ne. 0 )then
               evapfrom_nucto = evaplg(inucto,ignucto) .gt. 0.
             else
               evapfrom_nucto = .false.
             endif
            
             if( (supsatl3(ixyz,igas) .gt. scrit(iz,ibin,igroup)) .and.
     $           (.not. evapfrom_nucto) .and.
     $           (pc3(ixyz,ibin,iepart) .gt. SMALL_PC) )then

              rnuclg(ibin,igroup,ignucto) = 1.e3

              if( ibin .lt. inucstep(igroup) )then
                inucstep(igroup) = ibin
              endif

             endif

            enddo   ! ibin = 1,NBIN
           endif    ! inucproc(iepart,ienucto) .eq. I_DROPACT
          enddo     ! inuc = 1,nnuc2elem(iepart)
         endif      ! (igas = inucgas(igroup)) .ne. 0 
        endif       ! pconmax(ixyz,igroup) .gt. FEW_PC
        enddo       ! igs
       enddo        ! igroup = 1,NGROUP
c
c
c  Return to caller with particle loss rates due to nucleation evaluated.
c
      return
      end