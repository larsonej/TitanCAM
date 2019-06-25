       subroutine setupnuc
c
c
c  @(#) setupnuc.f  Ackerman  Dec-1995
c  This routine evaluates derived mapping arrays and calculates the critical
c  supersaturation <scrit> used to nucleate dry particles (CN) to droplets.
c
c
c  This routine requires that array <akelvin> is defined.
c  (i.e., setupgkern.f must be called before this)
c
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
      logical bad_grid
c
c
c  Define formats
c
    1 format(a,':  ',12i6)
    2 format(/,a,':  ',i6)
    3 format(a,a)
    4 format(a,':  ',1pe12.3)
    5 format(/,'Particle nucleation mapping arrays (setupnuc):')
    6 format(i4,5x,1p2e11.3)
    7 format(/,'Warning: nucleation cannot occur from group',i3,
     $ '   bin',i3,'   into group',i3,'   (<inuc2bin> is zero)')
    8 format(/,'Critical supersaturations for ',a,//,
     $ '   i        r [cm]     scrit',/)
    9 format(a,':  ',f5.3)
   10 format(i4,':  ',12i6)
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupnuc'
c
c-------------------------------------------------------------------------------
c
c  Bin mapping for nucleation : nucleation would transfer mass from particles
c  in <ifrom,igfrom> into target bin <inuc2bin(ifrom,igfrom,igto)> in group
c  <igto>.  The target bin is the smallest bin in the target size grid with
c  mass exceeding that of nucleated particle.
c  
      do igfrom = 1,NGROUP		! nucleation source group
        do igto = 1,NGROUP    		! nucleation target group

          do ifrom = 1,NBIN		! nucleation source bin

          inuc2bin(ifrom,igfrom,igto) = 0

            do ibto = NBIN,1,-1        ! nucleation target bin

              if( rmass(ibto,igto) .ge. rmass(ifrom,igfrom) )then
                inuc2bin(ifrom,igfrom,igto) = ibto
              endif

            enddo
          enddo
        enddo
      enddo
c
c
c  Mappings for nucleation sources: 
c
c   <nnucelem(ielem)> is the number of particle elements that nucleate to
c    particle element <ielem>.
c
c   <inucelem(jefrom,ielem)> are the particle elements that
c    nucleate to particle element <ielem>, where 
c    jefrom = 1,nnucelem(ielem).
c
c   <nnucbin(igfrom,ibin,igroup)> is the number of particle bins that nucleate
c    to particles in bin <ibin,igroup> from group <igfrom>.
c
c   <inucbin(jfrom,igfrom,ibin,igto)> are the particle bins 
c    that nucleate to particles in bin <ibin,igto>, where
c    jfrom = 1,nnucbin(igfrom,ibin,igto).
c  
c
c  (Moved from {setupaer} since these can be determined from <inuc2elem>
c
c  <if_nuc(iefrom,ieto)> is true if nucleation transfers mass from element
c  <iefrom> to element <ieto>.
c
c  <nnuc2elem(ielem)> is the total number of target elements from which
c  source particles in element <ielem> are nucleated.
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c
c  Set <nnuc2elem> to zero if mass is not transferred by nucleation
c  (in which case <inuc2elem> is automatically set to zero elsewhere).
c
      do iefrom = 1,NELEM
        do ieto = 1,NELEM
          if_nuc(iefrom,ieto) = .false.
        enddo
      enddo
      do iefrom = 1,NELEM
       nnuc2elem(iefrom) = 0
       do i = 1,5  !to match definition if {globaer.h}
        if( inuc2elem(i,iefrom) .ne. 0 ) then
          nnuc2elem(iefrom) = nnuc2elem(iefrom) + 1
          if_nuc(iefrom,inuc2elem(i,iefrom)) = .true.
          write(*,*) 'if_nuc = .true.',iefrom,inuc2elem(i,iefrom)
        endif
       enddo
      enddo
c
c
c  First enumerate and count elements that nucleate to cores.
c
      do igroup = 1,NGROUP

        ip = ienconc(igroup)    ! target particle number concentration element

        do jcore = 1,ncore(igroup)

          iecore = icorelem(jcore,igroup)    ! target core element 
          nnucelem(iecore) = 0

          do iefrom = 1,NELEM

            if( if_nuc(iefrom,iecore) ) then

              nnucelem(iecore) = nnucelem(iecore) + 1
              inucelem(nnucelem(iecore),iecore) = iefrom

            endif
          enddo      ! iefrom=1,NELEM
        enddo        ! jcore=1,ncore
      enddo          ! igroup=1,NGROUP
c
c  
c  Now enumerate and count elements nucleating to particle concentration
c  (itype=I_INVOLATILE and itype=I_VOLATILE) and core second moment
c  (itype=I_CORE2MOM).  Elements with itype = I_VOLATILE are special because all
c  nucleation sources for core elements in same group are also sources
c  for the itype = I_VOLATILE element.
c
      do igroup = 1,NGROUP
      
        ip = ienconc(igroup)    ! target particle number concentration element
        im = imomelem(igroup)   ! target core second moment element
 
        nnucelem(ip) = 0
        if( im .ne. 0 )then
          nnucelem(im) = 0
        endif
 
        do jcore = 1,ncore(igroup)
 
         iecore = icorelem(jcore,igroup)       ! target core mass element

          do jnucelem = 1,nnucelem(iecore)     ! elements nucleating to cores
 
            iefrom = inucelem(jnucelem,iecore)  ! source

c
c  For particle concentration target elements, only count source elements
c  that are also particle concentrations.
c
            if( is_grp_mixed_comp(igroup)   .and. 
     $          (( itype(iefrom) .eq. I_COREMASS  ) .or.
     $             ( itype(iefrom) .eq. I_VOLCORE  ))  ) then

              !For mixed clouds, don't want to double count

            else
             nnucelem(ip) = nnucelem(ip) + 1
             inucelem(nnucelem(ip),ip) = ienconc( igelem(iefrom) )
            endif      ! ienconc element

c    No core 2nd moments for GROWCOREs
            if( im .ne. 0 .and. itype(iecore) .ne. I_GROWCORE )then

c    Check to see if source of second moment element should be the
c    second moment element in the source group 
              if( if_sec_mom(igelem(iefrom)) ) then
               if( if_nuc(imomelem(igelem(iefrom)),im) )
     $           iefrom = imomelem(igelem(iefrom))
              endif

              nnucelem(im) = nnucelem(im) + 1
              inucelem(nnucelem(im),im) = iefrom
            endif

          enddo
        enddo  	    ! jcore=1,ncore
      enddo         ! igroup=1,NGROUP

       write(*,*) 'Nucleation Element Mapping:'
       do iefrom = 1,NELEM
         write(*,*) iefrom,' nnucelem:',nnucelem(iefrom)
         do i=1,nnucelem(iefrom)
           write(*,*) '     ',i,' inucelem(i,ielem):',inucelem(i,iefrom)
         enddo
       enddo
       !stop
c
c  Now enumerate and count nucleating bins.
c
      do igroup = 1,NGROUP      ! target group
       do ibin = 1,NBIN         ! target bin

        do igfrom = 1,NGROUP    ! source group

        nnucbin(igfrom,ibin,igroup) = 0

         do ifrom = 1,NBIN      ! source bin

          if( inuc2bin(ifrom,igfrom,igroup) .eq. ibin ) then
           nnucbin(igfrom,ibin,igroup) = nnucbin(igfrom,ibin,igroup) + 1
           inucbin(nnucbin(igfrom,ibin,igroup),igfrom,ibin,igroup)
     $            = ifrom
          endif
         enddo
        enddo           ! igfrom=1,NGROUP
       enddo            ! ibin=1,NBIN=1,NGROUP
      enddo             ! igroup=1,NGROUP
c
c  Set up <ipownuc> array for use in <up/downgxfer>
c   <ipow> is the power to which the source particle mass must
c   be taken to match the type of the target element.
c
      do iefrom = 1,NELEM
        do ieto = 1,NELEM

         if( itype(iefrom) .eq. I_INVOLATILE .or.
     $       itype(iefrom) .eq. I_VOLATILE )then
           ipow_from = 0
         elseif ( itype(iefrom) .eq. I_COREMASS .or.
     $       itype(iefrom) .eq. I_GROWCORE .or.
     $       itype(iefrom) .eq. I_VOLCORE )then
           ipow_from = 1
         else
           ipow_from = 2
         endif

         if( itype(ieto) .eq. I_INVOLATILE .or.
     $       itype(ieto) .eq. I_VOLATILE )then
           ipow_to = 0
         elseif ( itype(ieto) .eq. I_COREMASS .or.
     $       itype(ieto) .eq. I_GROWCORE .or.
     $       itype(ieto) .eq. I_VOLCORE )then
           ipow_to = 1
         else
           ipow_to = 2
         endif

         ipownuc(iefrom,ieto) = ipow_to - ipow_from

        enddo
      enddo
c
c
c  Report nucleation mapping arrays (should be 'write' stmts, of course)
c
c      do iefrom=1,NELEM
c       if( nnuc2elem(iefrom) .ne. 0 ) 
c    $    write(*,*) 'iefrom, nnuc2elem(iefrom)=',
c    $                iefrom,nnuc2elem(iefrom)
c       do ieto = 1,NELEM
c        if( if_nuc(iefrom,ieto) )
c    $    write(*,*) 'if_nuc(iefrom,ieto) =  .true.',
c    $                iefrom,ieto
c       enddo
c      enddo
c
c      print*,' '
c      print*,'Nucleation mapping arrays (setupnuc):'
c      print*,' '
c      print*,'Elements mapping:'
c      do ielem = 1,NELEM
c       print*,' '
c       print*,'ielem,nnucelem=',ielem,nnucelem(ielem)
c       if(nnucelem(ielem) .gt. 0) then
c        do jfrom = 1,nnucelem(ielem)
c          print*,'jfrom,inucelem = ',jfrom,inucelem(jfrom,ielem)
c        enddo
c       endif
c      enddo
c      print*,' '
c      print*,'Bin mapping:'
c      do igfrom = 1,NGROUP
c       do igroup = 1,NGROUP
c        print*,' '
c        print *,'Groups (from, to) = ', igfrom, igroup
c        do ibin = 1,NBIN
c         nnucb = nnucbin(igfrom,ibin,igroup)
c         if(nnucb .eq. 0) print*,'  None for bin ',ibin
c         if(nnucb .gt. 0) then
c          print*,'  ibin,nnucbin=',ibin,nnucb
c        print*,'   inucbin=',(inucbin(j,igfrom,ibin,igroup),j=1,nnucb)
c         endif
c        enddo
c       enddo
c      enddo
c
c
c-----Check that values are valid------------------------------------------
c
c
      do ielem = 1, NELEM
 
        if( isolelem(ielem) .gt. NSOLUTE )then
          write(LUNOPRT,'(/,a)') 'component of isolelem > NSOLUTE'
          call endcarma
        endif
 
        if( ievp2elem(ielem) .gt. NELEM )then
          write(LUNOPRT,'(/,a)') 'component of ievp2elem > NELEM'
          stop 1
        endif
c
c
c  Check that <isolelem> is consistent with <ievp2elem>.
c 
        if( ievp2elem(ielem) .ne. 0 .and.
     $      itype(ielem) .eq. I_COREMASS )then
         if( isolelem(ielem) .ne. isolelem( ievp2elem(ielem) ) )then
          write(LUNOPRT,'(/,a)') 
     $      'isolelem and ievp2elem are inconsistent'
          call endcarma
         endif
        endif
c
c
c  Check that <isolelem> is consistent with <inucgas>.
c
        igas = inucgas( igelem(ielem),1 )
        if( igas .ne. 0 )then
          if( itype(ielem) .eq. I_COREMASS .and.
     $        isolelem(ielem) .eq. 0 )then
            write(LUNOPRT,'(/,a)') 'inucgas ne 0 but isolelem eq 0'
c           call endcarma
          endif
        endif

      enddo
 
      do ielem = 1, NELEM
        if( nnuc2elem(ielem) .gt. 0 ) then
          do inuc2 = 1, nnuc2elem(ielem)
            if( inuc2elem(inuc2,ielem) .gt. NELEM )then
              write(LUNOPRT,'(/,a)') 'component of inuc2elem > NELEM'
              call endcarma
            endif
          enddo
        endif
      enddo
c
c
c  Particle grids are incompatible if there is no target bin with enough
c  mass to accomodate nucleated particle.
c
      bad_grid = .false.

      do iefrom = 1,NELEM       ! source element
        igfrom = igelem(iefrom)
        neto   = nnuc2elem(iefrom)
        if( neto .gt. 0 )then
          do inuc2 = 1,neto
            ieto = inuc2elem(inuc2,iefrom)
            igto = igelem(ieto) 
            do ifrom = 1,NBIN   ! source bin
              if( inuc2bin(ifrom,igfrom,igto) .eq. 0 )then
                write(LUNOPRT,7) igfrom,ifrom,igto
                bad_grid = .true.
              endif
            enddo
          enddo
        endif
      enddo

      if( bad_grid )then
        write(LUNOPRT,'(/,a)') 'incompatible grids for nucleation'
        call endcarma
      endif
c
c
c--------------------------------------------------------------------------
c
c 
c  Define critical supersaturation and target bin for each (dry) particle
c  size bin that is subject to nucleation.
c  (only for CN groups subject to nucleation)
c
      do igroup = 1,NGROUP

        igas = inucgas(igroup,1)

        if( igas .ne. 0 .and.
     $      itype( ienconc( igroup ) ) .eq. I_INVOLATILE )then

          isol = isolelem( ienconc( igroup ) )
 
          do ibin = 1,NBIN
c
c
c  This is term "B" in Pruppacher and Klett's eqn. 6-28.
c
c           bsol = 3.*sol_ions(isol)*rmass(ibin,igroup)*gwtmol(igas)
c    $           / ( 4.*PI*solwtmol(isol)*RHO_W )
            bsol = 0.
c
c
c  Loop over vertical grid layers because of temperature dependence
c  in solute term.
c
            do k = 1,NZ

              scrit(k,ibin,igroup) = 0.0
c             scrit(k,ibin,igroup) = sqrt( 4. * akelvin(k,igas)**3
c    $                                   / ( 27. * bsol ) )
            enddo
          enddo
        endif
      enddo
c     
c
c--------------------------------------------------------------------------
c
c  Define contact parameter (source group, nucleating gas)
c
      do igroup = 1,NGROUP
        do igas = 1,NGAS
          ct(igas,igroup) = 0.
        enddo
      enddo
        ct(1,1)=0.986     !scrit = 1.15 for ethane onto tholin

        if( NGAS .gt. 1 ) then
cc       ct(2,2)=0.979     !scrit = 1.094 for methane onto ethane (T=45.9)
         ct(2,2)=0.981     !scrit = 1.094 for methane onto ethane (T=45.9)
         ct(2,1)=0.981     !scrit = 1.1 for methane onto tholin (T=45.0)
cc       ct(2,1)=0.8336    !scrit = 1.5 for methane onto tholin (T=72.7)

cc       ct(2,1)=0.8712     !scrit = 1.4  for methane
cc       write(*,*) 'Setting scrit methane to 1.4 !!!!!!!'
        endif

        if( gasname(1) .eq. 'methane' ) ct(1,1) = 0.981 !!!0.8336

c       ct(2)=0.967      !scrit = 1.094 (T=72.7) CH4/C2H6
c       ct(1)=0.9475     !scrit = 1.15 for methane
c       ct(1)=0.7265     !scrit = 2.0  for methane
c       ct(1)=0.9736     !scrit = 1.3

c       ct(1,1)=0.966      !scrit = 1.4
c       write(*,*) 'Setting scrit ethane to 1.4 !!!!!!!'

c       ct(1)=0.9589     !scrit = 1.5
c
  
        if( NCCNBIN .gt. 35 ) write(*,*) 'NCCNBIN > 35 !!!!!',NCCNBIN
c--------------------------------------------------------------------------
c
c  Define parameters needed for freezing nucleation calculations.
c
      adelf = 1.29e-12
      bdelf = 0.05
      prenuc = 2.075e33 * RHO_W / RHO_I
      rmiv   = 0.6
c
c
c--------------------------------------------------------------------------
c
c  Report some initialization values
c
      write(LUNOPRT,5)
c     write(LUNOPRT,1) 
c    $ 'inucgas  ',(((inucgas(j,i),i=1,NGROUP)),j=1,NGAS)
c     write(LUNOPRT,1) 'inuc2elem',(inuc2elem(1,i),i=1,NELEM)
      write(LUNOPRT,*) 'inucgas (# x group)'
      do j=1,NGAS
        write(LUNOPRT,1) '    ',(inucgas(j,i),i=1,NGROUP)
      enddo
      write(LUNOPRT,*) 'element   inuc2elem'
      do i=1,NELEM
       write(LUNOPRT,10) i,(inuc2elem(j,i),j=1,nnuc2elem(i)) 
      enddo
      write(LUNOPRT,1) 'ievp2elem',(ievp2elem(i),i=1,NELEM)
      write(LUNOPRT,1) 'isolelem ',(isolelem(i),i=1,NELEM)

      do isol = 1,NSOLUTE

        write(LUNOPRT,2) 'solute number   ',isol
        write(LUNOPRT,3) 'solute name:    ',solname(isol)
        write(LUNOPRT,4) 'molecular weight',solwtmol(isol)
        write(LUNOPRT,4) 'mass density    ',rhosol(isol)
        write(LUNOPRT,9) 'contact param   ',(ct(i,1),i=1,NGAS)

        do igroup = 1,NGROUP

          if( isol .eq. isolelem(ienconc(igroup)) )then

            write(LUNOPRT,8) groupname(igroup)

            write(LUNOPRT,6) (i,r(i,igroup),scrit(1,i,igroup),i=1,NBIN)

          endif
        enddo
      enddo
c
c
c  Return to caller with nucleation mapping arrays and critical
c  supersaturations defined.
c
      return
      end
