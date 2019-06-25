      subroutine setuprest
c
c EJL 4-10-13 
c Putting the bottom half of setupaer here so that I can break up 
c the callse to the setups from carma.F90. I needed to do this in
c order to initialize that atmosphere before calling these routines.
c
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer'
c
c-------------------------------------------------------------------------------
c
c
c  Evaluate fall velocities.
c
      if (do_print_setup) print*, 'call setupvf'
      call setupvf
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for condensational growth and evaporation.
c
      if (do_print_setup) print*, 'call setupgrow'
      call setupgrow
      if (do_print_setup) print*, 'call setupgkern'
      call setupgkern
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for nucleation.
c
      if (do_print_setup) print*,'call setupnuc'
      call setupnuc

cc Print nucleation mapping arrays to screen for debugging, then stop
c     write(*,*) 'Liquid clouds - nucleation mapping'
c     do ielem = 1,NELEM
c       write(*,*) 'ielem,nnucelem',ielem,nnucelem(ielem)
c       do jfrom = 1,nnucelem(ielem)
c        write(*,*) 'jfrom,inucelem =',jfrom,inucelem(jfrom,ielem)
c       enddo
c     enddo
cc    stop
c
c  Set up array of element indexes for volatile elements corresponding
c  to each gas: <ivolelem>
      do igas = 1,NGAS 
        nvolelem(igas) = 0
        do n=1,NGROUP-1
          ivolelem(n,igas) = 0
        enddo
      enddo

      do igas = 1,NGAS
        n=1
        do ielem = 1,NELEM
          if( igrowgas(ielem) .eq. igas ) then
            ivolelem(n,igas) = ielem
            nvolelem(igas) = n
            n = n + 1

           !Don't include methane ice crystals if also droplets in this run
            if( icomp(ielem).eq.I_CH4_ICE .and. T0(igas) .lt. 81.)
     $         n = n-1
           !Don't include ethane growcore in mixed ice cloud if also mixed
           ! cloud droplets in this run
            if( itype(ielem).eq.I_GROWCORE .and.  
     $            is_grp_ice(igelem(ielem)) .and. T0(igas)
     $            .lt.81.) n=n-1
          endif
        enddo !elements
      enddo !igas
       if (do_print_setup) then 
        write(*,*) '# vol elements:',nvolelem(1),nvolelem(2)
       endif
c
c  Evaluate derived coagulation mapping arrays and kernels.
c
      if( do_coag ) then
      if (do_print_setup) print*,'call setupckern, setupcoag'
        call setupckern
        call setupcoag
      endif
c
c  Evaluate aerosol mass production
c
c      call setupmprod
c
c  Check particle setup for incompatibilities with other subroutines
c
      do igroup=1,NGROUP
        ncm = 0
        do ic=1,ncore(igroup)
          ie = icorelem(ic,igroup)
          if( itype(ie) .eq. I_COREMASS ) ncm = ncm + 1
        enddo
        if( ncm .gt. 1 ) then
          if( if_sec_mom(igroup) .and. is_grp_mixed_phase(igroup)) then
            write(*,*) '<setupaer> Particle group ',igroup,
     $         'incompatible with <evap_mono>.  Can only have one
     $          coremass element in a group with core 2nd moments
     $          and growcores'
            stop
          endif
        endif !more than one coremass element in group
      enddo
c
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined.
c
      return
      end
