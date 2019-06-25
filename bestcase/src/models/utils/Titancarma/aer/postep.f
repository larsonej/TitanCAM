      subroutine postep
c
c
c  @(#) postep.f  McKie  Oct-1995
c  This routine handles all post-timestep processing  at the end
c  of every timestep.  Things that would appropriately be done
c  here include:
c    Smoothing or other corrections to the latest state.
c    Possible output of print file information on current step.
c    Possible output of history file information on current step.
c    Possible output of restart file information on current step.
c    Debugging output of changes within the step.
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
c  Local declarations.
c
      logical nuc_update
      logical last_step
      logical do_rad_now
c
c
c  Define formats
c
    1 format(i5,2x,3(1pe11.4))
    2 format(a,2(i5),4(2x,1pe13.6))
c
c 
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter postep'
c
c
c  Turn yearly printing back on after an extrapolation
c
      prt_year = .true.
c
c
c  Rescale vertical winds and diffusion coefficients for
c  other than cartesian coordinate systems
c
      if( igridv .ne. I_CART) call rescale
c
c
c  Update timestep index counter to end of current timestep
c
      itime = itime + 1
c
c
c  Update simulation time to end of current timestep
c
      time = time + dtime

c     if( itime .gt. 13656 ) nhist=1

c     write(*,*) '------- ',time/60.**2/24.,' days --------',itime
      
c
c  Don't allow particle concentrations to get too small.
c 
      do ixyz = 1,NXYZ
        do ibin = 1,NBIN
          do ielem = 1,NELEM
            call smallconc(ibin,ielem)
          enddo
        enddo
      enddo
c
c
c  Store amount of rainfall from this timestep and set up evaporation
c  flux for next timestep
c
      do igas=1,NGAS

        do n=1,nvolelem(igas)

          ielem = ivolelem(n,igas)

          ig = igelem(ielem)
            ! Note sign of pcmflux is negative for particles falling
            ! out of model, so need to subtract clouds for puddle to
            ! be positive
          if( ielem .eq. ienconc(ig) ) then
           !Add mass of cloud and subtract all cores
            puddle(igas) = puddle(igas) - pcmflux2(1,1,ielem)*dtime

            do ic=1,ncore(ig)
              iecore = icorelem(ic,ig)
              puddle(igas) = puddle(igas) + pcmflux2(1,1,iecore)*dtime
            enddo
          else
           !Volatile element is a growcore
            puddle(igas) = puddle(igas) - pcmflux2(1,1,ielem)*dtime
          endif

         !Correct for pc ~ SMALL_PC (i.e. where coremass may be
         ! greater than volatile mass)
          if( puddle(igas) .lt. ZERO ) puddle(igas) = ZERO

        enddo
        Rterm = RGAS/gwtmol(igas) * t3(1)
        pp = gc3(1,igas)*Rterm
       !Only calculate flux if subsaturated in bottom layer
        if( pp .lt. pvapl3(1,igas) )
     $   fbotevap(igas) = 1./(4.*Rterm) * sqrt( 8.*Rterm/PI ) *
     $                    ( pvapl3(1,igas) - pp )

c     Make sure enough liquid at surface for evaporation
        if( fbotevap(igas)*dtime .gt. puddle(igas) )
     $    fbotevap(igas) = puddle(igas)/dtime

cc        write(*,2) 'postep',itime,igas,fbotevap(igas),puddle(igas),
cc   $                 pcmflux2(1,1,2),pcmflux2(1,1,9)
      enddo

c
c
c  Compute a small fraction of the current timestep
c
      dt_small = .001 * dtime
c
c
c  Define local flag to indicate if this was last timestep for this run
c
      last_step = (itime .eq. ietime) .or. (time.ge.(endtime-.5*dtime))
c
c
c  Calculate radiative transfer for one column at at time
c  (frequency based on nrad when > 0, otherwise based on prad)
c
      if( do_rad )then

        do_rad_now = .false.

        if( last_step )then

          do_rad_now = .true.

        elseif( nrad .gt. 0 )then

           if( mod( itime, nrad ) .eq. 0 )then
             do_rad_now = .true.
           endif

        else

          if( mod( time + dt_small, prad ) .lt. dtime )then
            do_rad_now = .true.
          endif

        endif

        if( do_rad_now )then
          do iy = 1,NY
            do ix = 1,NX

              ixy = NX * ( iy - 1 ) + ix

              call prerad
              call radtran
              call postrad

            enddo
          enddo
        endif

      endif
c
c
c  Management of index of smallest bin nucleated from a group <inucmin>
c  and nucleation update time <time_nuc>:
c
c   Whenever <inucstep> is less than <inucmin>, update <inucmin> and <time_nuc>.
c   Otherwise, update <inucmin> and <time_nuc> every <period_nuc> seconds.
c
c   <time_nuc> is the previous time that <inucmin> was updated.  
c   If <inucmin> is updated this time step, then <time_nuc>
c   becomes the next time that it will be updated.
c
      nuc_update = .false.

      do ig = 1,NGROUP

        if( ( inucstep(ig) .lt. inucmin(ig) ) .or.
     $      ( time .ge. time_nuc ) )then

          inucmin(ig) = inucstep(ig)
          nuc_update  = .true.

        endif
      enddo

      if( nuc_update )then
        time_nuc = time + period_nuc
      endif
c
c
c  Possibly output print info from current timestep
c  (based on nprint when > 0, otherwise based on pprint)
c
      if( do_print )then

        if( last_step )then
 
          call outprt
 
        else if( nprint .gt. 0 )then
  
          if( nprint .le. ietime )then
           if( mod( itime, nprint ) .eq. 0 )then
             call outprt
           endif
          endif
   
        else
  
          if( mod(time+dt_small,pprint) .lt. dtime )then
           call outprt
          endif
   
        endif

      endif
c
c
c  Possibly output model state at current timestep
c  (based on nhist when > 0, otherwise based on phist)
c
c      if( do_hist )then
c
c        if( last_step )then
c 
c          call outhis
c 
c        else if( nhist .gt. 0 )then
c  
c          if( nhist .le. ietime )then
c           if( mod( itime, nhist ) .eq. 0 )then
c             call outhis
c           endif
c          endif
c   
c        else
c  
c          if( mod(time+dt_small,phist) .lt. dtime )then
c           call outhis
c
c    Possibly adjust <phist> (comment out when phist = const)
c           phist = 8.64d4
c
c          endif
c   
c        endif
c
c      endif
c
c
c  Possibly output model restart info at current timestep
c  (based on nrest when > 0, otherwise based on prest)
c
c      if( do_rest )then
c
c        if( last_step )then
c 
c          call outres
c 
c        else if( nrest .gt. 0 )then
c  
c          if( nrest .le. ietime )then
c           if( mod( itime, nrest ) .eq. 0 )then
c             call outres
c           endif
c          endif
c   
c        else
c  
c          if( mod(time+dt_small,prest) .lt. dtime )then
c           call outres
c          endif
c   
c        endif
c
c      endif
c
c
c  Possibly extrapolate to next time increment
c
c      if( mod(time+dt_small,tinc) .lt. dtime ) then
c        call extrapolate
c        prt_year = .false.
c      endif
c
c  Return to caller with end-of-timestep things completed
c
      return
      end
