      subroutine newstate
c
c
c  @(#) newstate.f  McKie  Oct-1995
c  This routine manages the calculations that update state variables
c  of the model with new values at the current simulation time.
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
c
c  Define formats
c
    1 format(/,'Timestep      time     max_ntsub   avg_ntsub'/)
    2 format(i6,3x,1pe9.2,3x,i6,3x,1pe13.2)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter newstate'
c
c
c  Calculate changes due to horizontal transport
c
c      if( do_ew .or. do_ns )then
c         call horizont
c      endif
c
c
c  Calculate changes in particle concentrations due to microphysical
c  processes, part 1.  (potentially slower microphysical calcs)
c  All spatial points are handled by one call to this routine.
c
      call microslow   !EJL - put this before vertical. It used to
c                       be before advection. 6-10-10
c	  
c  Calculate changes due to vertical transport
c
      if( do_vtran )then
        call vertical
      endif
c
c
c  Calculate transport forcings for parcel model
c
c      if( do_parcel )then
c        call parcel
c      endif
c
c
c  Calculate new saturation ratios for use in nsubsteps()
c
      do ixyz = 1,NXYZ
        if( do_thermo )then
          pt = ptc3(ixyz) / rhoa3(ixyz)
          t3(ixyz) = pt * ( p3(ixyz) / PREF )**RKAPPA
        endif
        call supersat
      enddo
c
c
c  Calculate the changes in concentrations and supersaturation
c  due to transport
c
      do ixyz = 1,NXYZ

        d_ptc(ixyz) = ptc3(ixyz) - ptcl(ixyz)

        do igas = 1,NGAS
          d_gc(ixyz,igas) = gc3(ixyz,igas) - gcl(ixyz,igas)
        enddo
      enddo
c
c
c  Reset concentrations to their (pre-advection) values from the
c  beginning of the time-step
c
      do ixyz = 1,NXYZ

        ptc3(ixyz) = ptcl(ixyz)
        if( do_thermo )then
          pt = ptc3(ixyz) / rhoa3(ixyz)
          t3(ixyz) = pt * ( p3(ixyz) / PREF )**RKAPPA
        endif

        do igas = 1,NGAS
          gc3(ixyz,igas) = gcl(ixyz,igas)
        enddo
      enddo
c
c
!      call microslow !EJL 6-11-10
c
c
c  Save current timestep used for processes so far
c
      dtime_save = dtime
c
c
c  Initialize diagnostics for substepping
c
      max_ntsub = 0
      avg_ntsub = 0.
c
c  Set vertical loop index to increment downwards
c  (for substepping of sedimentation)
c
      if( igridv .eq. I_CART )then
        kb  = NZ
        ke  = 1
        idk = -1
      else
        kb  = 1
        ke  = NZ
        idk = 1
      endif
c
c
c  Loop over horizontal grid poins
c
      do ix = 1,NX
        do iy = 1,NY
c
c
c  Loop over vertical grid poins
c
          do iz = kb,ke,idk
            ixy = NX * ( iy - 1 ) + ix 
            ixyz = NXY * ( iz - 1 ) + ixy
c
c
c  Compute or specify number of sub-timestep intervals for current spatial point
c  (Could be same for all spatial pts, or could vary as a function of location)
c
            call nsubsteps
            max_ntsub = max( max_ntsub, ntsubsteps )
            avg_ntsub = avg_ntsub + float(ntsubsteps)/float(NXYZ)
c
c
c  Compute sub-timestep time interval for current spatial grid point
c
            dtime = dtime_save / ntsubsteps
c
c
c  Do sub-timestepping for current spatial grid point
c
            do isubstep = 1,ntsubsteps
c
c
c  Apply the fractional advective forcing to specify the values of
c  concentrations at this sub-timestep
c
              ptc3(ixyz) = ptc3(ixyz) + d_ptc(ixyz)/ntsubsteps

c  NOTE: The addition of this small number can create small mass
c  conservation errors, since gc3 will not have enough precision
c  to reatin all the digits of d_gc/ntsubsteps. To conserve mass
c  to machine precision, the remainder of each of these additions
c  should be retained and added back to gc3 after the ntsubstebs
c  have been processed.
              do igas = 1,NGAS
                gc3(ixyz,igas) = gc3(ixyz,igas) +
     $                           d_gc(ixyz,igas)/ntsubsteps
              enddo
c
c
c  Update temperature
c
              if( do_thermo )then
                pt = ptc3(ixyz) / rhoa3(ixyz)
                t3(ixyz) = pt * ( p3(ixyz) / PREF )**RKAPPA
              endif
c
c
c  Calculate changes in particle concentrations for current spatial point
c  due to microphysical processes, part 2.  (faster microphysical calcs)
c
              call microfast
 
            enddo
          enddo
        enddo
      enddo
c
c
c  Report substep diagnostics
c
c      if( itime .eq. 0 ) write(LUNOSTEP,1)
c      write(LUNOSTEP,2) itime, time, max_ntsub, avg_ntsub
c
c
c  Restore normal timestep
c
      dtime = dtime_save
c
c
c  Hydrostatically update atmospheric profiles
c
c      if( .not. do_parcel )then
c        call hydrostat
c      endif
c
c
c  Return to caller with new state computed 
c
      return
      end