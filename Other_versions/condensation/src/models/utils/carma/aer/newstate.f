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
c  Local declarations.
c
      logical last_step
c
c
c  Define formats
c
    1 format(/,'Timestep     time     max_ntsub   avg_ntsub'/)
    2 format(i10,3x,1pe8.2,3x,i6,3x,1pe13.2)
    5 format(a,i4,2x,'gc/vol:',f15.10,2x,'gc/cld:',f15.10,2x,
     $       'cmf:',f15.10,2x,'cloud(pc):',1pe13.6,i4)
    6 format(a,2(i4),2x,'cld:',1pe13.6,2x,'cor:',1pe13.6,2x,
     $       'gcr:',1pe13.6)
    7 format(a,1pe11.4,2x,i4,' bin ',i3)
    8 format(a,1pe11.4,' dpc: ',1pe11.4,' nsub: ',i3,' bin ',i3)
    9 format(a,2(i4),'  init: ',f15.10,'  vtran: ',f15.10,
     $       '  coag: ',f15.10,'  grow: ',f15.10,'  evap: ',
     $       f15.10,'  exc: ',f15.10)
   10 format(a,2(i4),'  init: ',1pe11.4,'  vtran: ',1pe11.4,
     $       '  coag: ',e11.4,'  grow: ',e11.4,'  evap: ',
     $       e11.4,'  exc: ',e11.4)
   11 format(a,2(i4),'  init: ',1pe11.4,'  vtran: ',1pe11.4,
     $       '  coag: ',e11.4,'  grow: ',e11.4,'  evap: ',
     $       e11.4,0p,i8,f10.5,' days')
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter newstate'

c----------------------------------------------------------
c
c  Calculate changes due to horizontal transport
c
c       if( do_ew .or. do_ns )then
c          call horizont
c       endif
c
cc      if( itime.gt.10000 ) do_vtran = .false.

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
c  Calculate new saturation ratios for use in nsubsteps.f
c
      do ixyz = 1,NXYZ
c       pt = ptc3(ixyz) / rhoa3(ixyz)
c       t3(ixyz) = pt * ( p3(ixyz) / PREF )**RKAPPA
        call supersat
      enddo
c
c
c  Calculate the changes in concentrations and supersaturation
c  due to transport
c
      if( do_ew .or. do_ns .or. do_vtran .or. do_parcel )then


        do ixyz = 1,NXYZ

          d_ptc(ixyz) = ptc3(ixyz) - ptcl(ixyz)

            do igas = 1,NGAS
            d_gc(ixyz,igas) = gc3(ixyz,igas) - gcl(ixyz,igas)
          enddo

c          do ielem = 1,NELEM
c            do ibin = 1,NBIN
c              d_pc(ixyz,ibin,ielem) = pc3(ixyz,ibin,ielem) -
c     $                                pcl(ixyz,ibin,ielem)
c
c            enddo
c          enddo

        enddo
c
c
c  Reset concentrations to their (pre-advection) values from the
c  beginning of the time-step
c
        do ixyz = 1,NXYZ

          ptc3(ixyz) = ptcl(ixyz)

          do igas = 1,NGAS
            gc3(ixyz,igas) = gcl(ixyz,igas)
          enddo

          do ielem = 1,NELEM
            do ibin = 1,NBIN
              pc3(ixyz,ibin,ielem) = pcl(ixyz,ibin,ielem)
            enddo
          enddo

        enddo

      else

        do ixyz = 1,NXYZ

          d_ptc(ixyz) = 0.

            do igas = 1,NGAS
            d_gc(ixyz,igas) = 0.
          enddo

          do ielem = 1,NELEM
            do ibin = 1,NBIN
              d_pc(ixyz,ibin,ielem) = 0.
            enddo
          enddo

        enddo

      endif
c
c
c  Calculate changes in particle concentrations due to microphysical
c  processes, part 1.  (potentially slower microphysical calcs)
c  All spatial points are handled by one call to this routine.
c
      call microslow

c----------------------------------------------------------
c  Save growcore mass fraction for debugging

cc    k = iwa
cc    do i=1,NBIN
cc      volpart = pc3(k,i,2)*rmass(i,2) -  pc3(k,i,3)
cc      rfrac_sav(i,3) = pc3(k,i,5) / volpart
cc      rfrac_sav(i,3) = pc3(k,i,2)*rmass(i,2)
cc   $                    -pc3(k,i,3)+pc3(k,i,8)
cc    enddo
c   Save tholin bin 1 for debugging
c     rfrac_sav(1,3) = pc3(21,1,1)
c----------------------------------------------------------
c
c
c  Save current timestep used for processes so far
c
      dtime_save = dtime
c
c
c  Initialize diagnostics for substeping
c
      max_ntsub = 0
      avg_ntsub = 0.
c
c
c  Visit each of the spatial grid points
c
      do iz = 1,NZ
        do iy = 1,NY
          do ix = 1,NX
c
c
c  Compute linearized spatial grid pt indices for horiz 2-D & overall 3-D
c   (Global vars <ix>, <iy>, <iz>, <ixy>, <ixyz> used by various microphysical routines)
c
            ixy = NX * ( iy - 1 ) + ix 
            ixyz = NXY * ( iz - 1 ) + ixy
c
c
c  Compute or specify number of sub-timestep intervals for current spatial point
c  (Could be same for all spatial pts, or could vary as a function of location)
c
            ntsubsteps = minsubsteps
            call nsubsteps
            max_ntsub = max( max_ntsub, ntsubsteps )
            avg_ntsub = avg_ntsub + float(ntsubsteps)/float(NXYZ)
c
c
c  Compute sub-timestep time interval for current spatial grid point
c
            dtime = dtime_save / float( ntsubsteps )
c
c
c  Do sub-timestepping for current spatial grid point
c
            do isubtim = 1,ntsubsteps
c
c
c  Apply the fractional advective forcing to specify the values of
c  concentrations at this sub-timestep
c
              ptc3(ixyz) = ptc3(ixyz) + d_ptc(ixyz)/ntsubsteps

              do igas = 1,NGAS
                gc3(ixyz,igas) = gc3(ixyz,igas) +
     $                           d_gc(ixyz,igas)/ntsubsteps

              enddo

              do ielem = 1,NELEM
                do ibin = 1,NBIN

                  pc3(ixyz,ibin,ielem) = pc3(ixyz,ibin,ielem) + 
     $                           d_pc(ixyz,ibin,ielem)/ntsubsteps
c                 
c    
c  Prevent particle concentrations from dropping below SMALL_PC
c  (In some cases, this prevention of zero concentrations
c  may be needed to prevent numerical problems)
c
                   call smallconc(ibin,ielem)

                enddo
              enddo
c
c
c  Update temperature
c
c             pt = ptc3(ixyz) / rhoa3(ixyz)
c             t3(ixyz) = pt * ( p3(ixyz) / PREF )**RKAPPA
c
c
c  Fix any bad growcore fractions from coagulation and transport
c
              call fixfrac
c
c
c  Calculate changes in particle concentrations for current spatial point
c  due to microphysical processes, part 2.  (faster microphysical calcs)
c
cc       if(itime.ge.14692.and. ixyz.eq.14)
cc   $    write(*,*) itime,' before microfast:',pc3(ixyz,48,7),
cc   $                 pc3(ixyz,48,8),pc3(ixyz,48,9)
              call microfast
cc       if(itime.ge.14692.and. ixyz.eq.14)
cc   $    write(*,*) itime,' after microfast:',pc3(ixyz,48,7),
cc   $                 pc3(ixyz,48,8),pc3(ixyz,48,9)
c
c
c  Fix any bad growcore fractions from growth and evaporation
c
              call fixfrac
c
c
c  Go do next sub-timestep
c
            enddo

cc       if(itime.gt.14690 .and. ixyz.eq.14)
cc   $    write(*,*) itime,' gc after fixfrac:',gc3(ixyz,1),
cc   $                 gc3(ixyz,2)
c
c
c  Define local flag to indicate if this was last timestep for this run
c
            last_step = (itime+1 .eq. ietime) .or.
     $                  (time+dtime.ge.(endtime-.5*dtime))
c
c  Small fraction of the current timestep
            dt_small = .001 * dtime    
c
c  Possibly output rates info from current timestep (uncomment one)
c
c    Print rates yearly
c
c           if( do_print )then

c             if( last_step )then

c               call prtrates

c             else
c               if( mod(time+dtime+dt_small,YEAR) .lt. dtime
c    $                                    .and. prt_year) then
c                 call prtrates
c               endif

c             endif

c           endif

c
c    Print rates with pprint or nprint frequency 
c
      if( do_print )then
 
c          call prtrates
 
        else if( nprint .gt. 0 )then
 
          if( nprint .le. ietime )then
           if( mod( itime, nprint ) .eq. 0 )then
c             call prtrates
           endif
          endif
 
        else
 
          if( mod(time+dt_small,pprint) .lt. dtime )then
c           call prtrates
          endif
 
        endif
c
c
c  Go do next spatial grid point
c
          enddo
        enddo
      enddo
c
c
c  Restore normal timestep
c
      dtime = dtime_save
c
c
c  Hydrostatically update atmospheric profiles
c
      if( .not. do_parcel )then
        call hydrostat
      endif

c
c
c  Return to caller with new state computed 
c
      return
      end
