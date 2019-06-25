      subroutine varstep
c
c
c  @(#) varstep.f  Jensen  Mar-1996
c  This routine adjust the time-step depending upon how much
c  the concentrations changed during the time-step.
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
c  Local variables
c
      logical negative_gc
c
c
c  Define formats
c
    1 format(/,'Timestep     dtime    time     dpcmax   ixyzmax',
     $         '  ielemax  ibinmax    dsmax   ixyzsmax'/)
    2 format(i10,3x,1pe8.2,3x,1pe8.2,3x,1pe9.2,3x,
     $       i3,5x,i4,5x,i4,5x,1pe8.2,4x,i4)
c
c---------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter varstep'
c
c
c  Set max changes to small values
c
      dpcmax = 0.
      dsmax = 0.
      courmax = 0.

      ielemax = 0
      ibinmax = 0
      ixyzpmax = 0
      ixyzsmax = 0
c
c
c  Find maximum changes in particle concentration.
c  Only consider bins where the concentration is greater than
c  <conmax> times the peak concentration <pcmax(ielem)>.
c  Also, only consider particle number density elements.
c
      do ielem = 1,NELEM
       if( ( itype(ielem) .eq. I_INVOLATILE ) .or.
     $     ( itype(ielem) .eq. I_VOLATILE ) )then
        do i = 1,NBIN
         do ixyz = 1,NXYZ

c         Possibly change dtmax if particles are falling too quickly
           if( pc3(ixyz,i,ielem) .gt. 1.d-25 .and. do_vtran ) then
            if(dz3(ixyz)/vf(ixyz,i,igelem(ielem)) .lt. dtmax ) 
     $       dtmax = 0.9d0 * dz3(ixyz)/vf(ixyz,i,igelem(ielem))
           endif

c          if( ( pc3(ixyz,i,ielem) .ge. conmax*pcmax(ielem) .or.
c     $        pcl(ixyz,i,ielem) .ge. conmax*pcmax(ielem) ) .and.
c     $        pcmax(ielem) .gt. 1.e-2 )then
          if( ( pc3(ixyz,i,ielem) .ge. conmax*pcmax(ielem) .and.
     $        pcl(ixyz,i,ielem) .ge. conmax*pcmax(ielem) ) .and.
     $        pcmax(ielem) .gt. 1.e-2 )then
c          if( ( ( pc3(ixyz,i,ielem) .ge. conmax*pcmax(ielem) .and.
c     $            pcl(ixyz,i,ielem) .ge. 1.e6*SMALL_PC ) .or.
c     $          ( pcl(ixyz,i,ielem) .ge. conmax*pcmax(ielem) .and.
c     $          ( pc3(ixyz,i,ielem) .ge. 1.e6*SMALL_PC ) ) ) .and.
c     $        pcmax(ielem) .gt. 2.e-2 )then

           pcrat = pcl(ixyz,i,ielem)/pc3(ixyz,i,ielem)

           if( pcrat .ge. 1. .and. i .gt. 1 )then
             crat = abs( pcrat - 1. )
           else
             crat = abs( 1./pcrat - 1. )
           endif

           if( crat .gt. dpcmax )then
            dpcmax = crat
            ielemax = ielem
            ibinmax = i
            ixyzpmax = ixyz
           endif

          endif
         enddo ! 1,NXYZ
        enddo  ! 1,NBIN
       endif   ! I_INVOLATILE, I_VOLATILE
      enddo
c
c
c  Find maximum changes in volatile core mass concentration.
c  Only consider bins where the concentration is greater than
c  <conmax> times the peak concentration <pcmax(ielem)>.
c  Also, only consider bins with core mass fraction > 1.e-4.
c
      do ielem = 1,NELEM
       if( itype(ielem) .eq. I_VOLCORE .or.
     $     itype(ielem) .eq. I_GROWCORE ) then

        igroup = igelem( ielem )
        iepart = ienconc( igroup )

        do i = 1,NBIN
         do ixyz = 1,NXYZ

          core_mass_frac = pc3(ixyz,i,ielem) /
     $            ( pc3(ixyz,i,iepart) *rmass(i,igroup) )

          if( ( pc3(ixyz,i,iepart) .ge. conmax*pcmax(iepart) .or.
     $        pcl(ixyz,i,iepart) .ge. conmax*pcmax(iepart) ) .and.
     $        pcmax(iepart) .gt. 1.e-2 .and.
     $        core_mass_frac .gt. 1.e-1 )then

           pcrat = pcl(ixyz,i,ielem)/pc3(ixyz,i,ielem)

           if( pcrat .ge. 1. .and. i .gt. 1 )then
             crat = abs( pcrat - 1. )
           else
             crat = abs( 1./pcrat - 1. )
           endif

           if( crat .gt. dpcmax )then
            dpcmax = crat
            ielemax = ielem
            ibinmax = i
            ixyzpmax = ixyz
           endif

          endif
         enddo
        enddo
       endif
      enddo
c
c
c  Find maximum changes in supersaturation and
c  check for negative gas concentrations
c  (note: need to specify ice or liquid).
c
      if( do_grow )then
        negative_gc = .false.
        do igas = 1,NGAS
          do ixyz = 1,NXYZ

cc          if( t3(ixyz) .ge. T0(igas) .and.        !! comment these out again
cc   $          gasname(igas) .eq. 'methane') then  !! if not looking at the
cc            supsatold = supsatlold(ixyz,igas)     !! CH4 + N2 case
cc            supsatnew = supsatl3(ixyz,igas)       !!
cc          else                                    !!
              supsatold = supsatiold(ixyz,igas)
              supsatnew = supsati3(ixyz,igas)
cc          endif                                   !!
            if( supsatold .ge. 1.d-4 )then
              srat1 = abs( supsatnew / supsatold - 
     $                1. )
            else
              srat1 = 0.
            endif

            if( supsatnew .ge. 1.d-4 )then
              srat2 = abs( supsatold / supsatnew -
     $                1. )
            else
              srat2 = 0.
            endif

            srat = max( srat1, srat2 )

            if( srat .gt. dsmax )then
              dsmax = srat
              ixyzsmax = ixyz
            endif

            if( gc3(ixyz,igas) .lt. 0.d0) then
              negative_gc = .true.
              ialt_neg = ixyz  !last alt bin where -gc's found
              igas_neg = igas
            endif

          enddo
        enddo
      endif

c     if( do_write(7) ) then
       if(negative_gc) then
c       write(*,*) 'Negative Gas Concentrations at ',
c    $       time/60.**2/24.,' days',ialt_neg,dtime
c      if(dtime .le. dtmin) write(*,*) 'reached dtmin'
       if( gcl(ialt_neg,igas_neg) .lt. 0.d0 ) 
     $   negative_gc = .false.
       endif

c     endif
c
c
c  Find maximum horizontal advection Courant number
c
      if( do_ew .or. do_ns ) then
       do iy = 1,NY
        do ix = 1,NX
         do k = 1,NZ

          courx = u(ix,iy,k) * dtime / dx(ix,iy,k)
          if( courx .gt. courmax ) courmax = courx
          coury = v(ix,iy,k) * dtime / dy(ix,iy,k)
          if( coury .gt. courmax ) courmax = coury

         enddo
        enddo
       enddo
      endif
c
c
c  Possibly adjust time-step 
c
      do_step = .true.
      if( ( dpcmax .lt. dpctol/5. .and. dsmax .lt. dgstol/4. 
     $    .and. (.not. negative_gc) ) 
     $    .and. dtime .lt. dtmax .and. courmax .lt. 0.5  )then 

        dtime = min( dtime * 2., dtmax )
c
c
c  Re-do current timestep with smaller timestep size
c
      else if( ( dpcmax .gt. dpctol .or. dsmax .gt. dgstol .or.
     $         courmax .gt. 0.5  .or. negative_gc ) .and. 
     $                                  dtime .gt. dtmin )then

        dtime = dtime / 2.03
        do_step = .false.
c
c
c  Restore old concentrations etc
c
c       write(*,*) 'Varstep restoring old concentrations'

        do ielem = 1,NELEM
          do i = 1,NBIN
            do ixyz = 1,NXYZ
              pc3(ixyz,i,ielem) = pcl(ixyz,i,ielem)
            enddo
          enddo
        enddo
        do igas = 1,NGAS
          do ixyz = 1,NXYZ
            gc3(ixyz,igas) = gcl(ixyz,igas)
            supsatl3(ixyz,igas) = supsatlold(ixyz,igas)
            supsati3(ixyz,igas) = supsatiold(ixyz,igas)
          enddo
        enddo
        do ixyz = 1,NXYZ
          ptc3(ixyz) = ptcl(ixyz)
          zmet3(ixyz) = zmetl3(ixyz)
          t3(ixyz) = told(ixyz)
          p3(ixyz) = pold(ixyz)
          rhoa3(ixyz) = rhoaold(ixyz)
        enddo
        if( do_parcel )then
          do ixyz = 1,NXYZP1
            zl3(ixyz) = zlold(ixyz)
          enddo
          do ixyz = 1,NXYZ
            zc3(ixyz) = zcold(ixyz)
          enddo
        endif

      endif
c
c
c  Check to see whether <endtime> will be exceeded in next time
c  step.  If so, then reduce <dtime> such that <endtime> is
c  reached exactly.
c
      if( time+dtime .gt. endtime )then
        dtime = endtime - time
      endif
c
c
c  Report timestep diagnostics
c
cc    if( itime .eq. 0 ) write(LUNOSTEP,1)
cc    write(LUNOSTEP,2) itime, dtime, time, dpcmax, ixyzpmax,
cc   $                  ielemax, ibinmax, dsmax, ixyzsmax
c
c
c  Return to caller with timestep modified.
c
      return
      end
