      subroutine prestep
c
c
c  @(#) prestep.f  McKie  Oct-1995
c  This routine handles all preliminary setup at the beginning
c  of every timestep.  Things that would appropriately be done
c  here include:
c    Input or otherwise define interface quantities from other submodels.
c    Save any model state that is needed to compute tendencies.
c    Save any model state that might be needed for comparison at end of step.
c    Update timestep counter and simulation time.
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
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter prestep'
c
c
c  Don't allow particle concentrations to get too small.
c  --- moved to postep ---  
c     do ixyz = 1,NXYZ
c       do ibin = 1,NBIN
c         do ielem = 1,NELEM
c           call smallconc(ibin,ielem)
c         enddo
c       enddo
c     enddo
c
c
c  Set <pcl> to <pc> from previous time step
c  and find maximum for each element.
c
      do ielem = 1,NELEM
        pcmax(ielem) = SMALL_PC
        do i = 1,NBIN
          do ixyz = 1,NXYZ
c            pcl(ixyz,i,ielem) = pc3(ixyz,i,ielem)
            pcmax(ielem) = max( pcmax(ielem), pc3(ixyz,i,ielem) )
          enddo
        enddo
      enddo
c
c
c  Find maximum particle concentration for each spatial grid box
c  (in units of cm^-3)
c
      do igroup = 1,NGROUP
        iep = ienconc(igroup)
        do ixyz = 1,NXYZ
          pconmax(ixyz,igroup) = SMALL_PC
          do i = 1,NBIN
            pconmax(ixyz,igroup) = max( pconmax(ixyz,igroup),
     $                                  pc3(ixyz,i,iep) )
          enddo
          pconmax(ixyz,igroup) = pconmax(ixyz,igroup) /
     $         ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )
        enddo
      enddo
c
c
c  Set <gcl> to <gc> from previous time step
c  and store old ice and liquid supersaturations for variable time stepping
c
c      do igas = 1,NGAS
c        do ixyz = 1,NXYZ
c          gcl(ixyz,igas) = gc3(ixyz,igas)
c          supsatlold(ixyz,igas) = supsatl3(ixyz,igas)
c          supsatiold(ixyz,igas) = supsati3(ixyz,igas)
c        enddo
c      enddo
c
c
c  Set <gc> so that bottom layer remains at 75% relative humidity (ethane)
c
c     gc3(1,1) = pvapi3(1,1) * 0.75 * gwtmol(1)/RGAS /t3(1)
c
c  Set <gc> (methane) so that it doesn't go below Lellouch's value
c
cc    gc3(1,2) = max( gc3(1,2) , 2.356d-4 )
c
c
c  Store old values of atmospheric properties for variable time stepping.
c
c      do ixyz = 1,NXYZ
c        zmetl3(ixyz) = zmet3(ixyz)
c      enddo
c
c      if( do_varstep )then
c
c        do ixyz = 1,NXYZ
c          ptcl(ixyz) = ptc3(ixyz)
c          told(ixyz) = t3(ixyz)
c          pold(ixyz) = p3(ixyz)
c          rhoaold(ixyz) = rhoa3(ixyz)
c        enddo
c
c        if( do_parcel )then
c          do ixyz = 1,NXYZP1
c            zlold(ixyz) = zl3(ixyz)
c          enddo
c          do ixyz = 1,NXYZ
c            zcold(ixyz) = zc3(ixyz)
c          enddo
c        endif
c
c      endif
c
c
c
c  Set production terms and loss rates due to slow microphysics
c  processes (coagulation) to zero.
c
      if( ibtime .eq. 0 .or. do_coag )then

        do i = 1,NBIN

          do ielem = 1,NELEM
            do ixyz = 1,NXYZ
              coagpe(ixyz,i,ielem) = 0.
            enddo
          enddo

          do igroup = 1,NGROUP
            do ixyz = 1,NXYZ
              coaglg(ixyz,i,igroup) = 0.
            enddo
          enddo

        enddo

      endif
c
c
c  Set production terms and loss rates due to fast microphysics
c  processes to zero if just starting a simulation.
c
      if( ibtime .eq. 0 )then
        call zeromicro
      endif

c  Adjust vertical wind speed profile (hack for test simulation) 

c     z_peak = 5.5e5 + 5.e2*time
c     do k=3,NZ
c       do ix=1,NX
c         do iy=1,NY
c           w(ix,iy,k) = 5.e2 *
c    $                   exp( -1.*(zc(ix,iy,k)-z_peak)**2 /
c    $                                 1.e5**2 )
c         enddo
c       enddo
c     enddo
c
c
c  Parcel Model Test -- 
c
c      if(do_parcel) then
c        w3(1) = 500.d0
c            if ( itime .eq. int(endtime/5./dtime) ) then
c              w3(ixyz) = -500.d0
c            endif
c      endif
c
c
c  Test Case -- Temperature perturbation for supersat when vertical
c               transport turned off
c
c     if( t3(21) .lt. 74.6 ) then
c       if( test_mass_cons .and. .not. do_vtran ) then
c        do ixyz = 4,12 !1,NXYZ
c       t3(ixyz) = tinit(ixyz) + 0.5d0*sin(1.26d-2 * time)   !p = 500 s 
c       t3(ixyz) = tinit(ixyz) + 1.0d0*sin(1.26d-2 * time)   !p = 500 s 
c       t3(ixyz) = tinit(ixyz) + 2.0d0*sin(6.28d-2 * time)   !p = 1000 s 
c       t3(ixyz) = tinit(ixyz) + 1.0d0*sin(7.27d-5 * time)   !1 deg / day
c       t3(ixyz) = tinit(ixyz) + 1.0d0*sin(1.04d-5 * time)   !1 deg / week
c       t3(ixyz) = tinit(ixyz) - 1.d0*sin(2.42d-6 * time)   !p = 1 month
c        t3(ixyz) = tinit(ixyz) - 2.d0*sin(1.21d-6 * time)   !p = 2 months
c       if( time/60.**2/24. .gt. 15. .and. 
c    $      t3(ixyz) .lt. told(ixyz)       )
c    $      t3(ixyz) = told(ixyz) + 0.5
c        enddo
c       endif
c
c  Return to caller with preliminary timestep things completed.
c
      return
      end
