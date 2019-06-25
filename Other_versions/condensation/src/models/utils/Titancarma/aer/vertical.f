      subroutine vertical
c
c
c  @(#) vertical.f  Jensen  Mar-1995
c  This routine drives the vertical transport calculations.
c
c  Modified  Sep-1997  (McKie)
c  Remove <ixy> from arg list of called routines.
c  <ixy> now available as a global var in common block.
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables.
c
      include 'globaer.h'
c
c
c  Declare local variables
c
      dimension drho_dt(NZ)

c
c  Define formats
c
    1 format(4e20.13)
    2 format(5e20.13)
    3 format(i3,3x,6e13.4)
    4 format(i3,3x,i3,3x,i3,3x,6f9.2)
    5 format(i4,i4,1pe13.6)
    6 format(i3,f7.1,i6,1pe13.4)
    7 format(3(1pe13.6))
    8 format(a,3(i4),2x,'mflux:',1pe13.6,2x,'advd:',
     $       1pe13.6,0p,2x,f4.2)
    9 format(a,2x,1pe13.6,2(i5))
   10 format(a,2x,1pe13.6,2(i5),i8,0p,f10.5)
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter vertical'
c 
c=======================================================================
c
c     Some values for mass flux output
c
      dt_small = .001 * dtime    !small fraction of the current timestep
c
c  Determine limits for vertical loops: For 1-D simulations, we do not
c  want the concentration to change due to flux out of a boundary
c  grid box (i.e., upward motion for the bottom boundary and
c  downward motion for the top boundary).  This constraint ensures
c  that advection in a 1-D simulation will not change concentrations
c  at the upwind boundary.
c
      nlowlim = 1
      nuplim = NZ
      if( NXY .eq. 1 ) then
        if( w2(1,2) .gt. 0. ) nlowlim = 2
        if( w2(1,NZ) .lt. 0. ) nuplim = NZ-1
      endif
c
c  Loop over horizontal grid points
c
      do ixy = 1,NXY
c
c=======================================================================
c
c  Treat vertical transport of air density 
c
c
c  Boundary fluxes are zero
c
c        itbnd = I_FLUX_SPEC
c        ibbnd = I_FLUX_SPEC
c
c        ftop = 0.
c        fbot = 0.
c
c
c  Store temporary (work) values in <cvert>, 
c  evaluate vertical velocities at layer boundaries,
c  and set divergence correction term to zero.
c
c        do k = 1,NZP1
c          vtrans(k) = w2(ixy,k)
c          if( k .le. NZ )then
c            cvert(k) = rhoa2(ixy,k)
c            divcor(k) = 0.
c          endif
c        enddo
c
c
c  Calculate particle transport rates due to vertical advection
c  (note: no diffusion for air density), and solve for air density
c  at end of time step.
c
c        call vertadv
c        call versol
c
c        if( NXY .eq. 1 )then
c
c
c  In 1D, assume assume any vertical divergence is
c  balanced by horizontal convergence: don't allow air density
c  to change in time, but calculate rate of change that would have
c  resulted from advection -- this tendency is then used below for
c  a divergence correction.
c
c          do k = 1,NZ
c            drho_dt(k) = ( rhoa2(ixy,k) - cvert(k) ) / 
c     $                   ( rhoa2(ixy,k) * dtime )
c          enddo
c
c        else
c
c
c  Update air density when not in 1D.
c
c          do k = nlowlim,nuplim
c            rhoa2(ixy,k) = cvert(k)
c          enddo
c
c        endif
c
c=======================================================================
c
c  Treat vertical transport of particles.
c
        itbnd = itbnd_pc
        ibbnd = ibbnd_pc


        do ielem = 1,NELEM          ! Loop over particle elements

          ig = igelem(ielem)        ! particle group

          do k = 1,NZ
            pcmflux2(ixy,k,ielem) = 0.
          enddo

          do ibin = 1,NBIN          ! Loop over particle mass bins
c
c
c  Store temporary (work) values in <cvert>, 
c  evaluate vertical velocities at layer boundaries, and
c  when 1D, assign divergence correction term
c
            do k = 1,NZP1
              vtrans(k) = w2(ixy,k) - vf(k,ibin,ig)
              if( k .le. NZ )then
                cvert(k) = pc2(ixy,k,ibin,ielem)
c                if( NXY. eq. 1 ) divcor(k) = cvert(k) * drho_dt(k)
c             endif
            enddo

c
c
c  Fluxes across top and bottom of model
c
c   Allow both tholin and cloud particles to deposit onto the surface
c   Flux of particles into top of model is only for tholin

            ftop = ftoppart(ixy,ibin,ielem)
            fbot = fbotpart(ixy,ibin,ielem)+vtrans(1)*cvert(1)

c           if(ielem .eq. 1) print *,'Fluxes',ielem,ibin,ftop, fbot
c           print *,'Flux',ielem,r(ibin,ig),fbot
c
c  Calculate particle transport rates due to vertical advection
c  and vertical diffusion, and solve for concentrations at end of time step.
c
            call vertadv
c            call vertdif
            call versol 

            if( test_mass_cons .and. itime.eq.ibtime+1
     $         .and. ibin.eq.1 ) write(*,*) ftop,fbot,ielem
c
c
c   Calculate particle mass flux
c
            if(ielem .eq. ienconc(ig)) then
c
c             Particle number concentration 
c
              pcmflux2(ixy,1,ielem) = pcmflux2(ixy,1,ielem) + 
     $                                 rmass(ibin,ig) * fbot
              ! print *,ibin,pcmflux2(ixy,1,ielem),fbot
              do k=2,NZ
               pcmflux2(ixy,k,ielem) = pcmflux2(ixy,k,ielem) + 
     $            rmass(ibin,ig) *( (vertdifu(k) + vertadvu(k)) *
     $            ((1-uc)*pc2(ixy,k-1,ibin,ielem) + uc*cvert(k-1)) 
     $            -  (vertdifd(k) + vertadvd(k)) *
     $            ((1-uc)*pc2(ixy,k,ibin,ielem) + uc*cvert(k)) ) 

              enddo
c
            endif
            if(itype(ielem) .eq. I_COREMASS .or.
     $         itype(ielem) .eq. I_GROWCORE      ) then
c
c             Core mass concentration
c
              pcmflux2(ixy,1,ielem) = pcmflux2(ixy,1,ielem) + fbot
              do k=2,NZ
               pcmflux2(ixy,k,ielem) = pcmflux2(ixy,k,ielem) + 
     $            (vertdifu(k) + vertadvu(k)) *
     $            ((1-uc)*pc2(ixy,k-1,ibin,ielem) + uc*cvert(k-1)) 
     $            -  (vertdifd(k) + vertadvd(k)) *
     $            ((1-uc)*pc2(ixy,k,ibin,ielem) + uc*cvert(k)) 

              enddo
c
            endif
c
c             (Not calculated for second moment of core mass distribution)
c
c  Update particle concentrations.
c
            do k = nlowlim,nuplim
              pc2(ixy,k,ibin,ielem) = cvert(k)
c!            if(pc2(ixy,k,ibin,ielem) .lt. 0. .and. ielem.eq.2 
c!   $            .and. abs(pc2(ixy,k,ibin,ielem)) .gt. FEW_PC   )
c!   $          write(*,9) 'Negative pc(cloud) from vtrans',
c!   $             pc2(ixy,k,ibin,ielem),k,ibin
c!            if(ielem.eq.1 .and. pc2(ixy,k,ibin,ielem).lt.0.)
c!   $         write(*,10) 'Negative pc(tholin) from vtrans',
c!   $             pc2(ixy,k,ibin,ielem),k,ibin,itime,time/8.64d4
              ixyz = k
              call smallconc(ibin,ielem)
            enddo

          enddo
        enddo

c=======================================================================
c
c  Treat vertical transport of gases 
c  (for comments, see above treatment of particles).
c
c
c        do k = 1,NZP1
c          vtrans(k) = w2(ixy,k)
c        enddo
c
c        itbnd = itbnd_gc
c        ibbnd = ibbnd_gc
c
c       do igas = 1,NGAS
c
c          do k = 1,NZ
c            gcmflux2(ixy,k,igas) = 0.
c          enddo
c
c          ftop = ftopgas(ixy,igas)
c          fbot = fbotgas(ixy,igas) !+ fbotevap(igas) -- Turn off puddle
c
c          do k = 1,NZP1
c            vtrans(k) = w2(ixy,k)
c            if( k .le. NZ )then
c              cvert(k) = gc2(ixy,k,igas)
c              if( NXY. eq. 1 ) divcor(k) = cvert(k) * drho_dt(k)
c            endif
c          enddo
c
c  Sensitivity Test:  Use normal eddy diffusion coefficient for particles,
c  but 10x for gas.  Diffusion coefficient initialized at normal 
c  value in <initatm>. 
c
c      Change dkz for gas
c         do k = 1,NZ
c           dkz3(k) = 10.d0 * dkz3(k)
c         enddo
c 
c          call vertadv
c          call vertdif
c          call versol
c
c            if( test_mass_cons .and. itime.eq.ibtime+1 ) 
c     $        write(*,*) ftop,fbot,igas
c
c      Change dkz back for particles
c         do k = 1,NZ
c           dkz3(k) = dkz3(k) / 10.d0
c         enddo
c 
c
c   Calculate gas mass flux
c
c              gcmflux2(ixy,1,igas) = gcmflux2(ixy,1,igas) + fbot
c              do k=2,NZ
c               gcmflux2(ixy,k,igas) = gcmflux2(ixy,k,igas) + 
c     $             (vertdifu(k) + vertadvu(k)) *
c     $            ((1-uc)*gc2(ixy,k-1,igas) + uc*cvert(k-1)) 
c     $            -  (vertdifd(k) + vertadvd(k)) *
c     $            ((1-uc)*gc2(ixy,k,igas) + uc*cvert(k)) 
c              enddo
c
c           ! Remove evaporated liquids from puddle reservoir
c              puddle(igas) = puddle(igas) - gcmflux2(ixy,1,igas)*dtime
c
c          do k = nlowlim,nuplim
c            gc2(ixy,k,igas) = cvert(k)
cc      if(k.eq.1) then
cc        write(*,*) 'gc from vertical',igas,gc2(ixy,k,igas),
cc   $               vertdifd(k),vertadvd(k),fbot
cc      endif
c          enddo
c
c        enddo
c
c=======================================================================
c
c  Treat vertical transport of potential temperature 
c  (for comments, see above treatment of particles).
c
c        itbnd = itbnd_ptc
c        ibbnd = ibbnd_ptc
 
c        ftop = 0.
c        fbot = 0.

c       do k = 1,NZP1
c         vtrans(k) = w2(ixy,k)
c         if( k .le. NZ )then
c           cvert(k) = ptc2(ixy,k)
c           if( NXY. eq. 1 ) divcor(k) = cvert(k) * drho_dt(k)
c         endif
c       enddo

c       call vertadv
c       call vertdif
c       call versol

c       do k = nlowlim,nuplim
c         ptc2(ixy,k) = cvert(k)
c       enddo
c
c=======================================================================
c
      enddo   ! <ixy=1,NXY>
c
c=======================================================================
c
c
c  Return to caller with new particle, gas, and potential temperature
c  concentrations and air density.
c
c     print *,'VERTICAL',pc3(6,1,1),dtime
      return
      end
