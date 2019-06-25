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
c    ibin,ielem
c
c
c  Include global constants and variables.
c
      include 'globaer.h'
c
c
c  Declare local variables
c
      dimension drho_dt_flux_spec(NZ), drho_dt_fixed_conc(NZ)
c
c
c  Define formats
c
    3 format(i3,3x,6e13.4)
    4 format(i3,3x,i3,3x,i3,3x,6f9.2)
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter vertical'
c 
c=======================================================================
c
c
!     open(unit=99, file='pcold.txt', status='unknown') !EJL
!	  write(99,*) pc2(1,:,:,1)
!	  close(unit=99)
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
c  First calculate change in density for the case when boundary
c  fluxes are zero
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
c        vtrans(1) = 0.
c        vtrans(NZP1) = 0.
c
c
c  Since diffusion does not occur for air density, set diffusion
c  fluxes to zero.
c
c        do k = 1,NZP1
c          vertdifu(k) = 0.
c          vertdifd(k) = 0.
c        enddo
c
c
c  Calculate particle transport rates due to vertical advection
c  (note: no diffusion for air density), and solve for air density
c  at end of time step.
c
c        call vertadv
c        call versol

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
c            drho_dt_flux_spec(k) = ( rhoa2(ixy,k) - cvert(k) ) / 
c     $                   ( rhoa2(ixy,k) * dtime )
c          enddo

c        else
c
c
c  Update air density when not in 1D.
c
c          do k = 1,NZ
c            rhoa2(ixy,k) = cvert(k)
c          enddo

c        endif
c
c
c  Next, calculate change in density for the case when the
c  concentration just beyond the boundary is fixed
c
c        itbnd = I_FIXED_CONC
c        ibbnd = I_FIXED_CONC

c        do k = 1,NZP1
c          vtrans(k) = w2(ixy,k)
c        enddo
c        do k = 1,NZ
c          cvert(k) = rhoa2(ixy,k)
c        enddo

c        cvert_tbnd = cvert(NZ)
c        cvert_bbnd = cvert(1)

c        call vertadv
c        call versol

c        if( NXY .eq. 1 )then
c          do k = 1,NZ
c            drho_dt_fixed_conc(k) = ( rhoa2(ixy,k) - cvert(k) ) /
c     $                   ( rhoa2(ixy,k) * dtime )
c          enddo
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

          do ibin = 1,NBIN          ! Loop over particle mass bins
c
c
c  Fluxes across top and bottom of model
c
            ftop = ftoppart(ixy,ibin,ielem)
            fbot = fbotpart(ixy,ibin,ielem)
c
c
c  Store temporary (work) values in <cvert>, 
c  evaluate vertical velocities at layer boundaries, and
c  when 1D, assign divergence correction term
c

            do k = 1,NZP1
            w2(ixy,k)=0.0
              vtrans(k) = w2(ixy,k)
!              if (w2(ixy,k) .ne. 0.0) then
!                write(LUNOPRT, *) "ERROR: k=", k, " w2 != 0, w2=",
!     $            w2(ixy,k)
!              endif
              vtrans(k) = vtrans(k) - vf(k,ibin,ig)
            enddo
            if( ibbnd .eq. I_FLUX_SPEC ) vtrans(1) = 0.
            if( itbnd .eq. I_FLUX_SPEC ) vtrans(NZP1) = 0.
            if( itbnd .eq. I_FIXED_CONC )
     $        cvert_tbnd = pc_topbnd(ixy,ibin,ielem)
            if( ibbnd .eq. I_FIXED_CONC )
     $        cvert_bbnd = pc_botbnd(ixy,ibin,ielem)
            do k = 1,NZ
              cvert(k) = pc2(ixy,k,ibin,ielem)
c            if( NXY. eq. 1 ) 
c     $        divcor(k) = cvert(k) * drho_dt_flux_spec(k)
            enddo
c            if( NXY .eq. 1 ) then
c              if( ibbnd .eq. I_FIXED_CONC )
c     $          divcor(1) = cvert(1) * drho_dt_fixed_conc(1)
c              if( itbnd .eq. I_FIXED_CONC )
c     $          divcor(NZ) = cvert(NZ) * drho_dt_fixed_conc(NZ)
c            endif
c
!!      if (ibin .eq. 20) then
!!      open(unit=99, file='vert.txt', status='unknown') !EJL
!!      write(99,*) vtrans
!!      write(99,*) ''
!!      write(99,*) cvert
!!      close(unit=99)
!!      endif
c
c  Calculate particle transport rates due to vertical advection
c  and vertical diffusion, and solve for concentrations at end of time step.
c
            call vertadv

c            call vertdif
!            call versol
!            write(6,*) 'call verpsolve', ielem, ibin, ig
            call verpsolve(ibin,ielem)
!	  open(unit=99, file='vert2.txt', status='unknown') !EJL
!	  write(99,*) cvert
!	  close(unit=99)
c
c
c  Update particle concentrations.
c
            do k = 1,NZ
              pc2(ixy,k,ibin,ielem) = cvert(k)
            enddo

          enddo  !end  bin
        enddo   !end  elem
c
c=======================================================================
c
c  Treat vertical transport of gases 
c  (for comments, see above treatment of particles).
c
c
c        itbnd = itbnd_gc
c        ibbnd = ibbnd_gc

c        do k = 1,NZP1
c          vtrans(k) = w2(ixy,k)
c        enddo
c        if( ibbnd .eq. I_FLUX_SPEC ) vtrans(1) = 0.
c        if( itbnd .eq. I_FLUX_SPEC ) vtrans(NZP1) = 0.
c
c        do igas = 1,NGAS
c
c          ftop = ftopgas(ixy,igas)
c          fbot = fbotgas(ixy,igas)
c          if( itbnd .eq. I_FIXED_CONC )
c     $      cvert_tbnd = gc_topbnd(ixy,igas)
c          if( ibbnd .eq. I_FIXED_CONC )
c     $      cvert_bbnd = gc_botbnd(ixy,igas)
c
c          do k = 1,NZ
c            cvert(k) = gc2(ixy,k,igas)
c            if( NXY. eq. 1 ) 
c     $       divcor(k) = cvert(k) * drho_dt_flux_spec(k)
c          enddo
c            if( NXY .eq. 1 ) then
c            if( ibbnd .eq. I_FIXED_CONC )
c     $        divcor(1) = cvert(1) * drho_dt_fixed_conc(1)
c            if( itbnd .eq. I_FIXED_CONC )
c     $        divcor(NZ) = cvert(NZ) * drho_dt_fixed_conc(NZ)
c          endif
c
c          call vertadv
c          call vertdif
c          call versol
c
c          do k = 1,NZ
c            gc2(ixy,k,igas) = cvert(k)
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
c 
c        ftop = 0.
c        fbot = 0.
c
c        do k = 1,NZP1
c          vtrans(k) = w2(ixy,k)
c        enddo
c        if( ibbnd .eq. I_FLUX_SPEC ) vtrans(1) = 0.
c        if( itbnd .eq. I_FLUX_SPEC ) vtrans(NZP1) = 0.
c        if( itbnd .eq. I_FIXED_CONC ) cvert_tbnd = ptc_topbnd(ixy)
c        if( ibbnd .eq. I_FIXED_CONC ) cvert_bbnd = ptc_botbnd(ixy)
c        do k = 1,NZ
c          cvert(k) = ptc2(ixy,k)
c          if(( NXY. eq. 1 ) .and. ( is_one_dim )) 
c     $      divcor(k) = cvert(k) * drho_dt_flux_spec(k)
c        enddo
c        if(( NXY. eq. 1 ) .and. ( is_one_dim )) then
c          if( ibbnd .eq. I_FIXED_CONC )
c     $      divcor(1) = cvert(1) * drho_dt_fixed_conc(1)
c          if( itbnd .eq. I_FIXED_CONC )
c     $      divcor(NZ) = cvert(NZ) * drho_dt_fixed_conc(NZ)
c        endif
c
c        call vertadv
c        call vertdif
c        call versol
c
c        do k = 1,NZ
c          ptc2(ixy,k) = cvert(k)
c        enddo
c
c=======================================================================
c
      enddo   ! <ixy=1,NXY>
!	  open(unit=99, file='pcnew.txt', status='unknown') !EJL
!	  write(99,*) pc2(1,:,:,1)
!	  close(unit=99)
c
c=======================================================================
c
c
c  Return to caller with new particle, gas, and potential temperature
c  concentrations and air density.
c
      return
      end
