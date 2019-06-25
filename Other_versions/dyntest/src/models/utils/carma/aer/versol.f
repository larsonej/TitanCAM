      subroutine versol
c
c
c  @(#) versol.f  Jensen  Dec-1996
c  This routine solves the vertical transport equation.
c  <cvert> is temporary storage for concentrations (particles,
c  gas, potential temperature) being transported.
c  New values of <cvert> are calculated.
c
c  Modified  Sep-1997  (McKie)
c  Remove <ixy> from arg list
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
      dimension al(NZ), bl(NZ), dl(NZ), ul(NZ),
     $          el(NZ), fl(NZ), ctempl(NZ), ctempu(NZ)
c 
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter versol'
c
c
c  Determine whether transport should be solved explicitly (uc=0)
c  or implicitly (uc=1).
c
      uc = 0.
      do k = 1,NZ
        cour = dz2(ixy,k)/dtime -
     $   ( vertdifu(k+1) + vertdifd(k) + vertadvu(k+1) + vertadvd(k) )
        if( cour .lt. 0. .and. uc .ne. 1. )then
          uc = 1.0
c          write(LUNOPRT,'(a,i3,7(1x,1pe8.1))')
c     $      'in versol: k dz/dt vdifd vdifu vadvd vadvu cour uc = ',
c     $       k, dz2(ixy,k)/dtime, vertdifd(k), vertdifu(k+1),
c     $       vertadvd(k), vertadvu(k+1), cour, uc
        endif
      enddo
c
c
c  Store concentrations in local variables (shifted up and down
c  a vertical level).
c
      do k = 2,NZ
        ctempl(k) = cvert(k-1)
        ctempu(k-1) = cvert(k)
      enddo

      if( ibbnd .eq. I_FIXED_CONC ) then
        ctempl(1) = cvert_bbnd
      else
        ctempl(1) = 0.
      endif

      if( itbnd .eq. I_FIXED_CONC ) then
        ctempu(NZ) = cvert_tbnd
      else
       ctempu(NZ) = 0.
      endif
c
c
c  Calculate coefficients of the transport equation:
c    al(k)c(k+1) + bl(k)c(k) + ul(k)c(k-1) = dl(k)
c
      do k = 1,NZ
        al(k) = uc * ( vertdifd(k+1) + vertadvd(k+1) )
        bl(k) = -( uc*(vertdifd(k)+vertdifu(k+1)+
     $                 vertadvd(k)+vertadvu(k+1))
     $             + dz2(ixy,k)/dtime )
        ul(k) = uc * ( vertdifu(k) + vertadvu(k) )
        dl(k) = cvert(k) *
     $     ( (1.-uc)*(vertdifd(k)+vertdifu(k+1)+
     $                vertadvd(k)+vertadvu(k+1))
     $     - dz2(ixy,k)/dtime ) -
     $     (1.-uc) * ( (vertdifu(k)+vertadvu(k))*ctempl(k) +
     $     (vertdifd(k+1)+vertadvd(k+1))*ctempu(k) )
!     $     - divcor(k) * dz2(ixy,k)
      enddo
c
c
c  Boundary fluxes: <ftop> is the downward flux across the
c  upper boundary; <fbot> is the upward flux across the
c  lower boundary.
c
       if( igridv .eq. I_SIG ) then
         if( itbnd .eq. I_FLUX_SPEC ) dl(1) = dl(1) - ftop
         if( ibbnd .eq. I_FLUX_SPEC ) dl(NZ) = dl(NZ) - fbot
       else
         if( itbnd .eq. I_FLUX_SPEC ) dl(NZ) = dl(NZ) - ftop
         if( ibbnd .eq. I_FLUX_SPEC ) dl(1) = dl(1) - fbot
       endif
c
c
c  Calculate recursion relations.
c
      el(1) = dl(1)/bl(1)
      fl(1) = al(1)/bl(1)
      do k = 2,NZ
        denom = bl(k) - ul(k) * fl(k-1)
        el(k) = ( dl(k) - ul(k)*el(k-1) ) / denom
        fl(k) = al(k) / denom
      enddo
c
c
c  Calculate new concentrations.
c
      cvert(NZ) = el(NZ)
      do k = NZ-1,1,-1
        cvert(k) = el(k) - fl(k)*cvert(k+1)
      enddo
c
c
c  Return to caller with new concentrations.
c
      return
      end
