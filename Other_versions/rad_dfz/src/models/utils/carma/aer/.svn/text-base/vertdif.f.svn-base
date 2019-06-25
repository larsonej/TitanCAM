       subroutine vertdif
c
c
c  @(#) vertdif.f  Jensen  Dec-1996
c  This routine calculates vertrical transport rates.
c  Currently treats diffusion only.
c  Not necessarily generalized for irregular grid.
c
c  <vertdifu(k)> is upward vertical transport rate into level k from level k-1.
c  <vertdifd(k)> is downward vertical transport rate into level k-1 from level k.
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
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter vertdif'
c
c
c  Set some constants
c
      nzm1 = max( 1, NZ-1 )
      itwo = min( 2, NZ   )
c
c
c  Loop over vertical levels.
c
      do k = 2, NZ

        dz_avg = dz2(ixy,k)                            ! layer thickness
c
c  Check the vertical coordinate
c  
        if( igridv .eq. I_CART ) then
          rhofact = log(  rhoa2(ixy,k)/rhoa2(ixy,k-1)
     $                  * zmet2(ixy,k-1)/zmet2(ixy,k) )
          xex = rhoa2(ixy,k-1)/rhoa2(ixy,k) *
     $          zmet2(ixy,k)/zmet2(ixy,k-1)
          vertdifu(k) = ( rhofact * dkz2(ixy,k)/dz_avg ) /
     $                  ( 1. - xex )

          vertdifd(k) = vertdifu(k) * xex
c
c  ...else you're in sigma coordinates...
c
        elseif(( igridv .eq. I_SIG ) .or. ( igridv .eq. I_HYB)) then
          vertdifu(k) = dkz2(ixy,k)/dz_avg
          vertdifd(k) = dkz2(ixy,k)/dz_avg
c
c  ...else write an error (maybe redundant)...
c  
        else
          write(LUNOPRT,*) 'Invalid vertical grid type of ', igridv 
          call endcarma
        endif

      enddo
c
c
c  Fluxes at boundaries specified by user
c
      if( ibbnd .eq. I_FLUX_SPEC ) then
        vertdifu(1) = 0.
        vertdifd(1) = 0.
      endif

      if( itbnd .eq. I_FLUX_SPEC ) then
        vertdifu(NZ+1) = 0.
        vertdifd(NZ+1) = 0.
      endif
c
c
c  Diffusion across boundaries using fixed boundary concentration:
c
      if( ibbnd .eq. I_FIXED_CONC ) then
        dz_avg = dz2(ixy,1)                            ! layer thickness
        rhofact = log( rhoa2(ixy,itwo)/rhoa2(ixy,1) )
        ttheta = rhofact
        if( ttheta .ge. 0. ) then
          ttheta = min(ttheta,POWMAX)
        else
          ttheta = max(ttheta,-POWMAX)
        endif

        xex = exp(-ttheta)
        if( abs(ONE - xex) .lt. ALMOST_ZERO ) xex = ALMOST_ONE

        vertdifu(1) = ( rhofact * dkz2(ixy,1)/dz_avg ) /
     $                ( 1. - xex )
        vertdifd(1) = vertdifu(1) * xex
      endif

      if( itbnd .eq. I_FIXED_CONC ) then
        dz_avg = dz2(ixy,NZ)                            ! layer thickness
        rhofact = log( rhoa2(ixy,NZ)/rhoa2(ixy,nzm1) )
        ttheta = rhofact
        if( ttheta .ge. 0. ) then
          ttheta = min(ttheta,POWMAX)
        else
          ttheta = max(ttheta,-POWMAX)
        endif

        xex = exp(-ttheta)
        if( abs(ONE - xex) .lt. ALMOST_ZERO ) xex = ALMOST_ONE

        vertdifu(NZ+1) = ( rhofact * dkz2(ixy,NZ+1)/dz_avg ) /
     $                   ( 1. - xex )
        vertdifd(NZ+1) = vertdifu(NZ+1) * xex
      endif
c
c
c  Return to caller with vertical diffusion rates.
c
      return
      end
