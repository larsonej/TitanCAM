       subroutine vaporp
c
c
c  @(#) vaporp.f  Ackerman  Dec-1995
c  This routine calculates the vapor pressure for all gases 
c  over the entire spatial grid:
c
c  <pvapl> and <pvapi> are vapor pressures in units of [dyne/cm^2]
c
c  Uses temperature <t> as input.
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c  (actually, only <ixyz> is defined -- the others are meaningless)
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
c  Define coefficients in Buck's formulation for saturation vapor pressures
c
      parameter( BAI = 6.1115e3 )
      parameter( BBI = 23.036 )              
      parameter( BCI = 279.82 )
      parameter( BDI = 333.7 )

      parameter( BAL = 6.1121e3 )           
      parameter( BBL = 18.729 )
      parameter( BCL = 257.87 )
      parameter( BDL = 227.3 )
c
c
c  Define formats
c
    1 format('T = ',1pe12.3,a,i6,a,1pe11.3)
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter vaporp'
c
c-------------------------------------------------------------------------------
c
c
c  Loop over all gases.
c
      do igas = 1, NGAS
c
c
c  Check for expected gas index
c
          if( igas .eq. 1 )then
c
c
c  Saturation vapor pressure over liquid water and water ice
c  (from Buck [J. Atmos. Sci., 20, 1527, 1981])
c
           tt = t3(ixyz) - 273.16

           pvapl3(ixyz,igas) = BAL * 
     $                         exp( (BBL - tt/BDL)*tt / (tt + BCL) )

           pvapi3(ixyz,igas) = BAI * 
     $                         exp( (BBI - tt/BDI)*tt / (tt + BCI) )
c
c  Check to see whether temperature is ouside range of validity
c  for parameterizations - changed for Titan
c
!           if( pvapl3(ixyz,igas) .le. 1.e-10 ) then
!            write(*,1) t3(ixyz), ' too small for ixyz = ',
!     $                       ixyz, ' time = ',time
!            call endcarma
!           endif
c
c
c  Report unexpected gas index
c
          else
            write(LUNOPRT,'(/,a)') 'invalid <igas> in vaporp.f'
            call endcarma
          endif
 
      enddo
c
c
c  Return to caller with vapor pressures evaluated.
c
      return
      end
