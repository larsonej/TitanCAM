subroutine d_ave(frac_day, day_in_year, clat, clon, diurnal_ave, ha, ncol)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute diurnally averaged cosine of solar zenith angle for chemistry computations.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Kiehl/E.Wilson
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_orb_mod
   use shr_const_mod, only: shr_const_pi

   implicit none

#include <crdcon.h>
#include <comsol.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of positions
   real(r8), intent(in) :: frac_day 
   real(r8), intent(in) :: day_in_year
   real(r8), intent(in) :: clat(ncol)          ! Current centered latitude (radians)
   real(r8), intent(in) :: clon(ncol)          ! Centered longitude (radians)
!
! Output arguments
!
   real(r8), intent(out) :: diurnal_ave(ncol)       ! Cosine solar zenith angle
   real(r8), intent(out) :: ha(ncol)           ! Hour angle
!
!---------------------------Local variables-----------------------------
!
   integer i         ! Position loop index
   real(r8) delta    ! Solar declination angle  in radians
   real(r8) eccf     ! Earth orbit eccentricity factor
   real(r8) coszrs
   real(r8), parameter :: pi = shr_const_pi
!
!-----------------------------------------------------------------------
!

   ! need to modify eccen, mvelpp, lambm0, and obliqr for titan
   call shr_orb_decl (day_in_year  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf      )
!
! Compute local cosine solar zenith angle,
!
   do i=1,ncol
      ha(i) = -tan(clat(i))*tan(delta)
      if (ha(i).gt.1.0d0) then
         ha(i) = 0
      elseif (ha(i).lt.-1.0) then
         ha(i) = pi
      else
         ha(i) = acos(ha(i))
      endif
      diurnal_ave(i) = sin(clat(i))*sin(delta) +    &
                       cos(clat(i))*cos(delta)*sin(ha(i))/ha(i)
   end do
!
   return
end subroutine d_ave
