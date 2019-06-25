#include <misc.h>
#include <preproc.h>

module QSatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: QSatMod
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: QSat
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! Modified for Titan by ajf, 3/10/06
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: QSat
!
! !INTERFACE:
  subroutine QSat (T, p, es, esdT, qs, qsdT)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  CRC Handbook for methane saturation
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    use shr_const_mod, only:  SHR_CONST_RWV, SHR_CONST_LATVAP, & 
     SHR_CONST_ZVIR
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
    real(r8), intent(out) :: es       ! vapor pressure (pa)
    real(r8), intent(out) :: esdT     ! d(es)/d(T)
    real(r8), intent(out) :: qs       ! humidity (kg/kg)
    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8) :: T0,x,x0,y,xlor,ec,c,rup,rdwn,es0
    real(r8) :: eovp,vp,eps_loc
!
! For methane vapor, 65 K to 270 K
!-----------------------------------------------------------------------


    ec=0.0183156                    !exp(-4.)
    es0=4.520e6                     !Pa, vapor pressure at 190 K over liquid
    eps_loc=1.0/(1.0+SHR_CONST_ZVIR)
    xlor=SHR_CONST_LATVAP/SHR_CONST_RWV
    T0=190.0
    x0=xlor/T0
    x= xlor/T-x0
    if ( x < 4.0) then 
       y=x
       c=1.0
    endif
    if (x >= 4.0 .and. x <= 8.0) then 
       y=x-4.0
       c=ec
    endif
    if (x > 8.0 .and. x <= 12.0) then 
       y=x-8.0
       c=ec*ec
    endif
    if (x > 12.0) then  !lowest vapor pressure corresponds to T=65K
       y=4.0
       c=ec*ec
    endif
       rup= 1680.+y*(-840.+y*(180.+y*(-20.+y)))
       rdwn=1680.+y*( 840.+y*(180.+y*( 20.+y)))
       es=es0*c*rup/rdwn             !rup/rdwn is Pade' approx for e(-y)
       esdT=es/T*(x+x0)
     eovp=es/p
     eovp=min(eovp,1.0)  !prevent vp, and therefore qs, from becoming  < 0
     vp    = 1.0   / (1. - (1.-eps_loc)*eovp)

    qs    = eps_loc*vp*eovp             ! kg/kg
    qsdT  = esdT/es  * vp * qs         ! 1 / K

  end subroutine QSat

end module QSatMod
