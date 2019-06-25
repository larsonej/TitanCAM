!===============================================================================
! CVS $Id: shr_const_mod.F90 61 2008-03-18 22:05:15Z cam_titan $
! CVS $Source: /home/cvsroot/cam_titan/src/models/csm_share/shr/shr_const_mod.F90,v $
! CVS $Name:  $
!===============================================================================

MODULE shr_const_mod

   use shr_kind_mod, only: SHR_KIND_R8

   !----------------------------------------------------------------------------
   ! physical constants (all data public)
   !----------------------------------------------------------------------------
   public
   real(SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_SHR_KIND_R8  ! pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_CDAY   = 86400.0_SHR_KIND_R8  ! sec in Earth calendar day 
   real(SHR_KIND_R8),parameter :: SHR_CONST_SDAY   = 86164.0_SHR_KIND_R8      ! sec in siderial day ~ sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/1.37808e6_SHR_KIND_R8 ! titan rot ~ rad/sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 2.57500e6_SHR_KIND_R8    ! radius of titan ~ m
   real(SHR_KIND_R8),parameter :: SHR_CONST_G      = 1.35000_SHR_KIND_R8      ! acceleration of gravity ~ m/s^2
   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 146640.0_SHR_KIND_R8     ! standard pressure ~ pascals

   real(SHR_KIND_R8),parameter :: SHR_CONST_STEBOL = 5.67e-8_SHR_KIND_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(SHR_KIND_R8),parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_SHR_KIND_R8  ! Boltzmann's constant ~ J/K/molecule
   real(SHR_KIND_R8),parameter :: SHR_CONST_AVOGAD = 6.02214e26_SHR_KIND_R8   ! Avogadro's number ~ molecules/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWDAIR = 28.000_SHR_KIND_R8       ! molecular weight dry air ~ kg/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWWV   = 16.000_SHR_KIND_R8       ! molecular weight methane vapor
   real(SHR_KIND_R8),parameter :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! CH4 - vapor gas constant J/K/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_ZVIR   = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0_SHR_KIND_R8   ! RWV/RDAIR - 1.0
   real(SHR_KIND_R8),parameter :: SHR_CONST_KARMAN = 0.4_SHR_KIND_R8          ! Von Karman constant
 
   real(SHR_KIND_R8),parameter :: SHR_CONST_TKFRZ  = 90.67_SHR_KIND_R8       ! freezing T of CH4 ~ K (intentionally made == to TKTRIP)
   real(SHR_KIND_R8),parameter :: SHR_CONST_TKTRIP = 90.67_SHR_KIND_R8       ! triple point of fresh water ~ K

   real(SHR_KIND_R8),parameter :: SHR_CONST_RHODAIR=SHR_CONST_PSTD/ &
     (SHR_CONST_RDAIR*SHR_CONST_TKFRZ)         ! density of dry air at STP   ~ kg/m^3  
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOFW  = 4.400e2_SHR_KIND_R8      ! density of liquid methane kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOSW  = 4.400e2_SHR_KIND_R8      ! density of   "       "   kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOICE = 0.440e3_SHR_KIND_R8      ! ~density of ch4 ice   ~ kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = 1.0714e3_SHR_KIND_R8     ! specific heat of dry air ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPFW   = 3.480e3_SHR_KIND_R8      ! specific heat of liquid ch4 ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPSW   = 3.480e3_SHR_KIND_R8      ! specific heat of   "     " ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPWV   = 2.162e3_SHR_KIND_R8      ! specific heat of ch4 vap ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPICE  = 3.480e3_SHR_KIND_R8      ! approx specific heat of ch4 ice ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_LATICE = 5.875e4_SHR_KIND_R8      ! latent heat of fusion ch4 ~ J/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_LATVAP = 5.466e5_SHR_KIND_R8      ! latent heat of evaporation ch4 ~ J/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_LATSUB = SHR_CONST_LATICE + SHR_CONST_LATVAP ! latent heat of sublimation ~ J/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_OCN_REF_SAL = 34.7_SHR_KIND_R8    ! ocn ref salinity (psu)
   real(SHR_KIND_R8),parameter :: SHR_CONST_ICE_REF_SAL =  4.0_SHR_KIND_R8    ! ice ref salinity (psu)

   real(SHR_KIND_R8),parameter :: SHR_CONST_SPVAL       = 1.0e30_SHR_KIND_R8  ! special missing value

   ! FAO:  Additional parameters for Titan
   real(SHR_KIND_R8),parameter :: TITAN_DAY_RATIO = 15.9454_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: TITAN_YEAR_RATIO = 29.46_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: TITAN_ECCENTRICITY = 0.0542_SHR_KIND_R8
   !real(SHR_KIND_R8),parameter :: TITAN_ECCENTRICITY = 0.056_SHR_KIND_R8  ! Tokano
   real(SHR_KIND_R8),parameter :: TITAN_OBLIQUITY = 0.4665_SHR_KIND_R8  ! = 26.73 degrees
   real(SHR_KIND_R8),parameter :: TITAN_EMISSIVITY = 0.86_SHR_KIND_R8

END MODULE shr_const_mod
