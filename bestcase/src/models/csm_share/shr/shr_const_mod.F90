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
   
   real(SHR_KIND_R8),parameter :: rad_sun = 6.96265e10/1.496e13 !solar radius, AU


#ifdef _TITAN

#ifndef _MV10   
   
   real(SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_SHR_KIND_R8  ! pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_CDAY   = 86400.0_SHR_KIND_R8  ! sec in Earth calendar day 
   real(SHR_KIND_R8),parameter :: SHR_CONST_SDAY   = 86164.0_SHR_KIND_R8      ! sec in siderial day ~ sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/1.37808e6_SHR_KIND_R8 ! titan rot ~ rad/sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 2.57500e6_SHR_KIND_R8    ! radius of titan ~ m
   real(SHR_KIND_R8),parameter :: SHR_CONST_G      = 1.35000_SHR_KIND_R8      ! acceleration of gravity ~ m/s^2
   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 146640.0_SHR_KIND_R8     ! standard pressure ~ pascals
   real(SHR_KIND_R8),parameter :: SHR_CONST_QH2    = 1.0e-3_SHR_KIND_R8      !mole fraction
   real(SHR_KIND_R8),parameter :: SHR_CONST_QHE    = 0.0_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QN2    = 1.0_SHR_KIND_R8 !Recomputed in process_cias for Titan

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
   real(SHR_KIND_R8),parameter :: PLANET_DAY_RATIO = 15.9454_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_YEAR_RATIO = 29.46_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_ECCENTRICITY = 0.0542_SHR_KIND_R8
   !real(SHR_KIND_R8),parameter :: PLANET_ECCENTRICITY = 0.056_SHR_KIND_R8  ! Tokano
   real(SHR_KIND_R8),parameter :: PLANET_OBLIQUITY = 0.4665_SHR_KIND_R8  ! = 26.73 degrees
   real(SHR_KIND_R8),parameter :: PLANET_MVELP = 90.643_SHR_KIND_R8  !as of 2000, JD 245 1800.5
   real(SHR_KIND_R8),parameter :: PLANET_EMISSIVITY = 0.86_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: shr_const_au = 9.539_SHR_KIND_R8   
   real(SHR_KIND_R8),parameter :: sun_size=shr_const_pi*(rad_sun/shr_const_au)*(rad_sun/shr_const_au) !angular size of sun at planet 
   real(SHR_KIND_R8),parameter :: shr_const_intrnlflux = 0.0  ! W m^-2      

#else ! _MV10 defined, do MV10 experiment with the following parameters
      !  Note _TITAN must also be defined to do this experiment

   real(SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_SHR_KIND_R8  ! pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_CDAY   = 86400.0_SHR_KIND_R8  ! sec in Earth calendar day 
   real(SHR_KIND_R8),parameter :: SHR_CONST_SDAY   = 86164.0_SHR_KIND_R8      ! sec in siderial day ~ sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/8.97600e4_SHR_KIND_R8 ! mv10 rot ~ rad/sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 2.80000e5_SHR_KIND_R8    ! radius of mv10 exper ~ m
   real(SHR_KIND_R8),parameter :: SHR_CONST_G      = 9.80000_SHR_KIND_R8      ! acceleration of gravity ~ m/s^2
   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 100000.0_SHR_KIND_R8     ! standard pressure ~ pascals
   real(SHR_KIND_R8),parameter :: SHR_CONST_QH2    = 1.0e-3_SHR_KIND_R8      !mole fraction
   real(SHR_KIND_R8),parameter :: SHR_CONST_QHE    = 0.0_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QN2    = 1.0_SHR_KIND_R8 !Recomputed in process_cias for Titan

   real(SHR_KIND_R8),parameter :: SHR_CONST_STEBOL = 5.67e-8_SHR_KIND_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(SHR_KIND_R8),parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_SHR_KIND_R8  ! Boltzmann's constant ~ J/K/molecule
   real(SHR_KIND_R8),parameter :: SHR_CONST_AVOGAD = 6.02214e26_SHR_KIND_R8   ! Avogadro's number ~ molecules/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWDAIR = 28.000_SHR_KIND_R8       ! molecular weight dry air ~ kg/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWWV   = 18.000_SHR_KIND_R8       ! molecular weight water vapor
   real(SHR_KIND_R8),parameter :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! CH4 - vapor gas constant J/K/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_ZVIR   = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0_SHR_KIND_R8   ! RWV/RDAIR - 1.0
   real(SHR_KIND_R8),parameter :: SHR_CONST_KARMAN = 0.4_SHR_KIND_R8          ! Von Karman constant
 
   real(SHR_KIND_R8),parameter :: SHR_CONST_TKFRZ  = 90.67_SHR_KIND_R8       ! freezing T of CH4 ~ K (intentionally made == to TKTRIP)
   real(SHR_KIND_R8),parameter :: SHR_CONST_TKTRIP = 90.67_SHR_KIND_R8       ! triple point of fresh water ~ K

   real(SHR_KIND_R8),parameter :: SHR_CONST_RHODAIR=SHR_CONST_PSTD/ &
     (SHR_CONST_RDAIR*SHR_CONST_TKFRZ)         ! density of dry air at STP   ~ kg/m^3  
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOFW   = 1.000e3_SHR_KIND_R8      ! density of fresh water     ~ kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOSW   = 1.026e3_SHR_KIND_R8      ! density of sea water       ~ kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOICE  = 0.917e3_SHR_KIND_R8      ! density of ice             ~ kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR  = 1.00464e3_SHR_KIND_R8    ! specific heat of dry air   ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPWV    = 1.810e3_SHR_KIND_R8      ! specific heat of water vap ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPVIR   = (SHR_CONST_CPWV/SHR_CONST_CPDAIR)-1.0_SHR_KIND_R8 ! CPWV/CPDAIR - 1.0
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPFW    = 4.188e3_SHR_KIND_R8      ! specific heat of fresh h2o ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPSW    = 3.996e3_SHR_KIND_R8      ! specific heat of sea h2o   ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPICE   = 2.11727e3_SHR_KIND_R8    ! specific heat of fresh ice ~ J/kg/K
   real(SHR_KIND_R8),parameter :: SHR_CONST_LATICE  = 3.337e5_SHR_KIND_R8      ! latent heat of fusion      ~ J/kg 
   real(SHR_KIND_R8),parameter :: SHR_CONST_LATVAP = 2.501e6_SHR_KIND_R8      ! latent heat of evaporation ch4 ~ J/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_LATSUB = SHR_CONST_LATICE + SHR_CONST_LATVAP ! latent heat of sublimation ~ J/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_OCN_REF_SAL = 34.7_SHR_KIND_R8    ! ocn ref salinity (psu)
   real(SHR_KIND_R8),parameter :: SHR_CONST_ICE_REF_SAL =  4.0_SHR_KIND_R8    ! ice ref salinity (psu)

   real(SHR_KIND_R8),parameter :: SHR_CONST_SPVAL       = 1.0e30_SHR_KIND_R8  ! special missing value

   ! FAO:  Additional parameters for Titan
   real(SHR_KIND_R8),parameter :: PLANET_DAY_RATIO = 15.9454_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_YEAR_RATIO = 29.46_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_ECCENTRICITY = 0.0542_SHR_KIND_R8
   !real(SHR_KIND_R8),parameter :: PLANET_ECCENTRICITY = 0.056_SHR_KIND_R8  ! Tokano
   real(SHR_KIND_R8),parameter :: PLANET_OBLIQUITY = 0.4665_SHR_KIND_R8  ! = 26.73 degrees
   real(SHR_KIND_R8),parameter :: PLANET_MVELP = 90.643_SHR_KIND_R8  !as of 2000, JD 245 1800.5
   real(SHR_KIND_R8),parameter :: PLANET_EMISSIVITY = 0.86_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: shr_const_au = 9.539_SHR_KIND_R8   
   real(SHR_KIND_R8),parameter :: sun_size=shr_const_pi*(rad_sun/shr_const_au)*(rad_sun/shr_const_au) !angular size of sun at planet 
   real(SHR_KIND_R8),parameter :: shr_const_intrnlflux = 0.0  ! W m^-2
#endif !end _MV10 block

#endif !end _TITAN block

#ifdef _JUPITER   
   
   real(SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_SHR_KIND_R8  ! pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_CDAY   = 86400.0_SHR_KIND_R8  ! sec in Earth calendar day 
   real(SHR_KIND_R8),parameter :: SHR_CONST_SDAY   = 86164.0_SHR_KIND_R8      ! sec in siderial day ~ sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/3.57273e4_SHR_KIND_R8 ! jupiter rot ~ rad/sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 6.9911e7_SHR_KIND_R8    ! mean radius of jupiter ~ m
   real(SHR_KIND_R8),parameter :: SHR_CONST_G      = 23.750_SHR_KIND_R8      ! acceleration of gravity ~ m/s^2
   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 200000.0_SHR_KIND_R8   ! standard pressure ~ pascals, nust match p0!
   real(SHR_KIND_R8),parameter :: SHR_CONST_QH2    = 0.864_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QHE    = 0.136_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QN2    = 1.0e-10_SHR_KIND_R8

   real(SHR_KIND_R8),parameter :: SHR_CONST_STEBOL = 5.67e-8_SHR_KIND_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(SHR_KIND_R8),parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_SHR_KIND_R8  ! Boltzmann's constant ~ J/K/molecule
   real(SHR_KIND_R8),parameter :: SHR_CONST_AVOGAD = 6.02214e26_SHR_KIND_R8   ! Avogadro's number ~ molecules/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWDAIR = 2.0*SHR_CONST_QH2+4.0*SHR_CONST_QHE       ! molecular weight dry air ~ kg/kmole, includes helium
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
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = 1.2654e4_SHR_KIND_R8     ! specific heat of normal hydrogen ~ J/kg/K
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

   ! FAO:  Additional parameters for Jupiter
   real(SHR_KIND_R8),parameter :: PLANET_DAY_RATIO = 0.4135_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_YEAR_RATIO = 11.86262_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_ECCENTRICITY = 0.0484_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_OBLIQUITY = 0.0545_SHR_KIND_R8  ! = 3.12 degrees
   real(SHR_KIND_R8),parameter :: PLANET_MVELP = 15.431_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_EMISSIVITY = 0.86_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: shr_const_au = 5.203_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: sun_size=shr_const_pi*(rad_sun/shr_const_au)*(rad_sun/shr_const_au) !angular size of sun at planet 
   real(SHR_KIND_R8),parameter :: shr_const_intrnlflux = 5.44  ! W m^-2
            
#endif

#ifdef _SATURN   
   
   real(SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_SHR_KIND_R8  ! pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_CDAY   = 86400.0_SHR_KIND_R8  ! sec in Earth calendar day 
   real(SHR_KIND_R8),parameter :: SHR_CONST_SDAY   = 86164.0_SHR_KIND_R8      ! sec in siderial day ~ sec
!!!   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/3.8340e4_SHR_KIND_R8 ! saturn rot ~ rad/sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/3.8053e4_SHR_KIND_R8 ! saturn rot ~ Sys w (Read et al)
   real(SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 6.0000e7_SHR_KIND_R8    ! ~ radius of saturn ~ m
   real(SHR_KIND_R8),parameter :: SHR_CONST_G      = 10.420_SHR_KIND_R8      ! acceleration of gravity ~ m/s^2
   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 200000.0_SHR_KIND_R8     ! standard pressure ~ pascals
   real(SHR_KIND_R8),parameter :: SHR_CONST_QH2    = 0.920_SHR_KIND_R8  
   real(SHR_KIND_R8),parameter :: SHR_CONST_QHE    = 0.080_SHR_KIND_R8 !From Leigh Fletcher, pers comm (Cassini results)
   real(SHR_KIND_R8),parameter :: SHR_CONST_QN2    = 1.0e-10_SHR_KIND_R8

   real(SHR_KIND_R8),parameter :: SHR_CONST_STEBOL = 5.67e-8_SHR_KIND_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(SHR_KIND_R8),parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_SHR_KIND_R8  ! Boltzmann's constant ~ J/K/molecule
   real(SHR_KIND_R8),parameter :: SHR_CONST_AVOGAD = 6.02214e26_SHR_KIND_R8   ! Avogadro's number ~ molecules/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWDAIR = 2.0*SHR_CONST_QH2+4.0*SHR_CONST_QHE       ! molecular weight dry air ~ kg/kmole, includes helium
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
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = 1.2654e4_SHR_KIND_R8     ! specific heat of normal hydrogen ~ J/kg/K
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

   ! FAO:  Additional parameters for Saturn
   real(SHR_KIND_R8),parameter :: PLANET_DAY_RATIO = 0.44375_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_YEAR_RATIO = 29.46_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_ECCENTRICITY = 0.0542_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_OBLIQUITY = 0.4665_SHR_KIND_R8  ! = 26.73 degrees
   real(SHR_KIND_R8),parameter :: PLANET_MVELP = 90.643_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_EMISSIVITY = 0.86_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: shr_const_au = 9.539_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: sun_size=shr_const_pi*(rad_sun/shr_const_au)*(rad_sun/shr_const_au) !angular size of sun at planet    
   real(SHR_KIND_R8),parameter :: shr_const_intrnlflux = 2.01  ! W m^-2
   
#endif

#ifdef _URANUS   
   
   real(SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_SHR_KIND_R8  ! pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_CDAY   = 86400.0_SHR_KIND_R8  ! sec in Earth calendar day 
   real(SHR_KIND_R8),parameter :: SHR_CONST_SDAY   = 86164.0_SHR_KIND_R8      ! sec in siderial day ~ sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/6.2064e4_SHR_KIND_R8 ! uranus rot ~ rad/sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 2.5559e7_SHR_KIND_R8    ! ~ radius of uranus ~ m (Hubbard, Neptune Book, Tab. 1)
   real(SHR_KIND_R8),parameter :: SHR_CONST_G      = 8.920_SHR_KIND_R8      ! acceleration of gravity ~ m/s^2
!!!   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 7.0e7_SHR_KIND_R8     ! standard pressure ~ pascals
   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 4.0e7_SHR_KIND_R8     ! standard pressure ~ pascals
!!!   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 1.0e7_SHR_KIND_R8     ! standard pressure ~ pascals
!!!   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 5.0e6_SHR_KIND_R8     ! standard pressure ~ pascals 
!!!   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 1.0e6_SHR_KIND_R8     ! standard pressure ~ pascals
!!!   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 2.0e5_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QH2    = 0.848_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QHE    = 0.152_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QN2    = 1.0e-10_SHR_KIND_R8

   real(SHR_KIND_R8),parameter :: SHR_CONST_STEBOL = 5.67e-8_SHR_KIND_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(SHR_KIND_R8),parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_SHR_KIND_R8  ! Boltzmann's constant ~ J/K/molecule
   real(SHR_KIND_R8),parameter :: SHR_CONST_AVOGAD = 6.02214e26_SHR_KIND_R8   ! Avogadro's number ~ molecules/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWDAIR = 2.0*SHR_CONST_QH2+4.0*SHR_CONST_QHE       ! molecular weight dry air ~ kg/kmole, includes helium
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
! The following expression for "frozen" specific heat approximates CP_H2 as 3.0*R, a compromise to true temperature/fp variation
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = (SHR_CONST_QH2*3.0+SHR_CONST_QHE*1.5)*SHR_CONST_RDAIR
!!!   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = 1.2654e4_SHR_KIND_R8     ! specific heat of normal hydrogen ~ J/kg/K
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

   ! FAO:  Additional parameters for Uranus
   real(SHR_KIND_R8),parameter :: PLANET_DAY_RATIO = 0.7183_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_YEAR_RATIO = 84.0119_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_ECCENTRICITY = 0.047_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_OBLIQUITY = 1.4329_SHR_KIND_R8  ! = 82.10 degrees (IAU)
   real(SHR_KIND_R8),parameter :: PLANET_MVELP = 169.440_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_EMISSIVITY = 0.86_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: shr_const_au = 19.182_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: sun_size=shr_const_pi*(rad_sun/shr_const_au)*(rad_sun/shr_const_au) !angular size of sun at planet 
   real(SHR_KIND_R8),parameter :: shr_const_intrnlflux = 0.042  ! W m^-2
     
#endif

#ifdef _NEPTUNE   
   
   real(SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_SHR_KIND_R8  ! pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_CDAY   = 86400.0_SHR_KIND_R8  ! sec in Earth calendar day 
   real(SHR_KIND_R8),parameter :: SHR_CONST_SDAY   = 86164.0_SHR_KIND_R8      ! sec in siderial day ~ sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/5.799e4_SHR_KIND_R8 ! uranus rot ~ rad/sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 2.4764e7_SHR_KIND_R8    ! ~ radius of uranus ~ m
   real(SHR_KIND_R8),parameter :: SHR_CONST_G      = 11.140_SHR_KIND_R8      ! acceleration of gravity ~ m/s^2
   real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 200000.0_SHR_KIND_R8     ! standard pressure ~ pascals
   real(SHR_KIND_R8),parameter :: SHR_CONST_QH2    = 0.810_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QHE    = 0.190_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: SHR_CONST_QN2    = 1.0e-10_SHR_KIND_R8

   real(SHR_KIND_R8),parameter :: SHR_CONST_STEBOL = 5.67e-8_SHR_KIND_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(SHR_KIND_R8),parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_SHR_KIND_R8  ! Boltzmann's constant ~ J/K/molecule
   real(SHR_KIND_R8),parameter :: SHR_CONST_AVOGAD = 6.02214e26_SHR_KIND_R8   ! Avogadro's number ~ molecules/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWDAIR = 2.0*SHR_CONST_QH2+4.0*SHR_CONST_QHE       ! molecular weight dry air ~ kg/kmole, includes helium
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
! The following expression for "frozen" specific heat approximates CP_H2 as 3.0*R, a compromise to true temperature/fp variation
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = (SHR_CONST_QH2*3.0+SHR_CONST_QHE*1.5)*SHR_CONST_RDAIR
!!!   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = 1.2654e4_SHR_KIND_R8     ! specific heat of normal hydrogen ~ J/kg/K
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

   ! FAO:  Additional parameters for Neptune
   real(SHR_KIND_R8),parameter :: PLANET_DAY_RATIO = 0.6712_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_YEAR_RATIO = 164.7935_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_ECCENTRICITY = 0.009_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_OBLIQUITY = 0.5166_SHR_KIND_R8  ! = 29.60 degrees (IAU)
   real(SHR_KIND_R8),parameter :: PLANET_MVELP = 46.981_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: PLANET_EMISSIVITY = 0.86_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: shr_const_au = 30.058_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: sun_size=shr_const_pi*(rad_sun/shr_const_au)*(rad_sun/shr_const_au) !angular size of sun at planet       
   real(SHR_KIND_R8),parameter :: shr_const_intrnlflux = 0.433  ! W m^-2
   
#endif


END MODULE shr_const_mod
