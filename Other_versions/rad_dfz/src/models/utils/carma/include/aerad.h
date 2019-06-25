!  @(#) aerad.h  McKie  Oct-1995
!  This is the include file for the interface between the aerosol
!  microphysics and radiative transfer components of CARMA.
!
!  Global symbolic constants are defined and common blocks are
!  declared.
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Start of user-defined symbolic constants 
!
!
!  Define # grid pts in x, y, z directions
!
      parameter( NX =  1 )
      parameter( NY =  1 )
!      parameter( NZ =  1 )
      parameter( NZ =  61 )
!
!
!  Define maximum of NX or NY
!
      parameter( NXORNY = NX )
!
!
!  Define # x, y direction grid box boundaries
!
      parameter( NXP1 = NX + 1 )
      parameter( NYP1 = NY + 1 )
!
!
!  Define maximum of NXP1 or NYP1
!
      parameter( NXORNYP1 = NXP1 )
!
!
!  Define # particle radius bins
!
      parameter( NBIN = 30) !EJL
!
!
!  Define # particle elements 
!
!     parameter( NELEM = 5 )
      parameter( NELEM = 1 )
!
!
!  Define # particle groups
!
!     parameter( NGROUP = 3 )
      parameter( NGROUP = 1)
!
!  Define optical groups EJL 5-13-13
!
      parameter( NOPT = 3 )
!
!  Define # solutes
!
      parameter( NSOLUTE = 1 )
!
!
!  Define # gases
!
      parameter( NGAS = 1 )
!
!
!  Define # solar wavelength bins
!
      parameter( NSOL = 26 )
!
!
!  Define # infrared wavelength bins
!
      parameter( NIR = 18 )
!
!
!  Define total # wavelength bins
!
      parameter( NWAVE = NSOL + NIR )
!
!
!  Define # layers in rad xfer model domain underlying aerosol model domain
!
      parameter( NZ_BELOW = 0 )
!
!
!  End of user-defined symbolic constants
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  The remaining symbolic constants will need no attention from most
!  users
!
!
!  Define # layers in radiation model
!
      parameter( NZ_RAD = NZ + NZ_BELOW )
!
!
!  Define # vertical grid boundaries
!
      parameter( NZP1 = NZ + 1 )
!
!
!  Define logical unit number for output print file
!
!      parameter( LUNOPRT = 10 )
      parameter( LUNOPRT = 6 )
!
!
!  Define logical unit number for time step info output
!
!      parameter( LUNOSTEP = 11 )
      parameter( LUNOSTEP = 6 )
!
!
!  Define logical unit number for input restart file
!
      parameter( LUNIRES = 12 )
!
!
!  Define logical unit number for output restart file
!
      parameter( LUNORES = 13 )
!
!
!  Define logical unit number for output history file
!
      parameter( LUNOHIS = 14 )
!
!
!  Define logical unit number for input and output of Mie coefficients
!
      parameter( LUNMIE = 15 )
!
!
!  Define logical unit number for print output from radiation submodel
!
      parameter( LUNORAD = 16 )
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare common blocks for input to radiative transfer code
!
!   is_grp_ice     =.true. means group is ice crystals
!   r_aerad        radius mid-pts from aerosol grids [cm]
!   rup_aerad      upper radii from aerosol grids [cm]
!   p_aerad        pressure [dyne/cm^2]
!   t_aerad        temperature [K]
!   pc_aerad       particle concentration [#/cm^3]
!   qv_aerad       water vapor mixing ratio [g/g]
!   tabove_aerad   blackbody temperature for downwelling IR flux into model [K]
!   ptop_aerad     pressure at top of aerosol model domain [dyne/cm^2]
!   pbot_aerad     pressure at bottom of aerosol model domain [dyne/cm^2]
!   u0_aerad       cosine of solar zenith angle [dimensionless]
!   sfc_alb_aerad  surface albedo [dimensionless]
!   emisir_aerad   surface IR emissivity [dimensionless]
!   tsfc_aerad     surface temperature [K]
!   h2ocol_aerad   water vapor column above model domain [g/cm^2]
!   isl_aerad      =1 means do solar calculations
!   ir_aerad       =1 means do infrared calculations
!   do_below       =.true. means include radiative layers below model domain
!   ir_above_aerad =1 means include downwelling flux into top of model domain
!
      logical is_grp_ice_aerad, do_below

      common / aerad1 /                                                 & 
     &  is_grp_ice_aerad(NGROUP), do_below,                             &
     &  isl_aerad, ir_aerad, ir_above_aerad,                            &
     &  r_aerad(NBIN,NGROUP),                                           &
     &  rup_aerad(NBIN,NGROUP),                                         &
     &  p_aerad(NZ_RAD), t_aerad(NZ_RAD),                               &
     &  pc_aerad(NZ_RAD,NBIN,NGROUP),                                   &
     &  qv_aerad(NZ_RAD), tabove_aerad,                                 &
     &  ptop_aerad, pbot_aerad, u0_aerad, sfc_alb_aerad,                &
     &  emisir_aerad, tsfc_aerad, h2ocol_aerad,                         &
     &  iaerad1
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare common blocks for output from radiative transfer code
!
!   heati_aerad    infrared heating rates [K/s]
!   heats_aerad    solar heating rates [K/s]
!   qrad_aerad     particle radiative heating rates [K/s]
!   alb_tomi_aerad spectrally-integrated albedo at top-of-model
!   alb_toai_aerad spectrally-integrated albedo at top-of-atmosphere
!   alb_toa_aerad  spectrally-resolved albedo at top-of-atmosphere
!   opd_aerad      spectrally-resolved optical depth
!   fsl_up_aerad   solar upwelling flux [W m^-2]
!   fsl_dn_aerad   solar downwelling flux [W m^-2]
!   fir_up_aerad   infrared upwelling flux [W m^-2]
!   fir_dn_aerad   infrared downwelling flux [W m^-2]
!   do_mie_aerad   =.true. means calculate and write Mie coefficients
!
      logical do_mie_aerad

      common / aerad2 /                                                 &
     &  wave_aerad(NWAVE+1),                                            &
     &  heati_aerad(NZ_RAD), heats_aerad(NZ_RAD),                       &
     &  qrad_aerad(NBIN,NZ_RAD,NGROUP),                                 &
     &  alb_tomi_aerad, alb_toai_aerad,                                 &
     &  alb_toa_aerad(NSOL), opd_aerad(NWAVE),                          &
     &  fsl_up_aerad(NZ_RAD+1), fsl_dn_aerad(NZ_RAD+1),                 &
     &  fir_up_aerad(NZ_RAD+1), fir_dn_aerad(NZ_RAD+1),                 &
     &  do_mie_aerad, iblackbody_above,                                 &
     &  iaerad2
