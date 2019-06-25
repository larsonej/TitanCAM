c  @(#) aerad.h  McKie  Oct-1995
c  This is the include file for the interface between the aerosol
c  microphysics and radiative transfer components of CARMA.
c
c  Global symbolic constants are defined and common blocks are
c  declared.
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c  Define array dimensions.
c
c
c  Define # grid pts in z direction
c
      parameter( NZ_AERAD = 50 )
c
c
c  Define # vertical grid boundaries
c
      parameter( NZP1_AERAD = NZ_AERAD + 1 )
c
c
c  Define # layers below model domain
c
      parameter( NZ_BELOW = 0 )
c
c
c  Define # layers in radiation model
c
      parameter( NZ_RAD = NZ_AERAD + NZ_BELOW )
c
c
c  Define # particle radius bins
c
      parameter( NBIN_AERAD = 41 )
c
c
c  Define # particle groups
c
      parameter( NGROUP_AERAD = 2 )
c
c
c  Define # solar wavelengths
c
      parameter( NSOL_AERAD = 26 )
c
c
c  Define # infrared wavelengths
c
      parameter( NIR_AERAD = 18 )
c
c
c  Define total # wavelengths
c
      parameter( NWAVE_AERAD = NSOL_AERAD + NIR_AERAD )
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c  Define logical unit numbers.
c
c
c  Define logical unit number for output print file
c
      parameter( LUNOPRT_AERAD = 25 )
c
c
c  Define logical unit number for time step info output
c
      parameter( LUNOSTEP_AERAD = 10 )
c
c
c  Define logical unit number for input restart file
c
      parameter( LUNIRES_AERAD = 11 )
c
c
c  Define logical unit number for output restart file
c
      parameter( LUNORES_AERAD = 12 )
c
c
c  Define logical unit number for output history file
c
      parameter( LUNOHIS_AERAD = 13 )
c
c
c  Define logical unit number for input and output of Mie coefficients
c
      parameter( LUNMIE_AERAD = 14 )
c
c
c  Define logical unit number for print output from radiation submodel
c
      parameter( LUNORAD_AERAD = 15 )
c
c
c  Define logical unit number for input temperature values
c
      parameter( LUNITAEM = 17 )
c
c
c  Define logical unit number for temporary output file
c
      parameter( LUNOTEMP = 18 )
c
c
c  Define logical unit number for aerosol fluxes input file
c
      parameter( LUNIFLX = 19 )
c
c
c  Define logical unit number for particle concentration input file
c
      parameter( LUNIAPC = 23 )
c
c
c  Define logical unit number for timestep output file
c
      parameter( LUNOTME = 21 )
c
c
c  Define logical unit number for growth rate output file
c
      parameter( LUNOGRT = 22 )
c
c
c  Define logical unit number for extrapolation input file
c
      parameter( LUNIEXT = 26 )
c
c
c  Define logical unit number for mass flux output file
c
      parameter( LUNOMFLX = 27 )
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c  Declare common blocks for input to radiative transfer code
c
c   is_grp_ice     =.true. means group is ice crystals
c   r_aerad        radius mid-pts from aerosol grids [cm]
c   rup_aerad      upper radii from aerosol grids [cm]
c   p_aerad        pressure [dyne/cm^2]
c   t_aerad        temperature [K]
c   pc_aerad       particle concentration [#/cm^3]
c   qv_aerad       water vapor mixing ratio [g/g]
c   ptop_aerad     pressure at top of aerosol model domain [dyne/cm^2]
c   pbot_aerad     pressure at bottom of aerosol model domain [dyne/cm^2]
c   u0_aerad       cosine of solar zenith angle [dimensionless]
c   sfc_alb_aerad  surface albedo [dimensionless]
c   emisir_aerad   surface IR emissivity [dimensionless]
c   tsfc_aerad     surface temperature [K]
c   h2ocol_aerad   water vapor column above model domain [g/cm^2]
c   isl_aerad      =1 means do solar calculations
c   ir_aerad       =1 means do infrared calculations
c   do_below       =.true. means include radiative layers below model domain
c
      logical is_grp_ice_aerad, do_below

      common / aerad1 /
     $  r_aerad(NBIN_AERAD,NGROUP_AERAD), 
     $  rup_aerad(NBIN_AERAD,NGROUP_AERAD),
     $  p_aerad(NZ_RAD), t_aerad(NZ_RAD), 
     $  pc_aerad(NZ_RAD,NBIN_AERAD,NGROUP_AERAD),
     $  qv_aerad(NZ_RAD),
     $  ptop_aerad, pbot_aerad, u0_aerad, sfc_alb_aerad, 
     $  emisir_aerad, tsfc_aerad, h2ocol_aerad,
     $  is_grp_ice_aerad(NGROUP_AERAD), do_below,
     $  isl_aerad, ir_aerad,
     $  iaerad1
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c  Declare common blocks for output from radiative transfer code
c
c   heati_aerad    infrared heating rates [K/s]
c   heats_aerad    solar heating rates [K/s]
c   qrad_aerad     particle radiative heating rates [K/s]
c   alb_tomi       spectrally-integrated albedo at top-of-model
c   alb_toai       spectrally-integrated albedo at top-of-atmosphere
c   alb_toa        spectrally-resolved albedo at top-of-atmosphere
c   opd_aerad      spectrally-resolved optical depth
c   fsl_up_aerad   solar upwelling flux [W m^-2]
c   fsl_dn_aerad   solar downwelling flux [W m^-2]
c   fir_up_aerad   infrared upwelling flux [W m^-2]
c   fir_dn_aerad   infrared downwelling flux [W m^-2]
c
      common / aerad2 /
     $  wave_aerad(NWAVE_AERAD+1),
     $  heati_aerad(NZ_RAD), heats_aerad(NZ_RAD),
     $  qrad_aerad(NBIN_AERAD,NZ_RAD,NGROUP_AERAD),
     $  alb_tomi_aerad, alb_toai_aerad,
     $  alb_toa_aerad(NSOL_AERAD), opd_aerad(NWAVE_AERAD),
     $  fsl_up_aerad(NZ_RAD+1), fsl_dn_aerad(NZ_RAD+1),
     $  fir_up_aerad(NZ_RAD+1), fir_dn_aerad(NZ_RAD+1),
     $  iaerad2