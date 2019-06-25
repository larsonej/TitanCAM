!  @(#) aerad.h  McKie  Oct-1995
!  This is the include file for the interface between the aerosol
!  microphysics and radiative transfer components of CARMA.
!
!  Global symbolic constants are defined and common blocks are
!  declared.
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Define array dimensions.
!
!
!  Define # grid pts in z direction
!
!c    parameter( NZ_AERAD = 60 ) ! haze model
      parameter( NZ_AERAD = 61 ) ! cloud model  EJL-setting this to CAM levels
!
!
!  Define # vertical grid boundaries
!
      parameter( NZP1_AERAD = NZ_AERAD + 1 )
!
!
!  Define # layers below model domain
!
      parameter( NZ_BELOW = 0 )
!
!
!  Define # layers in radiation model
!
      parameter( NZ_RAD = NZ_AERAD + NZ_BELOW )
!
!
!  Define # particle radius bins
!
!c    parameter( NBIN_AERAD = 35) ! haze model 
      parameter( NBIN_AERAD = 48) ! cloud model
!c    parameter( NBIN_AERAD = 57) ! cloud model cm size
!
!
!  Define # particle groups
!
      parameter( NGROUP_AERAD = 2 )
!
!
!  Define # solar wavelengths
!
      parameter( NSOL_AERAD = 26 )
!
!
!  Define # infrared wavelengths
!
      parameter( NIR_AERAD = 18 )
!
!
!  Define total # wavelengths
!
      parameter( NWAVE_AERAD = NSOL_AERAD + NIR_AERAD )
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Define logical unit numbers.
!
!
!  Define logical unit number for output print file
!
      parameter( LUNOPRT_AERAD = 25 )
!
!
!  Define logical unit number for time step info output
!
      parameter( LUNOSTEP_AERAD = 10 )
!
!
!  Define logical unit number for input restart file
!
      parameter( LUNIRES_AERAD = 11 )
!
!
!  Define logical unit number for output restart file
!
      parameter( LUNORES_AERAD = 12 )
!
!
!  Define logical unit number for output history file
!
      parameter( LUNOHIS_AERAD = 13 )
!
!
!  Define logical unit number for input and output of Mie coefficients
!
      parameter( LUNMIE_AERAD = 14 )
!
!
!  Define logical unit number for print output from radiation submodel
!
      parameter( LUNORAD_AERAD = 15 )
!
!
!  Define logical unit number for input temperature values
!
      parameter( LUNITRHO1 = 17 )
      parameter( LUNITRHO2 = 18 )
!
!
!  Define logical unit number for aerosol fluxes input file
!
      parameter( LUNIFLX = 19 )
!
!
!  Define logical unit number for particle concentration input file
!
      parameter( LUNIAER1 = 23 )
      parameter( LUNIAER2 = 24 )
!
!
!  Define logical unit number for timestep output file
!
      parameter( LUNOTME = 21 )
!
!
!  Define logical unit number for growth rate output file
!
      parameter( LUNOGRT = 22 )
!
!
!  Define logical unit number for extrapolation input file
!
      parameter( LUNIEXT = 26 )
!
!
!  Define logical unit number for mass flux output file
!
      parameter( LUNOMFLX = 27 )
!
!
!  Define logical unit number for temporary output file
!
      parameter( LUNOTEMP = 28 )
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
!
      logical is_grp_ice_aerad, do_below

      common / aerad1 / &
        r_aerad(NBIN_AERAD,NGROUP_AERAD), &
        rup_aerad(NBIN_AERAD,NGROUP_AERAD),&
        p_aerad(NZ_RAD), t_aerad(NZ_RAD), &
        pc_aerad(NZ_RAD,NBIN_AERAD,NGROUP_AERAD),&
        qv_aerad(NZ_RAD),&
        ptop_aerad, pbot_aerad, u0_aerad, sfc_alb_aerad, &
        emisir_aerad, tsfc_aerad, h2ocol_aerad,&
        is_grp_ice_aerad(NGROUP_AERAD), do_below,&
        isl_aerad, ir_aerad,&
        iaerad1
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare common blocks for output from radiative transfer code
!
!   heati_aerad    infrared heating rates [K/s]
!   heats_aerad    solar heating rates [K/s]
!   qrad_aerad     particle radiative heating rates [K/s]
!   alb_tomi       spectrally-integrated albedo at top-of-model
!   alb_toai       spectrally-integrated albedo at top-of-atmosphere
!   alb_toa        spectrally-resolved albedo at top-of-atmosphere
!   opd_aerad      spectrally-resolved optical depth
!   fsl_up_aerad   solar upwelling flux [W m^-2]
!   fsl_dn_aerad   solar downwelling flux [W m^-2]
!   fir_up_aerad   infrared upwelling flux [W m^-2]
!   fir_dn_aerad   infrared downwelling flux [W m^-2]
!
      common / aerad2 / &
        wave_aerad(NWAVE_AERAD+1),&
        heati_aerad(NZ_RAD), heats_aerad(NZ_RAD),&
        qrad_aerad(NBIN_AERAD,NZ_RAD,NGROUP_AERAD),&
        alb_tomi_aerad, alb_toai_aerad,&
        alb_toa_aerad(NSOL_AERAD), opd_aerad(NWAVE_AERAD),&
        fsl_up_aerad(NZ_RAD+1), fsl_dn_aerad(NZ_RAD+1),&
        fir_up_aerad(NZ_RAD+1), fir_dn_aerad(NZ_RAD+1),&
        iaerad2
