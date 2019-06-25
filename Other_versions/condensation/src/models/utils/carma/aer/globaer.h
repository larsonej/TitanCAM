!  @(#) globaer.h  McKie  Oct-1995
!  This is the global include file for the Ames Aerosol model.
!  This file is intended to be included in all major model
!  source code modules that need access to global variables.
!  All constants and variables in this file are assumed to be
!  available for use in all major source code modules.
!
!  Note:  Some support source code modules that do not need any
!         model variables might not include this file.
!
!  This file contains the following:
!
!   Symbolic constant declarations & definitions.
!   Common block declarations.
!
!  Note:  When the brief definitions given here are not
!         entirely clear, pointers to the source code routines
!         with more complete definitions are given in brackets: {}
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare and define symbolic constants
!   (Implicit typing in effect unless explicitly specified)
!
!
!  Include implicit declarations
!
      include 'precision.h'
!
!
!  Include symbolic constants and common blocks shared between aerosol and
!  radiation models
!
      include 'aerad.h'
!
!
!  Define text string name of this model
!
      character*(*) PROGNAM
      parameter( PROGNAM = 'CARMA' )
!
!
!  Define version tag string for this version of the model
!
      character*(*) PROGTAG
      parameter( PROGTAG = '2.3' )
!
!
!  Define flag to indicate if debugging mode is active
!
      logical DEBUG
      parameter( DEBUG = .false. )
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!  Start of user-defined symbolic constants 
!
!
      parameter( NELPGS = NELEM + NGAS) !EJL 1-31-13
!
!  Define particle number concentration [ # / cm^3 ]
!  used to decide whether to bypass microphysical processes:
!  set it to SMALL_PC to never bypass the calculations.
!
!     parameter( FEW_PC = SMALL_PC*1e5 )
      parameter( FEW_PC = 0. )
!
!
!  Define core fraction (for core mass and second moment) used
!  when particle number concentrations are limited to SMALL_PC
!
      parameter( FIX_COREF = ONE * 0.001 )
!
!
!  End of user-defined symbolic constants
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!  The remaining symbolic constants will need no attention from most
!  users (unless extending the model capabilities or simulating other than
!  a terrestrial atmosphere)
!  
!
!  Define # vertical grid box boundaries, including top & bottom
!
      parameter( NZRADP1 = NZ_RAD + 1 )
!
!
!  Define # components of variables with 3 spatial dimensions
!   (used for collapsing first few dimensions in some calcs)
!
      parameter( NXY = NX * NY )
      parameter( NXYZ = NXY * NZ )
      parameter( NXYZP1 = NXY * NZP1 )
      parameter( NXYZRAD = NXY * NZ_RAD )
      parameter( NXYZRADP1 = NXY * NZRADP1 )
!
!
!  Define # components of pc() in first 4 and 5 dimensions
!   (used for collapsing dimensions of <pc> in some calcs)
!
      parameter( NPC4 = NXYZ * NBIN )
      parameter( NPC5 = NPC4 * NELEM )
!
!
!  Define integer safety marker value placed at end of common blocks
!
      parameter( ISAFETY = 12345 )
!
!
!  Define character safety marker value placed at end of common blocks
!
      character*(5) CSAFETY
      parameter( CSAFETY = '12345' )
!
!
!  Define values of flag used for specification of
!  horizontal transport algorithm
!
      parameter( I_PPM = 0 )
      parameter( I_GALERKIN = 1 )
!
!
!  Define values of flag used for vertical transport
!  boundary conditions
!
      parameter( I_FIXED_CONC = 0 )
      parameter( I_FLUX_SPEC = 1 )
!
!
!  Define values of flag used for particle element
!  composition specification
!
      parameter( I_H2SO4 = 0 )
      parameter( I_WATER = 1 )
      parameter( I_ICE = 2 )
      parameter( I_MIXEDWAT = 3 )
      parameter( I_DUST = 4 )
      parameter( I_C2H6_ICE = 5 )
      parameter( I_C2H6_LIQ = 6 )
      parameter( I_CH4_ICE = 7 )
      parameter( I_CH4_LIQ = 8 )
!
!
!  Define values of flag used for particle element
!  type specification
!
      parameter( I_INVOLATILE = 0 )
      parameter( I_VOLATILE = 1 )
      parameter( I_COREMASS = 2 )
      parameter( I_VOLCORE = 3 )
      parameter( I_CORE2MOM = 4 )
      parameter( I_GROWCORE = 5 )
!
!
!  Define values of flag used for nucleation process
!  specification
!
      parameter( I_DROPACT = 0 )
      parameter( I_AERFREEZE = 1 )
      parameter( I_DROPFREEZE = 2 )
      parameter( I_MIXEDFREEZE = 3 )
      parameter( I_MIXEDMELT = 4 )
      parameter( I_ICEMELT = 5 )
      parameter( I_VAPORDEP = 6) !EJL 1-31-13 this differs from Erika
!
!
!  Define values of flag used specify direction in
!  horizontal transport calculations
!
      parameter( IDIRX = 0 )
      parameter( IDIRY = 1 )
!
!
!  Define values of symbols used to specify horizontal & vertical grid type.
!   Grid selection is made by defining each of the variables
!   <igridv> and <igridh> to one of the grid types known to the model.
!
!   Possible values for igridv:
!       I_CART    cartesian
!       I_SIG     sigma
!
!    Possible values for igridh:
!       I_CART   cartesian
!       I_LL     longitude_latitude
!       I_LC     lambert_conformal
!       I_PS     polar_stereographic
!       I_ME     mercator
!
      parameter( I_CART = 1 )
      parameter( I_SIG = 2 )
      parameter( I_LL = 3 )
      parameter( I_LC = 4 )
      parameter( I_PS = 5 )
      parameter( I_ME = 6 )
      parameter( I_HYB = 7 )
!
!
!  Define values of flag used to specify calculation of solar zenith angle
!
      parameter( I_FIXED = 0 )
      parameter( I_DIURNAL = 1 )
!
!
!  Define triple-point temperature (K)
!
!      parameter( T0 = 90.348d+0 ) !EJL - now function of igas in initgas
!
!
!  Define circle constant [ unitless ]
!
      parameter( PI = 3.14159265358979d+0 )
!
!
!  Define degrees to radian & radian to degrees factors [unitless]
!
      parameter( DEG2RAD = PI / 180. )
      parameter( RAD2DEG = 180. / PI )
!
!
!  Define acceleration of gravity near Titan surface [ cm/s^2 ]
!
      parameter( GRAV = 135.0d+0 )
!
!
!  Define planet equatorial radius [ cm ]
!
      parameter( REARTH = 2.575d+8 )   !EJL - Titan
!
!
!  Define avogadro's number [ # particles / mole ]
!
      parameter( AVG = 6.02252d+23 )
!      parameter( AVOGAD = 6.02252d+23 )
!
!
!  Define Boltzmann's constant [ erg / deg_K ]
!
      parameter( BK = 1.38054d-16 )
!
!
!  Define Loschmidt's number [ mole / cm^3, @ STP ]
!
      parameter( ALOS = 2.68719d+19 )
!
!
!  Define molecular weight of dry air [ g / mole ]
!
      parameter( WTMOL_AIR = 28.0d+0 )
!
!
!  Define reference pressure, e.g. for potential temp calcs [ dyne / cm^2 ]
!
      parameter( PREF = 1013.d+3 ) ! Should this be changed? EJL
!
!
!  Define conversion factor for mb to cgm [ dyne / cm^2 ] units
!
      parameter( RMB2CGS = 1000.d+0 )
!
!
!  Define universal gas constant [ erg / deg_K / mole ]
!
      parameter( RGAS = 8.31430d+07 )
!
!
!  Define gas constant for dry air [ erg / deg_K / mole ]
!
      parameter( R_AIR = RGAS / WTMOL_AIR )
!
!
!  Define number of seconds per the planet's day [ s / d ]
!
      parameter( SCDAY = 1.378080d+6 )
!
!
!  Define specific heat at constant pres of dry air [ cm^2 / s^2 / deg_K ]
!
      parameter( CP = 1.039d+7 )
!
!
!  Define ratio of gas constant for dry air and specific heat
!
      parameter( RKAPPA = R_AIR / CP )
!
!
!  Define mass density of liquid water [ g / cm^3 ]
!  EJL 1-31-13 defined for methane and ethane
      parameter( RHO_W = 0.5446 ) !ethane
!      parameter(RHO_W = 0.4228) !methane
!
!
!  Define mass density of water ice [ g / cm^3 ]
! EJL 1-31-13 changed to ethane and methane
      parameter( RHO_I = 0.713 ) !ethane
!      parameter( RHO_I = 0.519 ) !methane
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for model startup control
!   (These variables are defined fresh for each new run, both
!    for cold starts and restarts, and do not get dumped into
!    the output restart file at the end of a run.  They are
!    all defined at the beginning of the init routine.  They
!    control the type of initialization that will be performed,
!    as well as control things that can change from run to run
!    within a single simulation.)
!
!   ibtime    Beginning timestep index for this run 
!   ietime    Ending timestep index for this run 
!   endtime   Total simulation time for this run
!   nprint    # of timesteps between print reports (used when > 0)
!   nhist     # of timesteps between history output (used when > 0)
!   nrest     # of timesteps between restart output (used when > 0)
!   pprint    time period between print reports  (used when nprint < 0)
!   phist     time period between history outputs (used when nhist < 0)
!   prest     time period between restart outputs (used when nrest < 0)
!   khist     Counter for # of history timepoints output this run
!   prtofil   Name of output print file
!   resifil   Name of input restart file
!   resofil   Name of output restart file
!   hisofil   Name of output history file
!   stepofil  Name of time-step diagnostics file
!   radofil   Name of radiation submodel print output file
!   do_print  .t. if print output during timestepping is desired
!   do_print_setup  .t. if print output during setup is desired
!   do_hist   .t. if history output during timestepping is desired
!   do_rest   .t. if restart output during timestepping is desired
!
      character*(50) prtofil
      character*(50) resifil
      character*(50) resofil
      character*(50) hisofil
      character*(50) stepofil
      character*(50) radofil
      logical do_print
      logical do_print_setup
      logical do_hist
      logical do_rest
!
      common /aer0/                                                     &
     &  pprint, phist, prest,                                           &
     &  ibtime, ietime, endtime,                                        &
     &  nprint, nhist, nrest,                                           &
     &  do_print, do_hist, do_rest, do_print_setup,                     &
     &  khist 
!
      common /aer0s/                                                    &
     &  prtofil, resifil, resofil, hisofil, stepofil, radofil
!
! EJL 1-31-13 unsure about steps above
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for model grid
!
!   dom_llx    Domain limit, lower left x coordinate                 {initatm}
!   dom_lly    Domain limit, lower left y coordiante                 {initatm}
!   dom_urx    Domain limit, upper right x coordinate                {initatm}
!   dom_ury    Domain limit, upper right y coordinate                {initatm}
!   rlon0      center longitude for LC, PS, LL projections           {initatm}
!   rlat0      true latitude for ME projection                       {initatm}
!   rlat1      #1 true latitude for LC projection                    {initatm}
!   rlat2      #2 true latitude for LC projection                    {initatm}
!   hemisph    +1.=southern, -1.=northern hemisphere for PS, LC proj {initatm}
!   igridv     flag to specify desired vertical grid coord system    {initatm}
!   igridh     flag to specify desired horizontal grid coord system  {initatm}
!   zl         Altitude at top of layer                              {initatm}
!   zlold      Altitude at top of layer at start of time step
!   zc         Altitude at layer mid-point                           {initatm}
!   zcold      Altitude at layer mid-point at start of time step
!   xc         Horizontal position at center of box                  {initatm}
!   yc         Horizontal position at center of box                  {initatm}
!   xl         Horizontal position at lower edge of box              {initatm}
!   yl         Horizontal position at lower edge of box              {initatm}
!   xu         Horizontal position at upper edge of box              {initatm}
!   yu         Horizontal position at upper edge of box              {initatm}
!   dx         Horizontal grid spacing                               {initatm}
!   dy         Horizontal grid spacing                               {initatm}
!   dz         Thickness of vertical layers                          {initatm}
!   xmet       Horizontal ds/dx (ds is metric distance)              {initatm}
!   ymet       Horizontal ds/dy (ds is metric distance)              {initatm}
!   zmet       Vertical ds/dz (ds is metric distance)                {initatm}
!   zmetl      Vertical ds/dz at beginning of time-step              
!   rlon,rlat  Longitude, latitude [deg] at each horiz grid point    {initatm}
!   gridname   Text description of horiz & vert grid coord system    {initatm}
!   iaer1      Safety marker for common block aer1
!
      logical do_huygens_Tprofile
      character*(80) gridname
      character*(5) caer1s
!
      common /aer1/                                                     & 
     &  dom_llx, dom_lly, dom_urx, dom_ury,                             &
     &  rlon0, rlat0, rlat1, rlat2, hemisph,                            &
     &  zl(NX,NY,NZP1), zc(NX,NY,NZ),                                   &
     &  zlold(NXYZP1), zcold(NXYZ),                                     &
     &  xc(NX,NY,NZ), xl(NX,NY,NZ), xu(NX,NY,NZ),                       &
     &  yc(NX,NY,NZ), yl(NX,NY,NZ), yu(NX,NY,NZ),                       &
     &  dx(NX,NY,NZ), dy(NX,NY,NZ), dz(NX,NY,NZ),                       &
     &  xmet(NX,NY,NZ), ymet(NX,NY,NZ), zmet(NX,NY,NZ),                 &
     &  zmetl(NX,NY,NZ),                                                &
     &  rlon(NX,NY), rlat(NX,NY),                                       &
     &  igridv, igridh, do_huygens_Tprofile,                            &
     &  iaer1
!
      common /aer1s/                                                    &
     &  gridname,                                                       &
     &  caer1s
!
!
!  Declare alias names for grid stuff with first 2, 3 dimensions treated linearly
!
      dimension zl2(NXY,NZP1), zl3(NXYZP1)
      dimension zc2(NXY,NZ), zc3(NXYZ)
      dimension dz2(NXY,NZ), dz3(NXYZ)
      dimension xc2(NXY,NZ), xc3(NXYZ)
      dimension yc2(NXY,NZ), yc3(NXYZ)
      dimension xl2(NXY,NZ), xl3(NXYZ)
      dimension yl2(NXY,NZ), yl3(NXYZ)
      dimension xu2(NXY,NZ), xu3(NXYZ)
      dimension yu2(NXY,NZ), yu3(NXYZ)
      dimension dx2(NXY,NZ), dx3(NXYZ)
      dimension dy2(NXY,NZ), dy3(NXYZ)
      dimension xmet2(NXY,NZ), xmet3(NXYZ)
      dimension ymet2(NXY,NZ), ymet3(NXYZ)
      dimension zmet2(NXY,NZ), zmet3(NXYZ)
      dimension zmetl2(NXY,NZ), zmetl3(NXYZ)
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for model option & control variables
!
!   time        Simulation time at end of current timestep [s]
!   dtime       Timestep size [s]
!   dtmin       Minimum time-step
!   dtmax       Maximum time-step
!   dpctol      Maximum change in particle concentrations that will be tolerated
!   dgstol      Maximum change in gas concentrations that will be tolerated
!   conmax      Minumum relative concentration to consider in varstep   {prestep}
!   itime       Timestep index at end of current timestep 
!   igelem      Groups to which elements belong                     {setupbins}
!   itype       Particle type specification array                   {setupbins}
!   icomp       Particle compound specification array               {setupbins}
!   nelemg      Number of elements in group           
!   ncore       Number of core elements (itype = 2) in group           
!   icorelem    Core elements (itype = 2) in group           
!   ienconc     Particle number conc. element for group             {setupbins}
!   ishape      Describes particle shape for group
!   icoag       Coagulation mapping array                           {setupcoag}
!   icoagelem   Coagulation element mapping array                   {setupcoag}
!   icoagelem_cm Coagulation element mapping array for second mom   {setupcoag}
!   ifall       Fall velocity options                               {setupvfall}
!   icoagop     Coagulation kernel options                          {setupckern}
!   icollec     Gravitational collection options                      {setupckern}
!   itbnd       Top boundary condition                                {vertical}
!   ibbnd       Bottom boundary condition                             {vertical}
!   itbnd_pc    Top boundary condition flag for particles             {init}
!   ibbnd_pc    Bottom boundary condition flag for particles          {init}
!   itbnd_gc    Top boundary condition flag for gas                   {init}
!   ibbnd_gc    Bottom boundary condition flag for gas                {init}
!   itbnd_ptc   Top boundary condition flag for potential temp.       {init}
!   ibbnd_ptc   Bottom boundary condition flag for potential temp.    {init}
!   ihoradv     Specification of horizontal advection algorithm       {init}
!   do_coag     If .true. then do coagulation                         {init}
!   do_grow     If .true. then do condensational growth and evap.     {init}
!   do_thermo   If .true. then do solve thermodynamics equation       {init}
!   do_vtran    If .true. then do vertical transport                  {init}
!   do_ew       If .true. then do east-west transport                 {init}
!   do_ns       If .true. then do north-south transport               {init}
!   do_ccoef    If .true. then calculate coefficients for htran       {htranglk}
!   do_varstep  If .true then use variable time-step                  {init}
!   do_step     If .false. then step was too big, so dont advance    {varstep}
!   do_error    If .true. then do error trapping for debugging        {init}
!   do_netcdf   If .true. then output history in netcdf file format   {init}
!   do_parcel   If .true. then do parcel simulation                   {init}
!   ncdf_file   Netcdf file handle integer for internal netcdf use    {outhis_ncdf}
!   sec_mom     If .true. then core second moment (itype = 3) used    {setupgrow}
!   igrowgas    Gas that condenses into a particle element            {setupgrow}
!   inucgas     Gas that nucleates a particle group                   {setupnuc}
!   if_nuc      Nucleation conditional array                          {setupaer}
!   inucproc    Nucleation conditional array                          {setupaer}
!   nnuc2elem   Number of elements that nucleate to element           {setupnuc}
!   inuc2elem   Nucleation transfers particles into element inuc2elem {setupnuc}
!   ievp2elem   Total evap. transfers particles into group ievp2elem  {setupnuc}
!   ievp2bin    Total evap. transfers particles into bin ievp2bin     {setupnuc}
!   inuc2bin    Nucleation transfers particles into bin inuc2bin      {setupnuc}
!   ipownuc     Mass exponent to match source/target elems for nuc    {setupnuc}
!   isolelem    Index of solute for each particle element             {setupnuc}
!   ix          Current index for spatial grid, general east-west direction
!   iy          Current index for spatial grid, general north-south direction
!   iz          Current index for spatial grid, vertical direction
!   ixy         Current index for spatial grid, linearized 2-D ix,iy (horizontal)
!   ixyz        Current index for spatial grid, linearized 3-D ix,iy,iz 
!   ntsubsteps  Number of time substeps for fast microphysics processes
!   maxsubsteps Maximum number of time substeps allowed
!   minsubsteps Maximum number of time substeps allowed
!   ibinver     Bin number in vertical subroutine
!   ielemver    Element number in vertical subroutine
!   iaer2       Safety marker for common block aer2
!
!
!   simtitle   Model simulation title string
!   elemname   Names of particle elements
!   groupname  Names of particle groups
!   gasname    Names of gas species
!   solname    Names of solutes
!   caer2s     Safety marker for common block aer2s
!   elemsname  Short names of particle elements
!   gassname   Short names of gas species
!    
      logical do_coag, do_grow, do_thermo, do_ew, do_ns, is_grp_ice
      logical is_grp_mixed, do_vtran, do_varstep, do_step, do_ccoef
      logical do_error, if_sec_mom, if_nuc(NELEM,NELEM)
      logical if_nuc_lh(NELEM,NELEM)
      logical do_netcdf, do_parcel
      logical is_grp_mixed_comp, is_grp_mixed_phase
      character*(80) simtitle
      character*(50) elemname(NELEM), groupname(NGROUP)
      character*(20) gasname(NGAS), solname(NSOLUTE)
      character*(5) caer2s
      character*(6) elemsname(NELEM)
      character*(8) gassname(NGAS)
!
      common /aer2/                                                     &
     &  time, dtime, dtmin, dtmax, dpctol, dgstol, conmax,              &
     &  time_nuc(NGROUP),                                               &
     &  period_nuc, maxsubsteps, minsubsteps, dtime_save,               &
     &  if_sec_mom(NGROUP),                                             &
     &  itime, ix, iy, iz, ixy, ixyz, ntsubsteps,                       &
     &  igelem(NELEM), itype(NELEM), icomp(NELEM), nelemg(NGROUP),      &
     &  ncore(NGROUP),                                                  &
     &  ishape(NGROUP), ienconc(NGROUP), icoag(NGROUP,NGROUP),          &
     &  icoagelem(NELEM,NGROUP), icoagelem_cm(NELEM,NGROUP),            &
     &  ifall, icoagop, icollec, itbnd,                                 &
     &  ibbnd, itbnd_pc, ibbnd_pc, itbnd_gc, ibbnd_gc, itbnd_ptc,       &
     &  ibbnd_ptc, ihoradv, do_coag, do_grow, do_thermo, do_vtran,      &
     &  do_ew, do_ns, do_varstep, do_step, do_ccoef, do_error,          &
     &  do_netcdf, do_parcel, if_nuc, if_nuc_lh,                        &
     &  ncdf_file,                                                      &
     &  imomelem(NGROUP), inucproc(NELEM,NELEM),                        &
     &  igrowgas(NELEM), inucgas(NGAS,NGROUP), nnuc2elem(NELEM),        &
     &  inuc2elem(NELEM,NELEM), ievp2elem(NELEM),                       &
     &  is_grp_ice(NGROUP), is_grp_mixed(NGROUP),                       &
     &  is_grp_mixed_phase(NGROUP), is_grp_mixed_comp(NGROUP),          &
     &  inuc2bin(NBIN,NGROUP,NGROUP),                                   &
     &  isolelem(NELEM), ievp2bin(NBIN,NGROUP,NGROUP),                  &
     &  icorelem(NELEM,NELEM), nnucelem(NELEM),                         &
     &  nnucbin(NGROUP,NBIN,NGROUP), ipownuc(NELEM,NELEM),              &
     &  inucelem(NELEM,NELEM*NGROUP),                                   &
     &  inucbin(NBIN*NGROUP,NGROUP,NBIN,NGROUP),                        &
     &  ifractal(NGROUP), ibinver, ielemver,                            &
     &  iaer2
!
      common /aer2s/                                                    &
     &  simtitle, elemname, groupname,                                  &
     &  gasname, solname,                                               &
     &  caer2s, elemsname, gassname      
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for particle grid structure
!
!   rmin      Radius of particle in first bin [g] 
!   rmassmin  Mass of particle in first bin [g]
!   rmrat     Ratio of masses of particles in consecutive bins  {setupaer}
!   r         Radius bins [cm]
!   rmass     Mass bins [g]
!   vol       Particle volume [cm^3]
!   dr        Width of bins in radius space [cm]
!   dm        Width of bins in mass space [g]
!   dv        Width of bins in volume space [cm^3]
!   rmassup   Upper bin boundary mass [g]
!   rup       Upper bin boundary radius [cm]
!   rlow      Lower bin boundary radius [cm]
!   diffmass  Difference between <rmass> values
!   small_val Small values used in smallconc()
!   rhop      Mass density of particle groups [g/cm^3]
!   rhoelem   Mass density of particle elements [g/cm^3]
!   eshape    Ratio of particle length / diameter 
!   rf	      Fractal radius
!   df	      fractal dimension
!   umon      number of monomers in aggregate
!   rm        mobility radius
!   rmon      radius of monomers
!   porosity  porosity
!   iaer3     Safety marker for common block aer3
!
      common /aer3/                                                     &
     &  rmin(NGROUP), rmassmin(NGROUP), rmrat(NGROUP),                  &
     &  r(NBIN,NGROUP), rmass(NBIN,NGROUP),                             &
     &  vol(NBIN,NGROUP), dr(NBIN,NGROUP),                              &
     &  dm(NBIN,NGROUP), dv(NBIN,NGROUP),                               &
     &  rmassup(NBIN,NGROUP), rup(NBIN,NGROUP),                         &
     &  rlow(NBIN,NGROUP), small_val(NBIN,NELEM),                       &
     &  diffmass(NBIN,NGROUP,NBIN,NGROUP),                              &
     &  rhop(NX,NY,NZ,NBIN,NGROUP),                                     &
     &  rhoelem(NELEM), eshape(NGROUP),                                 &
     &  rf(NBIN,NGROUP), df(NBIN,NGROUP),                               &
     &  umon(NBIN,NGROUP),                                              &
     &  rm(NBIN,NGROUP),rmon(NGROUP), porosity(NBIN,NGROUP),                    &
     &  iaer3
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for primary model state variables
!
!   pc          Particle concentration                             {initaer}
!   gc          Gas concentration [g/x_units/y_units/z_units]      {initgas}
!   ptc         Potential temperature concentration [g K/x_units/y_units/z_units]
!   iaer4       Safety marker for common block aer4
!
      common /aer4/                                                     &
     &  pc(NX,NY,NZ,NBIN,NELEM),                                        &
     &  gc(NX,NY,NZ,NGAS),                                              &
     &  ptc(NX,NY,NZ),                                                  &
     &  iaer4
!
!
!  Declare alias names for pc() and rhop() with first 2, 3, 4, 5 dimensions
!  treated linearly
!
      dimension pc2(NXY,NZ,NBIN,NELEM)
      dimension pc3(NXYZ,NBIN,NELEM)
      dimension pc4(NPC4,NELEM)
      dimension pc5(NPC5)
      dimension rhop2(NXY,NZ,NBIN,NGROUP)
      dimension rhop3(NXYZ,NBIN,NGROUP)
!
      equivalence( pc2, pc )
      equivalence( pc3, pc )
      equivalence( pc4, pc )
      equivalence( pc5, pc )
      equivalence( rhop2, rhop )
      equivalence( rhop3, rhop )
!
!
!  Declare alias name for gc() and ptc() with first 2, 3, 4, 5 dimensions treated linearly
!
      dimension gc2(NXY,NZ,NGAS)
      dimension gc3(NXYZ,NGAS)
      dimension ptc2(NXY,NZ)
      dimension ptc3(NXYZ)
!
      equivalence( gc2, gc )
      equivalence( gc3, gc )
      equivalence( ptc2, ptc )
      equivalence( ptc3, ptc )
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for secondary model variables
!
!   pcl         Particle concentration at beginning of time-step
!   pcmax       Maximum concentration for each element             {prestep}
!   pconmax     Maximum particle concentration for each grid point
!   gcl         Gas concentration at beginning of time-step
!   ptcl        Potential temperature concentration at beginning of time-step
!   d_pc        Change in particle concentration due to transport
!   d_gc        Change in gas concentration due to transport
!   d_ptc       Change in potential temperature concentration due to transport
!   cvert       Temporary storage for vertical transport 
!   cvert_tbnd  Temporary storage for vertical transport
!   cvert_bbnd  Temporary storage for vertical transport
!   divcor      Correction term for vertical divergence
!   chor        Temporary storage for horizontal transport
!   dhor        Temporary storage of grid spacing for horizontal transport
!   coaglg      Total particle loss rate due to coagulation for group
!   coagpe      Particle production due to coagulation 
!   rnuclg      Total particle loss rate due to nucleation for group
!   rnucpe      Particle production due to nucleation 
!   growlg      Total particle loss rate due to growth for group
!   growle      Partial particle loss rate due to growth for element 
!   growpe      Particle production due to growth 
!   evaplg      Total particle loss rate due to evaporation for group
!   evapls      Partial particle loss rate due to evaporation for element
!   evappe      Particle production due to evaporation
!   coreavg     Average total core mass in bin
!   coresig     logarithm^2 of std dev of core distribution
!   evdrop      Particle production of droplet number
!   evcore      Particle production of core elements
!   gasprod     Gas production term
!   rlheat      Latent heating rate [deg_K/s]
!   vertdifu    Upward vertical flux at grid boundary
!   vertdifd    Downward vertical flux at grid boundary
!   ftopgas     Downward gas flux across top boundary of model
!   fbotgas     Upward gas flux across bottom boundary of model
!   ftoppart    Downward particle flux across top boundary of model
!   fbotpart    Upward flux particle across bottom boundary of model
!   ftop        Downward flux across top boundary of model
!   fbot        Upward flux across bottom boundary of model
!   pc_topbnd   Particle concentration assumed just above the top boundary
!   pc_botbnd   Particle concentration assumed just below the bottom boundary
!   gc_topbnd   Gas concentration assumed just above the top boundary
!   gc_botbnd   Gas concentration assumed just below the bottom boundary
!   ptc_topbnd  Thermodynamic variable value assumed just above the top boundary
!   ptc_botbnd  Thermodynamic variable value assumed just below the bottom boundary
!   cmf         Core mass fraction in a droplet 
!   totevap     .true. if droplets are totally evaporating to CN
!   too_small   .true. if cores are smaller than smallest CN
!   too_big     .true. if cores are larger than largest CN
!   nuc_small   .true. if cores are smaller than smallest nucleated CN
!   inucmin     Index of smallest particle nucleated; used for tot. evap.
!   inucstep    Index of smallest particle nucleated during time step
!   iaer5       Safety marker for common block aer5
!
      logical totevap, too_small, too_big, nuc_small

      common /aer5/                                                     &
     &  pcl(NXYZ,NBIN,NELEM),                                           &
     &  gcl(NXYZ,NGAS), ptcl(NXYZ), d_pc(NXYZ,NBIN,NELEM),              &
     &  d_gc(NXYZ,NGAS), d_ptc(NXYZ),                                   &
     &  pcmax(NELEM), cvert(NZ), divcor(NZ), chor(NXORNY),              &
     &  cvert_tbnd, cvert_bbnd, dhor(NXORNY), pconmax(NXYZ,NGROUP),     &
     &  coaglg(NXYZ,NBIN,NGROUP), coagpe(NXYZ,NBIN,NELEM),              & 
     &  rnuclg(NBIN,NGROUP,NGROUP), rnucpe(NBIN,NELEM),                 &
     &  growlg(NBIN,NGROUP), growpe(NBIN,NELEM),                        &
     &  evaplg(NBIN,NGROUP), evappe(NBIN,NELEM),                        &
     &  gasprod(NGAS), rlheat, evdrop, evcore(NELEM),                   &
     &  vertdifd(NZP1), vertdifu(NZP1),                                 &
     &  ftopgas(NXY,NGAS), fbotgas(NXY,NGAS),                           &
     &  ftoppart(NXY,NBIN,NELEM), fbotpart(NXY,NBIN,NELEM),             &
     &  ftop, fbot, cmf(NBIN,NGROUP), coreavg, coresig,                 &
     &  pc_topbnd(NXY,NBIN,NELEM), pc_botbnd(NXY,NBIN,NELEM),           &
     &  gc_topbnd(NXY,NGAS), gc_botbnd(NXY,NGAS),                       &
     &  ptc_topbnd(NXY), ptc_botbnd(NXY),                               &
     &  totevap(NBIN,NGROUP),                                           &
     &  too_small, too_big, nuc_small,                                  &
     &  inucmin(NGROUP), inucstep(NGROUP),                              &
     &  ppd, rprod(NZ), pls, pcmflux(NX,NY,NZ,NELEM),                   &
     &  dmdte_gro(NBIN,NELEM), dmdt_gro(NBIN,NGROUP),                   &
     &  gprod_grow(NGROUP,NGAS), gprod_drop(NGROUP, NGAS),              &
     &  gprod_mono(NGROUP,NGAS), gprod_poly(NGROUP,NGAS), uc,           &
     &  iaer5
!
!
!  Declare alias name for <evappe> to collapse into one dimension
!
      dimension evappe5(NPC5,NGAS)
      dimension pcmflux2(NXY,NZ,NELEM)
!
      equivalence( evappe5, evappe )
      equivalence( pcmflux2, pcmflux)
!
! EJL 1-31-13 Again, the above may need to change or have things added
! same below
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for coagulation kernels and bin
!  pair mapping
!
!   ck0           Constant coagulation kernel           {setupaer}
!   grav_e_coll0  Constant value for collection effic.  {setupaer}
!   ckernel       Coagulation kernels [/cm^3/s]         {setupckern}
!   pkernel       Coagulation production variables      {setupcoag}
!   volx          Coagulation subdivision variable      {setupcoag}
!   ilow          Bin pairs for coagulation production  {setupcoag}
!   jlow          Bin pairs for coagulation production  {setupcoag}
!   iup           Bin pairs for coagulation production  {setupcoag}
!   jup           Bin pairs for coagulation production  {setupcoag}
!   npairl        Bin pair indices                      {setupcoag}
!   npairu        Bin pair indices                      {setupcoag}
!   cbr_term0     numerator of brownian coag. kern      {setupckern}
!   cbr_term1     from denominator of br. coag. kern    {setupckern}
!   cbr_term2     multiplied by sticking coeff in cbr   {setupckern}
!   pkern0        Coagulation production variables      {setupcoag}
!   prrat         Charge to radius ratio                {charge}
!   iaer6         Safety marker for common block aer6
!
      common /aer6/                                                     &
     &  ck0, grav_e_coll0,                                              &
     &  ckernel(NZ,NBIN,NBIN,NGROUP,NGROUP),                            &
     &  pkernel(NZ,NBIN,NBIN,NGROUP,NGROUP,NGROUP,6),                   &
     &  volx(NGROUP,NGROUP,NGROUP,NBIN,NBIN),                           &
     &  ilow(NGROUP,NBIN,NBIN*NBIN),                                    &
     &  jlow(NGROUP,NBIN,NBIN*NBIN),                                    &
     &  iup(NGROUP,NBIN,NBIN*NBIN),                                     &
     &  jup(NGROUP,NBIN,NBIN*NBIN),                                     &
     &  npairl(NGROUP,NBIN),                                            &
     &  npairu(NGROUP,NBIN),                                            &
     &  cbr_term0(NZ,NBIN,NBIN,NGROUP,NGROUP),                          &
     &  cbr_term1(NZ,NBIN,NBIN,NGROUP,NGROUP),                          &
     &  cbr_term2(NZ,NBIN,NBIN,NGROUP,NGROUP),                          &
     &  pkern0(NZ,NBIN,NBIN,NGROUP,NGROUP,NGROUP,6),                    &
     &  prrat(NZ),                                                      &
     &  iaer6
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for coagulation group pair mapping
!
!   iglow      Group pairs for coagulation production  {setupcoag}
!   iglow      Group pairs for coagulation production  {setupcoag}
!   jglow      Group pairs for coagulation production  {setupcoag}
!   igup       Group pairs for coagulation production  {setupcoag}
!   jgup       Group pairs for coagulation production  {setupcoag}
!   iaer7      Safety marker for common block aer7
!
      common /aer7/                                                     &
     &  iglow(NGROUP,NBIN,NBIN*NBIN),                                   &
     &  jglow(NGROUP,NBIN,NBIN*NBIN),                                   &
     &  igup(NGROUP,NBIN,NBIN*NBIN),                                    &
     &  jgup(NGROUP,NBIN,NBIN*NBIN),                                    &
     &  iaer7
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for particle fall velocities, transport
!  rates, and coagulation kernels
!
!   bpm       Corrects for non-sphericity and non-continuum effects {setupvfall}
!   vf        Fall velocities at layer mid-pt                       {setupvfall}
!   vtrans    Net vertical transport rate at layer boundary         {vertical}
!   re        Reynolds number based on <vfall>                     {setupvfall}
!   vf_const  Constant vertical fall velocity when ifall=0          {setupaer}
!   vertadvu  Upward mass transport rate due to advection	    {vertran}
!   vertadvd  Downward mass transport rate due to advection	    {vertran}
!   htrans    Net horizontal transport rate at layer boundary       {horizont}
!   hdiff     Horizontal diffusion coefficient at layer boundary    {horizont}
!   ca        Coefficient for Galerkin horizontal transport         {glkcoef}}
!   cb        Coefficient for Galerkin horizontal transport         {glkcoef}}
!   cc        Coefficient for Galerkin horizontal transport         {glkcoef}}
!   cd        Coefficient for Galerkin horizontal transport         {glkcoef}}
!   ce        Coefficient for Galerkin horizontal transport         {glkcoef}}
!   cg        Coefficient for Galerkin horizontal transport         {glkcoef}}
!   iaer8     Safety marker for common block aer8
!
      common /aer8/                                                     &
     &  bpm(NZ,NBIN,NGROUP), vf(NZP1,NBIN,NGROUP),                      &
     &  vtrans(NZP1), re(NZ,NBIN,NGROUP), vf_const,                     &
     &  vertadvu(NZP1), vertadvd(NZP1), htrans(NXORNYP1),               &
     &  hdiff(NXORNY),                                                  &
     &  ca(2,NX,NY), cb(2,NX,NY), cd(2,NX,NY),                          &
     &  ce(2,NX,NY), cf(2,NX,NY), cg(2,NX,NY),                          &
     &  iaer8
!
!
!  Declare alias names for <rhostar> with first 2, 3 dimensions treated linearly
!
      dimension rhostar2(NXY,NZ), rhostar3(NXYZ)
!
      equivalence( rhostar2, rhostar )
      equivalence( rhostar3, rhostar )
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for atmospheric structure.
!  Note: air density, winds, and diffusion coefficients are in scaled units
!
!   zbot      height of the bottom of the model [cm]                  {initatm}
!   p_surf    surface pressure [dyne/cm^2]                            {initatm}
!   p_top     Atmospheric pressure at top of model domain [dyne/cm^2] {initatm}
!   p         Atmospheric pressure at layer mid-pt [dyne/cm^2]        {initatm}
!   pl        Atmospheric pressure at layer edge [dyne/cm^2]          {initatm}
!   rhoa      Air density at layer mid-pt [g/x_units/y_units/z_units] {initatm}
!   t         Air temperature at layer mid-pt [deg_K]                 {initatm}
!   t_surf    Air temperature at surface [deg_K]                      {initatm}
!   rmu       Air viscosity at layer mid-pt [g/cm/s]                  {initatm}
!   thcond    Thermal conductivity of dry air [erg/cm/sec/deg_K]      {initatm}
!   w         Vertical wind speed at layer boundary [z_units/s]       {initatm}
!   u         East-west wind speed at layer center [x_units/s]        {initatm}
!   v         North-south wind speed at layer center [y_units/s]      {initatm}
!   dkz       Vert diffusion coef at layer boundary [z_units^2/s]     {initatm}
!   dkx       Horiz x diffusion coeff at layer boundary [x_units^2/s] {initatm}
!   dky       Horiz y diffusion coeff at layer boundary [y_units^2/s] {initatm}
!   told      Temperature at beginning of time-step
!   pold      Pressure at beginning of time-step
!   rhoaold   Air density at beginning of time-step
!   iaer9     Safety marker for common block aer9
!
      common /aer9/                                                     &
     &  zbot, p(NX,NY,NZ), pl(NX,NY,NZP1), rhoa(NX,NY,NZ),              &
     &  p_surf(NX,NY), p_top(NX,NY),                                    &
     &  t(NX,NY,NZ), told(NXYZ), pold(NXYZ), rhoaold(NXYZ),             &
     &  rmu(NZ), thcond(NZ), w(NX,NY,NZP1), zmetold(NXYZ),              &
     &  u(NX,NY,NZ), v(NX,NY,NZ), t_surf(NX,NY),                        &
     &  dkz(NX,NY,NZP1), dkx(NX,NY,NZP1), dky(NX,NY,NZP1),              &
     &  iaer9
!
!
!  Declare alias names for atm stuff with first 2, 3 dimensions treated linearly
!
      dimension p2(NXY,NZ), p3(NXYZ)
      dimension pl2(NXY,NZP1), pl3(NXYZP1)
      dimension t2(NXY,NZ), t3(NXYZ)
      dimension u2(NXY,NZP1), u3(NXYZP1)
      dimension v2(NXY,NZP1), v3(NXYZP1)
      dimension w2(NXY,NZP1), w3(NXYZP1)
      dimension dkx2(NXY,NZP1), dkx3(NXYZP1)
      dimension dky2(NXY,NZP1), dky3(NXYZP1)
      dimension dkz2(NXY,NZP1), dkz3(NXYZP1)
      dimension rhoa2(NXY,NZ), rhoa3(NXYZ)
      dimension p_surf2(NXY)
      dimension p_top2(NXY)
!
      equivalence( zl2, zl )
      equivalence( zl3, zl )
      equivalence( zc2, zc )
      equivalence( zc3, zc )
      equivalence( dz2, dz )
      equivalence( dz3, dz )
      equivalence( xc2, xc )
      equivalence( xc3, xc )
      equivalence( xl2, xl )
      equivalence( xl3, xl )
      equivalence( xu2, xu )
      equivalence( xu3, xu )
      equivalence( yc2, yc )
      equivalence( yc3, yc )
      equivalence( yl2, yl )
      equivalence( yl3, yl )
      equivalence( yu2, yu )
      equivalence( yu3, yu )
      equivalence( dx2, dx )
      equivalence( dx3, dx )
      equivalence( dy2, dy )
      equivalence( dy3, dy )
      equivalence( p2, p )
      equivalence( p3, p )
      equivalence( pl2, pl )
      equivalence( pl3, pl )
      equivalence( rhoa2, rhoa )
      equivalence( rhoa3, rhoa )
      equivalence( t2, t )
      equivalence( t3, t )
      equivalence( u2, u )
      equivalence( u3, u )
      equivalence( v2, v )
      equivalence( v3, v )
      equivalence( w2, w )
      equivalence( w3, w )
      equivalence( dkx2, dkx )
      equivalence( dkx3, dkx )
      equivalence( dky2, dky )
      equivalence( dky3, dky )
      equivalence( dkz2, dkz )
      equivalence( dkz3, dkz )
      equivalence( xmet2, xmet )
      equivalence( xmet3, xmet )
      equivalence( ymet2, ymet )
      equivalence( ymet3, ymet )
      equivalence( zmet2, zmet )
      equivalence( zmet3, zmet )
      equivalence( zmetl2, zmetl )
      equivalence( zmetl3, zmetl )
      equivalence( zmetold2, zmetold )
      equivalence( zmetold3, zmetold )
      equivalence( p_surf2, p_surf )
      equivalence( p_top2, p_top )
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for condensational growth parameters
!
!   gwtmol    Molecular weight for gases [g/mol]                  {setupgrow}
!   diffus    Diffusivity of gas in air [cm^2/s]                  {setupgrow}
!   rlhe      Latent heat of evaporation for gas [cm^2/s^2]       {setupgrow}
!   rlhe      Latent heat of ice melting for gas [cm^2/s^2]       {setupgrow}
!   pvapl     Saturation vapor pressure over water [dyne/cm^2]    {vaporp}   
!   pvapi     Saturation vapor pressure over ice [dyne/cm^2]      {vaporp}   
!   surfctwa  Surface tension of water-air interface              {setupgkern}
!   surfctiw  Surface tension of water-ice interface              {setupgkern}
!   surfctia  Surface tension of ice-air interface                {setupgkern}
!   akelvin   Exponential arg. in curvature term for growth       {setupgkern}
!   akelvini  Curvature term for ice                              {setupgkern}
!   ft        Ventilation factor                                  {setupgkern}
!   gro       Growth kernel [UNITS?]                              {setupgkern}
!   gro1      Growth kernel conduction term [UNITS?]              {setupgkern}
!   gro2      Growth kernel radiation term [UNITS?]               {setupgkern}
!   gvrat     For converting particle growth to gas loss [UNITS?] {setupgkern}
!   pratt     Terms in PPM advection scheme for condensation      {setupgkern}
!   prat
!   pden1
!   palr
!   supsatl   Supersaturation of vapor w.r.t. liquid water [dimless]
!   supsati   Supersaturation of vapor w.r.t. ice [dimless]                  
!   supsatlold Supersaturation (liquid) before time-step    {prestep}
!   supsatiold Supersaturation (ice) before time-step    {prestep}
!   scrit     Critical supersaturation for nucleation [dimless]   {setupnuc}
!   sol_ions  Number of ions solute dissociates into              {setupnuc}
!   solwtmol  Molecular weight of solute                          {setupnuc}
!   T0        Triple point temp for pure substance, EJL 2-12-13   {initgas}
!   Tfreez    Temperature of freezing point                       {microfast}
!   rhosol    Mass density of solute                              {setupnuc}
!   rlh_nuc   Nucleation latent heat                              {setupaer}
!   iaer10   Safety marker for common block aer10
!
      common /aer10/                                                    &
     &  gwtmol(NGAS), diffus(NZ,NGAS),                                  &
     &  rlhe(NZ,NGAS), rlhm(NZ,NGAS),                                   &
     &  pvapl(NX,NY,NZ,NGAS),                                           &
     &  pvapi(NX,NY,NZ,NGAS), vp_ch4_adjust(NZ),                        &
     &  surfctwa(NZ,NGAS), surfctiw(NZ,NGAS), surfctia(NZ,NGAS),        &
     &  akelvin(NZ,NGAS), akelvini(NZ,NGAS),                            &
     &  ft(NZ,NBIN,NGROUP),                                             &
     &  gro(NZ,NBIN,NGROUP),                                            &
     &  gro1(NZ,NBIN,NGROUP),                                           &
     &  gro2(NZ,NGROUP),                                                &
     &  gvrat(NBIN,NELEM,NGAS),                                         &
     &  supsatl(NX,NY,NZ,NGAS),                                         &
     &  supsati(NX,NY,NZ,NGAS),                                         &
     &  supsatlold(NXYZ,NGAS), supsatiold(NXYZ,NGAS),                   &
     &  scrit(NZ,NBIN,NGROUP),                                          &
     &  sol_ions(NSOLUTE), solwtmol(NSOLUTE), rhosol(NSOLUTE),          &
     &  rlh_nuc(NELEM,NELEM), T0(NGAS), Tfreez(NGAS),                   &
     &  pratt(3,NBIN,NGROUP), prat(4,NBIN,NGROUP), pden1(NBIN,NGROUP),  &
     &  palr(4,NGROUP),                                                 &
     &  iaer10
!
!
!  Declare alias names for pvapl(), pvapi(), supsatl(), and supsati() with
!  first 2, 3 dimensions treated linearly
!
      dimension pvapl2(NXY,NZ,NGAS), pvapl3(NXYZ,NGAS)
      dimension pvapi2(NXY,NZ,NGAS), pvapi3(NXYZ,NGAS)
      dimension supsatl2(NXY,NZ,NGAS)
      dimension supsatl3(NXYZ,NGAS)
      dimension supsati2(NXY,NZ,NGAS)
      dimension supsati3(NXYZ,NGAS)
      dimension supsatlold3(NXYZ,NGAS)
      dimension supsatiold3(NXYZ,NGAS)
!
      equivalence( pvapl2, pvapl )
      equivalence( pvapl3, pvapl )
      equivalence( pvapi2, pvapi )
      equivalence( pvapi3, pvapi )
      equivalence( supsatl2, supsatl )
      equivalence( supsatl3, supsatl )
      equivalence( supsati2, supsati )
      equivalence( supsati3, supsati )
      equivalence( supsatiold3, supsatiold )
      equivalence( supsatlold3, supsatlold )
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for nucleation parameters
!
!   ct       Contact parameter
!   adelf    Coefficient for phase change activation energy  {freznuc}
!   bdelf    Coefficient for phase change activation energy  {freznuc}
!   prenuc   Pre-exponential factore for frezzing nucleation {freznuc}
!   rmiv     Contact angle for ice/nucleus interface         {freznuc}
!   iaer11   Safety marker for common block rad1 
!
      common /aer11/                                                    &
     &  ct(NGAS,NGROUP),adelf, bdelf, prenuc, rmiv,                     &
     &  iaer11
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!  Declare global common blocks for extrapolation parameters
!
!   gcb         Gas concentration at beginning of extrapolation time
!   tbegin      Time at beginning of extrapolation time
!   tinc        Time increment for extrapolation over                 {init}
!   trecov      Time increment between each extrapolation             {init}
!   tlex        Time of last extrapolation                            {init}
!   prt_year    If .true. then print data that is only output yearly
!   iaer12     Safety marker for common block aer12
!
      logical prt_year

      common /aer12/                                                    &
     &   gcb(NXYZ,NGAS),tbegin,tinc,trecov,tlex,prt_year,               &
     &   iaer12
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare global common blocks for surface evaporation parameters
!
!   puddle      [g/cm2] of liquid at surface                        {postep}
!   fbotevap    Flux of igas into bottom layer from evaporation   {vertical}
!   ivolelem    Index array for volatile elements corresponding to igas
!   nvolelem    Number of volatile elements corresponding to igas {setupaer}
!
      common /aer13/                                                    &
     &   puddle(NGAS),fbotevap(NGAS),ivolelem(NGROUP-1,NGAS),           &
     &   nvolelem(NGAS),                                                &
     &   iaer13
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!  Declare global common blocks for radiative transfer parameters
!
!   do_rad      If .true. then do radiative transfer                  {init}
!   do_solar    If .true. then do solar calculations                  {init}
!   do_ir       If .true. then do infrared calculations               {init}
!   nrad        # of timesteps between radiation calcs (used when > 0)
!   prad        Time period between radiation calcs (used when nrad < 0)
!   isolar_zen  =I_FIXED: fixed, =I_DIURNAL: calculated every time step
!   u0          cos( solar_zenith_angle )
!   u0_fixed    Fixed value of cos( solar_zenith_angle )
!   rad_start   solar time corresponding to <time> = 0 [s]
!   zsin,zcos   sin and cos terms for solar zenith angle precalculation
!   wave        Bin-center wavelengths [micron]
!   radheat     Radiative heating rate [deg_K/s]
!   qrad        Particle heating rate [deg_K/s]
!   alb_tomi    Spectrally-integrated albedo at top-of-model
!   alb_toai    Spectrally-integrated albedo at top-of-atmosphere
!   alb_toa     Spectrally-resolved albedo at top-of-atmosphere
!   opd         Spectrally-resolved optical depth
!   fsl_up      Solar upwelling flux [W m^-2]
!   fsl_dn      Solar downwelling flux [W m^-2]
!   fir_up      Infrared upwelling flux [W m^-2]
!   fir_dn      Infrared downwelling flux [W m^-2]
!   irad1       Safety marker for common block rad1 
!
      logical do_rad, do_solar, do_ir

      common /rad1/                                                     &
     &  u0, u0_fixed, rad_start,                                        &
     &  zsin(NX,NY), zcos(NX,NY), wave(NWAVE),                          &
     &  qrad(NX,NY,NZ,NBIN,NGROUP), radheat(NX,NY,NZ),                  &
     &  alb_tomi(NX,NY), alb_toai(NX,NY),                               &
     &  alb_toa(NX,NY,NSOL), opd(NX,NY,NWAVE),                          &
     &  fsl_up(NX,NY,NZRADP1), fsl_dn(NX,NY,NZRADP1),                   &
     &  fir_up(NX,NY,NZRADP1), fir_dn(NX,NY,NZRADP1), prad,             &
     &  do_rad, do_solar, do_ir, nrad, isolar_zen,                      &
     &  irad1
!
!  EJL - here is where I need to put any variables I add to RT
!
!  Declare alias names for radiative transfer parameters
!  with first 2, 3 dimensions treated linearly
!
      dimension radheat2(NXY,NZ), radheat3(NXYZ)
      dimension qrad2(NXY,NZ,NBIN,NGROUP), qrad3(NXYZ,NBIN,NGROUP)
      dimension alb_tomi2(NXY)
      dimension alb_toai2(NXY)
      dimension alb_toa2(NXY,NSOL)
      dimension opd2(NXY,NWAVE)
      dimension fsl_up2(NXY,NZRADP1), fsl_up3(NXYZRADP1)
      dimension fsl_dn2(NXY,NZRADP1), fsl_dn3(NXYZRADP1)
      dimension fir_up2(NXY,NZRADP1), fir_up3(NXYZRADP1)
      dimension fir_dn2(NXY,NZRADP1), fir_dn3(NXYZRADP1)
!
      equivalence( radheat2, radheat )
      equivalence( radheat3, radheat )
      equivalence( qrad2, qrad )
      equivalence( qrad3, qrad )
      equivalence( alb_tomi2, alb_tomi )
      equivalence( alb_toai2, alb_toai )
      equivalence( alb_toa2, alb_toa )
      equivalence( opd2, opd )
      equivalence( fsl_up2, fsl_up )
      equivalence( fsl_up3, fsl_up )
      equivalence( fsl_dn2, fsl_dn )
      equivalence( fsl_dn3, fsl_dn )
      equivalence( fir_up2, fir_up )
      equivalence( fir_up3, fir_up )
      equivalence( fir_dn2, fir_dn )
      equivalence( fir_dn3, fir_dn )

!
!
!  Some extra definitions to use with error messages
!
      common /debug1/                                                   &
     &  isubstep,                                                       &
     &  parent_gc(NXYZ, NGAS), parent_gcl(NXYZ, NGAS),                  &
     &  parent_t(NXYZ),                                                 &
     &  parent_supsati(NXYZ,NGAS), parent_supsatl(NXYZ,NGAS),           &
     &  parent_supsatiold(NXYZ,NGAS), parent_supsatlold(NXYZ,NGAS)
     
    
