       subroutine setupaer
c
c
c  @(#) setupaer.f  Ackerman  Jan-1996
c  This master routine sets up user-defined mapping arrays and parameters
c  and calls all the other setup routines to calculate other time-independent
c  parameters for aerosol and cloud microphysics.
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
      include 'globals.h'
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer'
c
c-------------------------------------------------------------------------------
c
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Sample setup of particle types and the corresponding mapping arrays:
c
c    CN: sulfate particles, salt particles, mixed salt/sulfate particles
c    Cloud drops: liquid water drops, ice crystals
c    Core masses: mass of sulfate in water drops, salt in water drops
c
c    group 1   sulfate particles                  element 1   (#/cm^3)
c    group 2   salt particles                     element 2   (#/cm^3)
c    group 3   mixed salt/sulfate particles       element 3   (#/cm^3)
c              mass of salt in mixed particles    element 4   (g/cm^3)
c    group 4   water drops                        element 5   (#/cm^3)
c              mass of sulfate in water drops     element 6   (g/cm^3)
c              mass of salt in water drops        element 7   (g/cm^3)
c    group 5   ice crystal                        element 8   (#/cm^3)
c
c    ngroups = 5
c    nelem   = 8
c    nelemg  = (1,1,2,3,1)
c    itype   = (0,0,0,2,1,2,2,1)
c    ienconc  = (1,2,3,5,8)
c    igelem  = (1,2,3,3,4,4,4,5)
c
c
c  Number of particle mass bins 
c
      nbins = 10
c
c
c  Number of particle groups (each group represents a different
c  type of particle)
c
      ngroups = 5
c
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 1
      nelemg(3) = 1
      nelemg(4) = 3
      nelemg(5) = 4
c
c
c  This array specifies the type of each element:
c      0 is CN (involatile) number concentration [#/cm^3]
c      1 is water droplet (volatile) number concentration [#/cm^3]
c      2 is core mass concentration [g/cm^3]
c      3 is core second moment (mass^2) [g^2/cm^3]
c
      itype(1) = 0
      itype(2) = 0
      itype(3) = 0
      itype(4) = 0
      itype(5) = 2
      itype(6) = 2
      itype(7) = 1
      itype(8) = 2
      itype(9) = 2
      itype(10) = 2
c
c
c  Mass density for each particle element.  For elements of <itype> = 1,
c  <rhoelem> is the density of the shell.
c
      rhoelem(1) = 1.
      rhoelem(2) = 1.
      rhoelem(3) = 1.
      rhoelem(4) = 1.
      rhoelem(5) = 1.
      rhoelem(6) = 1.
      rhoelem(7) = 1.
      rhoelem(8) = 1.
      rhoelem(9) = 1.
      rhoelem(10) = 1.
c
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.e-4 
      rmin(2) = 1.e-4 
      rmin(3) = 1.e-4 
      rmin(4) = 1.e-4 
      rmin(5) = 1.e-4 
c
c
c  Ratio of particle mass between successive bins (one for each group)
c
      rmrat(1) = 2.0
      rmrat(2) = 2.0
      rmrat(3) = 2.0
      rmrat(4) = 2.0
      rmrat(5) = 2.0
c
c
c  The values of <ishape> and <eshape> determine particle geometry
c  (one for each group):
c
c    <ishape> = 1: spherical
c    <ishape> = 2: hexagonal prisms or plates
c    <ishape> = 3: circular disks, cylinders, or spheroids
c
      ishape(1) = 1
      ishape(2) = 1
      ishape(3) = 1
      ishape(4) = 1
      ishape(5) = 1
c
c    <eshape> = particle length/diameter
c
      eshape(1) = 1.
      eshape(2) = 1.
      eshape(3) = 1.
      eshape(4) = 1.
      eshape(5) = 1.
c
c
c  Evaluate derived bin mapping arrays and set up the particle size bins.
c
      call setupbins(rmin)
c
c-------------------------------------------------------------------------------
c
c
c==Set options for particle fall velocities.
c
c
c    <ifall> = 0: use constant fall velocity <vf_const>, but still calculate
c                 <vf> for use in Reynolds' number <re>
c            = 1: use calculated value of <vf>
c
      ifall = 1
c
c
c  <vf_const> is only used when <ifall> = 0
c
      vf_const = 1.0
c
c
c  Evaluate fall velocities.
c
      call setupvfall
c
c-------------------------------------------------------------------------------
c
c
c==Set options for coagulation kernel:
c
      if( do_coag )then
c
c
c   <icoagop> = 0: use fixed coagulation kernel <ck0>
c             = 1: calculate the coagulation kernel
c
        icoagop = 1
c
c
c   <icollec> determines gravitational collection efficiencies
c     = 0: use constant value <grav_e_coll0>
c     = 1: use binwise maxima of Fuchs' and Langmuir's efficiencies
c     = 2: use input data
c
        icollec = 2
c
c
c  Fixed kernel value <ck0> only used when <icoagop> = 0
c
        ck0 = 1.
c
c
c  <grav_e_coll0> only used when <icollec> = 0
c
        grav_e_coll0 = 1.
c
c
c  This <ngroups> by <ngroups> array maps aerosol groups for coagulation
c  <icoag(i,j)> is the particle group resulting from a collision
c  between particles in groups <i> and <j>.
c  [This array must be diagonal, so the user need only specify
c  elements (i=1:ngroups,j=1:i). ]
c
c  *** Coagulation does not yet treat multi-cores. ***
c
        icoag(1,1) = 1
        icoag(2,1) = 4
        icoag(2,2) = 2
        icoag(3,1) = 4
        icoag(3,2) = 4
        icoag(3,3) = 4
        icoag(4,1) = 4
        icoag(4,2) = 4
        icoag(4,3) = 4
        icoag(4,4) = 4
        icoag(5,1) = 5
        icoag(5,2) = 5
        icoag(5,3) = 5
        icoag(5,4) = 5
        icoag(5,5) = 5

c
c  Evaluate derived coagulation mapping arrays and kernels.
c
        call setupckern
        call setupcoag

      endif
c
c-------------------------------------------------------------------------------
c
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
      if( do_grow )then
c
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c
c  Total number of gases used in simulation
c
        ngas = 1
c
c
c  Names of gas species
c
        gasname(1) = 'water vapor'
c
c
c  Molecular weights of gas species
c
        gwtmol(1) = 18.
c
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Set to zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
        igrowgas(1) = 0
        igrowgas(2) = 0
        igrowgas(3) = 0
        igrowgas(4) = 0
        igrowgas(5) = 0
        igrowgas(6) = 0
        igrowgas(7) = 1
        igrowgas(8) = 0
        igrowgas(9) = 0
        igrowgas(10) = 0
c
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
        is_grp_ice(1) = .false.
        is_grp_ice(2) = .false.
        is_grp_ice(3) = .false.
        is_grp_ice(4) = .false.
        is_grp_ice(5) = .false.
c
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Set to zero if particles are not subject to nucleation.
c
      inucgas(1) = 1
      inucgas(2) = 1
      inucgas(3) = 1
      inucgas(4) = 1
      inucgas(5) = 0
c
c
c  Nucleation mapping: nucleation transfers particle mass from element <ielem>
c  to element <inuc2elem(ielem)>.
c  Set to zero if mass is not transferred by nucleation.
c
      inuc2elem(1) = 8
      inuc2elem(2) = 9
      inuc2elem(3) = 10
      inuc2elem(4) = 8
      inuc2elem(5) = 9
      inuc2elem(6) = 10
      inuc2elem(7) = 0
      inuc2elem(8) = 0
      inuc2elem(9) = 0
      inuc2elem(10) = 0
c
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Set to zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(1) = 0
      ievp2elem(2) = 0
      ievp2elem(3) = 0
      ievp2elem(4) = 0
      ievp2elem(5) = 0
      ievp2elem(6) = 0
      ievp2elem(7) = 0
      ievp2elem(8) = 4
      ievp2elem(9) = 5
      ievp2elem(10) = 6
c
c
c==Set up solute properties and mapping arrays.
c
c
c  Total number of solutes used in simulation
c
        nsolute = 2
c
c
c  Solute name (one for each solute)
c
        solwtmol(1) = 115.
        solwtmol(2) = 98.
c
c
c  Solute molecular weights (one for each solute)
c
        solname(1) = 'ammonium bisulfate'
c       solname(2) = 'sulfuric acid'
c
c
c  Solute mass densities (one for each solute) [g/cm^3]
c
        rhosol(1) = 1.78
        rhosol(2) = 1.38
c
c
c  Number of ions that solute dissociates into (one for each solute)
c
        sol_ions(1) = 2.
        sol_ions(2) = 2.
c
c
c  Solute mapping: particle element <ielem> is composed of solute
c  <isolelem(ielem)> (should only be non-zero for elements of
c  itype = 0 or 2).
c
        isolelem(1) = 1
        isolelem(2) = 1
        isolelem(3) = 2
        isolelem(4) = 1
        isolelem(5) = 1
        isolelem(6) = 2
        isolelem(7) = 0
        isolelem(8) = 1
        isolelem(9) = 1
        isolelem(10) = 2
c
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for condensational growth, evaporation, and nucleation.
c
        call setupgrow
        call setupgkern
        call setupsol
        call setupnuc

      endif
c
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined.
c
      return
      end
