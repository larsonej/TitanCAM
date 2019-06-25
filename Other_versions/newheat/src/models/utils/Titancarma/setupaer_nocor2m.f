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
      include 'globaer.h'
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer_multi'
c
c-------------------------------------------------------------------------------
c
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Sample setup of particle types and the corresponding mapping arrays:
c
c    CN: sulfate particles
c    Cloud hydrometeors: liquid water drops and frogs
c    Core masses: mass of sulfate in water drops and mass of worms in frogs
c
c    group 1   sulfate particles                  element 1   (#/cm^3)
c    group 2   liquid drops                       element 2   (#/cm^3)
c              mass of sulfate in cloud drops     element 3   (g/cm^3)
c    group 3   frogs                              element 4   (#/cm^3)
c              mass of worms in frogs             element 5   (g/cm^3)
c
c    NGROUP  = 3
c    NELEM   = 5
c    nelemg  = (1,2,2)
c    itype   = (0,1,2,1,2)
c    ienconc = (1,2,4)
c    igelem  = (1,2,2,3,3)
c
c
c  Name for each group
c
      groupname(1) = 'Aerosol group'
      groupname(2) = 'Hydrocarbon Cloud group'
c
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 3
c
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'C2H6 ice crystal'
      elemname(3) = 'tholin core mass'
      elemname(4) = 'CH4 core mass (growcore)'
c
c
c  This array specifies the composition of each element:
c     I_H2SO4    is sulfuric acid  
c     I_WATER    is liquid water   
c     I_ICE      is ice water      
c     I_MIXEDWAT is mixed phase (ice/liquid) water
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_C2H6_ICE
      icomp(3) = I_CxHyNz
      icomp(4) = I_CH4_ICE
c
c
c  This array specifies the type of each element:
c     I_INVOLATILE is CN (involatile) number concentration [#/cm^3]
c     I_VOLATILE   is water droplet (volatile) number concentration [#/cm^3]
c     I_COREMASS   is core mass concentration [g/cm^3]
c     I_VOLCORE    is core mass concentration [g/cm^3] of a volatile core
c     I_CORE2MOM   is core second moment (mass^2) [g^2/cm^3]
c
      itype(1) = I_INVOLATILE
      itype(2) = I_VOLATILE
      itype(3) = I_COREMASS
      itype(4) = I_GROWCORE
c
c
c  Mass density for each particle element.  For elements of
c  <itype> = I_VOLATILE, <rhoelem> is the density of the shell.
c
      rhoelem(1) = RHO_I
      rhoelem(2) = RHO_I
      rhoelem(3) = RHO_I
      rhoelem(4) = 0.519d0
c
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3e-7
      rmin(2) = 1.048e-5
c
c
c  Ratio of particle mass between successive bins (one for each group)
c
      rmrat(1) = 2.0
      rmrat(2) = 2.0
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
c
c    <eshape> = particle length/diameter
c
      eshape(1) = 1.
      eshape(2) = 1.
c
c
c  Evaluate derived bin mapping arrays and set up the particle size bins.
c
      call setupbins
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
      vf_const = 0.0
c
c
c  Evaluate fall velocities.
c
      call setupvf
c
c
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c
c  Names of gas species
c
      gasname(1) = 'ethane'
      gasname(2) = 'methane'
c     gasname( ) = 'CH4 + N2'
c
c
c  Molecular weights of gas species
c
      gwtmol(1) = 30.07  !ethane
      gwtmol(2) = 16.04  !methane
c     gwtmol( ) = 20.24  !mixture of CH4 and N2
c
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Set to zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
c  *** Condensational growth of cores is not presently treated. ***
c
      igrowgas(1) = 0
      igrowgas(2) = 1
      igrowgas(3) = 0
      igrowgas(4) = 2
c
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .true.
c
c
c  If <is_grp_mixed> = .true. then the particle group is a mixed ice/liquid
c  hydrometeor.  This array is used to select processes (such as core
c  melting) that only occur in mixed particles.
c
      is_grp_mixed(1) = .false.
      is_grp_mixed(2) = .false.
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
      inucgas(2) = 0
c
c
c  Nucleation mapping:
c
c  <nnuc2elem(ielem)> is the total number of target elements from which
c  source particles in element <ielem> are nucleated.
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c
c  Set <nnuc2elem> to zero if mass is not transferred by nucleation
c  (in which case <inuc2elem> is automatically set to zero elsewhere).
c
      nnuc2elem(1) = 1
      nnuc2elem(2) = 0
      nnuc2elem(3) = 0
      nnuc2elem(4) = 0

      inuc2elem(1,1) = 3
c
c
c  <if_nuc(iefrom,ieto)> is true if nucleation transfers mass from element
c  <iefrom> to element <ieto>.
c
      do iefrom = 1,NELEM
        do ieto = 1,NELEM
          if_nuc(iefrom,ieto) = .false.
        enddo
      enddo
      if_nuc(1,3) = .true.
c
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c   I_DROPACT:  Aerosol activation to droplets
c   I_AERFREEZE: Aerosol homogeneous freezing
c   I_DROPFREEZE: Droplet homogeneous freezing
c   I_MIXEDFREEZE: Mixed total freezing
c   I_MIXEDMELT: Mixed total melting
c  Set to zero if particles are not subject to nucleation.
c
      do iefrom = 1,NELEM
        do ieto = 1,NELEM
          inucproc(iefrom,ieto) = 0
        enddo
      enddo
      inucproc(1,3) = I_DROPACT
c
c  Initialize nucleation update time [s].
c
      time_nuc = 0.
c
c
c  Initialize nucleation update time interval [s].
c
      period_nuc = 30.
c
c
c  Inititalize index of smallest bin in each group from which a particle
c  has nucleated.
c
      do igroup = 1,NGROUP
        inucmin(igroup) = NBIN
      enddo
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
      ievp2elem(3) = 1
      ievp2elem(4) = 0
c
c
c==Set up solute properties and mapping arrays.
c
c
c  Solute name (one for each solute)
c
c     solname(1) = 'sulfuric acid'
c
c
c  Solute molecular weights (one for each solute)
c
c     solwtmol(1) = 98.
c
c
c  Solute mass densities (one for each solute) [g/cm^3]
c
c     rhosol(1) = 1.38
c
c
c  Number of ions that solute dissociates into (one for each solute)
c
c     sol_ions(1) = 2.
c
c
c  Solute mapping: particle element <ielem> is composed of solute
c  <isolelem(ielem)>.  Should only be non-zero for elements of
c  itype = I_INVOLATILE [involatile number concentration] or 
c  itype = I_COREMASS [core mass concentration]).
c
c     isolelem(1) = 1
c     isolelem(2) = 0
c     isolelem(3) = 1
c     isolelem(4) = 0
c     isolelem(5) = 1
c     isolelem(6) = 0
c     isolelem(7) = 0
c     isolelem(8) = 1
c
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for condensational growth and evaporation.
c
      call setupgrow
      call setupgkern
c
c
c  <if_nuc_lh(iefrom,ieto)> is true if nucleation from element
c  <iefrom> to element <ieto> releases latent heat.
c
      do iefrom = 1,NELEM
        do ieto = 1,NELEM
          if_nuc_lh(iefrom,ieto) = .false.
        enddo
      enddo
c
c
c  <rlh_nuc(iefrom,ieto)> is the latent heat released by nucleation
c  from element <iefrom> to element <ieto> [cm^2/s^2].
c
      do iefrom = 1,NELEM
        do ieto = 1,NELEM
          rlh_nuc(iefrom,ieto) = 0.
        enddo
      enddo
c
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for nucleation.
c
      call setupnuc
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
        ck0 = 1.e-5
c
c
c  <grav_e_coll0> only used when <icollec> = 0
c
        grav_e_coll0 = 1.
c
c
c  This <NGROUP> by <NGROUP> array maps aerosol groups for coagulation
c  <icoag(i,j)> is the particle group resulting from a collision
c  between particles in groups <i> and <j>.
c  [This array must be diagonal, so the user need only specify
c  elements (i=1:NGROUP,j=1:i). ]
c
        icoag(1,1) = 1
        icoag(2,1) = 2
        icoag(2,2) = 2
c
c
c  Evaluate derived coagulation mapping arrays and kernels.
c
        call setupckern
        call setupcoag

      endif
c
c
c-------------------------------------------------------------------------------
c
c
c==Set up mass production values:
c
c
        if(zl3(NZ) .gt. 100.d5)
     $    call setupmprod
c
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined.
c
      return
      end
