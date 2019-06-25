       subroutine setupaer
c
c
c  @(#) setupaer.f  Ackerman  Jan-1996
c  This master routine sets up user-defined mapping arrays and parameters
c  and calls all the other setup routines to calculate other time-independent
c  parameters for aerosol and cloud microphysics.
c
c  ELB Aug 2003
c  Modified to include all possible setups by calling different setupaer
c  subroutines (similar to what is done in {initatm})
c  ELB Jan 2004
c  Condensed by moving common setups to front, and only calling rest of setup
c  routines once.
c
c  (1) List some sample setup values and complete options for some variables
c  (2) Initialize variables that don''t generally change for each case
c      (eg. rmrat, ishape, etc.), but could still be changed later if necessary
c  (3) Call appropriate setupaer routine
c  (4) Call other subroutines (ie. setupvf, setupnuc, etc.)
c
c
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer'
c
c-------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------
c 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
c-------------------------------------------------------------------------------
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
c  --------------------------------------------------------------------
c
c  Element composition <icomp> options:
c     I_H2SO4    is sulfuric acid  
c     I_WATER    is liquid water   
c     I_ICE      is ice water      
c     I_MIXEDWAT is mixed phase (ice/liquid) water
c     I_CxHyNz   is tholin
c     I_C2H6_ICE is ethane ice
c     I_CH4_ICE  is methane ice
c     I_CH4_LIQ  is liquid methane 
c
c  Nucleation process <inucproc> options:
c   I_DROPACT:  Aerosol activation to droplets
c   I_AERFREEZE: Aerosol homogeneous freezing
c   I_DROPFREEZE: Droplet homogeneous freezing
c   I_MIXEDFREEZE: Mixed total freezing
c   I_MIXEDMELT: Mixed total melting
c   I_VAPORDEP: Heterogeneous nucl by vapor deposition
c
c
c-------------------------------------------------------------------------------
c 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
c-------------------------------------------------------------------------------
c  Initializations and set some variables that don''t vary from case to case
c  Note that these can still be changed later (in setupaer_ ...) if necessary
c
c  Ratio of particle mass between successive bins (one for each group)
c
      do ig = 1,NGROUP
        rmrat(ig) = 2.0
      enddo
c
c  The values of <ishape> and <eshape> determine particle geometry
c  (one for each group):
c
c    <ishape> = 1: spherical
c    <ishape> = 2: hexagonal prisms or plates
c    <ishape> = 3: circular disks, cylinders, or spheroids
c
c    <eshape> = particle length/diameter
c
      do ig = 1,NGROUP
        ishape(ig) = 1
        eshape(ig) = 1.
      enddo
c
c
c-------------------------------------------------------------------------------
c==Set options for particle fall velocities.
c
c    <ifall> = 0: use constant fall velocity <vf_const>, but still calculate
c                 <vf> for use in Reynolds'' number <re>
c            = 1: use calculated value of <vf>
c
      ifall = 1
c
c  <vf_const> is only used when <ifall> = 0
c
      vf_const = 1.d0
c-------------------------------------------------------------------------------
c
c  Initialize mapping arrays for condensational growth, evaporation, and nucleation.
c  (Descriptions of these arrays are given with each setupaer_... routine)
c
      do ielem = 1,NELEM
        igrowgas(ielem) = 0
        ievp2elem(ielem) = 0
      enddo
 
      do igroup = 1,NGROUP
        do igas = 1,NGAS
          inucgas(igas,igroup) = 0
        enddo
      enddo
 
      do i = 1,5           !5 is array size defined in {globaer.h}
        do iefrom = 1,NELEM
          inuc2elem(i,iefrom) = 0
        enddo
      enddo
 
      do iefrom = 1,NELEM
        do ieto = 1,NELEM
          inucproc(iefrom,ieto) = 0
        enddo
      enddo
c
c  Initialize nucleation update time [s].
c
      time_nuc = endtime
c
c  Initialize nucleation update time interval [s].
c
      period_nuc = endtime
c
c  Inititalize index of smallest bin in each group from which a particle
c  has nucleated.
c
      do igroup = 1,NGROUP
        inucmin(igroup) = NBIN
      enddo
c
c==Set up solute properties and mapping arrays.
c   (Note: All solute code is here since not using it for Titan models.
c    Will need to move some of this back to individual setupaer... routines
c    if want to include solutes)
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
c  itype = I_COREMASS [core mass concentration].
c
c     isolelem(1) = 1
c     isolelem(2) = 0
c
c
c  (Note: No latent heat release from nucleation.  Will need to move some of 
c   this back into individual setupaer_... routines if want to model latent
c   heat release from nucleation)
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
c-------------------------------------------------------------------------------
c
c==Set options for coagulation kernel:
c
c   <icoagop> = 0: use fixed coagulation kernel <ck0>
c             = 1: calculate the coagulation kernel
c
        icoagop = 1
c
c   <icollec> determines gravitational collection efficiencies
c     = 0: use constant value <grav_e_coll0>
c     = 1: use binwise maxima of Fuchs' and Langmuir's efficiencies
c     = 2: use input data
c
        icollec = 2
c
c  Fixed kernel value <ck0> only used when <icoagop> = 0
c
        ck0 = 1.e-5
c
c  <grav_e_coll0> only used when <icollec> = 0
c
        grav_e_coll0 = 1.
c
c-------------------------------------------------------------------------------
c 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
c-------------------------------------------------------------------------------
c
c  Choose appropriate {setupaer} routine based on number of elements
c
c    Subroutine name specifies types of clouds in simulation
c    (e.g. ethane, methane, mixed, water) and phase of condensate in
c    cloud (ice, liquid, both=2phase).  Mixed clouds are methane ice
c    crystal or droplet with an ethane growcore.  Mixed cloud models
c    always include pure ethane clouds.
c
       if( NELEM .eq. 2 ) then
         call setupaer_haze
       elseif( NELEM .eq. 4 ) then
c        call setupaer_ethane_ice
         call setupaer_methane_ice
c        call setupaer_water_ice
       else if( NELEM .eq. 5 ) then
         call setupaer_methane_volcore
       else if( NELEM .eq. 7 ) then
         call setupaer_methane_2phase
       else if( NELEM .eq. 8 ) then
         call setupaer_mixed_ice
       else if( NELEM .eq. 11 ) then
        !call setupaer_methane_ice_mixed_ice
         call setupaer_mixed_volcore
       else if( NELEM .eq. 12 ) then
         call setupaer_mixed_2phase
       else
         write(*,*) 'Wrong number of elements',NELEM
         stop
       endif
c
c-------------------------------------------------------------------------------
c 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
c-------------------------------------------------------------------------------
c
c  Automatically set some variables used in the subroutines called below.
c
c  Mass density for each particle element.  For elements of
c  <itype> = I_VOLATILE, <rhoelem> is the density of the shell.
c
      do ielem = 1,NELEM
        if( icomp(ielem) .eq. I_CxHyNz ) then
         !rhoelem(ielem) = 0.713d0
          rhoelem(ielem) = 1.0d0 !0.713d0
          write(*,*) 'setupaer: setting tholin rhoelem=1'
        elseif( icomp(ielem) .eq. I_C2H6_ICE ) then
          rhoelem(ielem) = 0.713d0
        elseif( icomp(ielem) .eq. I_C2H6_LIQ ) then
          rhoelem(ielem) = 0.5446d0
        elseif( icomp(ielem) .eq. I_CH4_ICE ) then
          rhoelem(ielem) = 0.519d0
        elseif( icomp(ielem) .eq. I_CH4_LIQ ) then
          rhoelem(ielem) = 0.4228d0
        elseif( icomp(ielem) .eq. I_ICE ) then
          rhoelem(ielem) = 0.930d0
        else
          write(*,*) '<setupaer> Undefined element composition'
          stop
        endif
      enddo
c
c  Molecular weights of gas species
c
      do igas = 1,NGAS
        if( gasname(igas) .eq. 'ethane' ) then
          gwtmol(igas) = 30.07 
        elseif( gasname(igas) .eq. 'methane' ) then
          gwtmol(igas) = 16.04
        elseif( gasname(igas) .eq. 'nitrogen' ) then
          gwtmol(igas) = 28.01
        elseif( gasname(igas) .eq. 'water' ) then
          gwtmol(igas) = 18.02
        else
          write(*,*) '<setupaer> Undefined gasname'
          stop
        endif
      enddo
c    
c-------------------------------------------------------------------------------
c 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
c-------------------------------------------------------------------------------
c
c  Call the rest of the setup subroutines
c
c  Evaluate derived bin mapping arrays and set up the particle size bins.
c
      call setupbins
c
c  Automatically set some variables used in the subroutines called below.
c  (Note: Must be called after {setupbins} because need <igelem> defined)
c
c  If <is_grp_mixed_phase> = .true. then the particle group is a mixed 
c  ice/liquid hydrometeor.  This array is used to select processes 
c  (such as core melting) that only occur in mixed-phase particles.
c
c  If <is_grp_mixed_comp> = .true. then the particle group is a cloud
c  composed of multiple volatile species
c 
      do igroup = 1,NGROUP
        is_grp_mixed_phase(igroup) = .false.
        is_grp_mixed_comp(igroup) = .false.
      enddo
      do ielem = 1,NELEM
        igroup = igelem(ielem)
        if( itype(ielem) .eq. I_VOLCORE ) then
                is_grp_mixed_phase(igroup) = .true.
          ! Also check that element preceeding VOLCORE is a numconc or
          ! growcore, as this is assumed in <upgxfer>
            if( .not. ( ( itype(ielem-1) .eq. I_VOLATILE ) .or.
     $                  ( itype(ielem-1) .eq. I_GROWCORE ) )   ) then
              write(*,*) '<setupaer> VOLCORE elements must follow their
     $                    associated droplet or growcore element'
              stop
            endif
        endif
        if( itype(ielem) .eq. I_GROWCORE ) 
     $          is_grp_mixed_comp(igroup) = .true.
      enddo
c
c  Evaluate fall velocities.
c
      call setupvf
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for condensational growth and evaporation.
c
      call setupgrow
      call setupgkern
c
c  Evaluate time-independent parameters and derived mapping arrays
c  used for nucleation.
c
      call setupnuc

cc Print nucleation mapping arrays to screen for debugging, then stop
c     write(*,*) 'Liquid clouds - nucleation mapping'
c     do ielem = 1,NELEM
c       write(*,*) 'ielem,nnucelem',ielem,nnucelem(ielem)
c       do jfrom = 1,nnucelem(ielem)
c        write(*,*) 'jfrom,inucelem =',jfrom,inucelem(jfrom,ielem)
c       enddo
c     enddo
cc    stop
c
c  Set up array of element indexes for volatile elements corresponding
c  to each gas: <ivolelem>
      do igas = 1,NGAS
        nvolelem(igas) = 0
        do n=1,NGROUP-1
          ivolelem(n,igas) = 0
        enddo
      enddo

      do igas = 1,NGAS
        n=1
        do ielem = 1,NELEM
          if( igrowgas(ielem) .eq. igas ) then
            ivolelem(n,igas) = ielem
            nvolelem(igas) = n
            n = n + 1

           !Don't include methane ice crystals if also droplets in this run
            if( icomp(ielem).eq.I_CH4_ICE .and. T0(igas).lt.81.) n = n-1
           !Don't include ethane growcore in mixed ice cloud if also mixed
           ! cloud droplets in this run
            if( itype(ielem).eq.I_GROWCORE .and. 
     $            is_grp_ice(igelem(ielem)) .and. T0(igas).lt.81.) n=n-1
          endif
        enddo !elements
      enddo !igas
        write(*,*) '# vol elements:',nvolelem(1),nvolelem(2)
c
c
c  Evaluate derived coagulation mapping arrays and kernels.
c
      if( do_coag ) then
        call setupckern
        call setupcoag
      endif
c
c  Evaluate aerosol mass production
c
      call setupmprod
c
c  Check particle setup for incompatibilities with other subroutines
c
      do igroup=1,NGROUP
        ncm = 0
        do ic=1,ncore(igroup)
          ie = icorelem(ic,igroup)
          if( itype(ie) .eq. I_COREMASS ) ncm = ncm + 1
        enddo
        if( ncm .gt. 1 ) then
          if( if_sec_mom(igroup) .and. is_grp_mixed_phase(igroup)) then
            write(LUNOPRT,*) '<setupaer> Particle group ',igroup,
     $         'incompatible with <evap_mono>.  Can only have one
     $          coremass element in a group with core 2nd moments
     $          and growcores'
            stop
          endif
        endif !more than one coremass element in group
      enddo
c
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined.
c
      return
      end
c
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c                        Begin specific case setupaer routines
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_haze
c                                                                              
c  Ethane clouds                                                              
c                                                                           
c  Include global constants and variables                                  
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer_haze'
      write(*,'(/,a)') 'Haze only model'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Name for each group
c
      groupname(1) = 'Aerosol Group'
      groupname(2) = 'Place Holder Group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 1
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = ' '
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_C2H6_ICE
c
c  This array specifies the type of each element:
c     I_INVOLATILE is CN (involatile) number concentration [#/cm^3]
c     I_VOLATILE   is water droplet (volatile) number concentration [#/cm^3]
c     I_COREMASS   is core mass concentration [g/cm^3]
c     I_CORE2MOM   is core second moment (mass^2) [g^2/cm^3]
c
      itype(1) = I_INVOLATILE
      itype(2) = VOLATILE
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3e-7
      rmin(2) = 1.3e-7
c
c  Defaults: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'ethane'
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .false.
c
c  Coagulation mapping:
c
c  This <NGROUP> by <NGROUP> array maps aerosol groups for coagulation
c  <icoag(i,j)> is the particle group resulting from a collision
c  between particles in groups <i> and <j>.
c  [This array must be diagonal, so the user need only specify
c  elements (i=1:NGROUP,j=1:i). ]
c
        icoag(1,1) = 1
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for the ethane cloud case.
c
      return
      end
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_ethane_ice
c                                                                              
c  Ethane clouds                                                              
c                                                                           
c  Include global constants and variables                                  
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer_ethane_ice'
      write(*,'(/,a)') '2 Groups - Ethane Clouds'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Name for each group
c
      groupname(1) = 'Aerosol Group'
      groupname(2) = 'Ethane Cloud Group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 3
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'C2H6 ice crystal'
      elemname(3) = 'tholin core mass'
      elemname(4) = '2nd moment of tholin mass in ice'
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_C2H6_ICE
      icomp(3) = I_CxHyNz
      icomp(4) = I_CxHyNz
c
c  This array specifies the type of each element:
c     I_INVOLATILE is CN (involatile) number concentration [#/cm^3]
c     I_VOLATILE   is water droplet (volatile) number concentration [#/cm^3]
c     I_COREMASS   is core mass concentration [g/cm^3]
c     I_CORE2MOM   is core second moment (mass^2) [g^2/cm^3]
c
      itype(1) = I_INVOLATILE
      itype(2) = I_VOLATILE
      itype(3) = I_COREMASS
      itype(4) = I_CORE2MOM
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3e-7
      rmin(2) = 2.08e-6  !1.048e-5
c
c  Defaults: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'ethane'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Leave this zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .true.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
      inucgas(1,1) = 1
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c
      inuc2elem(1,1) = 3
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c                                           (See top for more options)
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,3) = I_VAPORDEP
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(3) = 1
c
c  Coagulation mapping:
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
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for the ethane cloud case.
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_methane_ice
c
c  Methane clouds 
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer_methane_ice'
      write(*,'(/,a)') '2 Groups - Methane Clouds'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c
c  Name for each group
c
      groupname(1) = 'Aerosol Group'
      groupname(2) = 'Methane Cloud Group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 3
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'CH4 ice crystal'
      elemname(3) = 'tholin core mass'
      elemname(4) = '2nd moment of tholin mass in ice'
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice

      icomp(1) = I_CxHyNz
      icomp(2) = I_CH4_ICE
      icomp(3) = I_CxHyNz
      icomp(4) = I_CxHyNz
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
      itype(4) = I_CORE2MOM
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3e-7
      rmin(2) = 2.08e-6  !1.048e-5
c
c  Defaults:  <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'methane'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Leave this zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .true.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
      inucgas(1,1) = 1
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c
      inuc2elem(1,1) = 3
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c                                            (See top for more options)
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,3) = I_VAPORDEP
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(3) = 1
c
c-------------------------------------------------------------------------------
c
c  Coagulation mapping:
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
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for methane cloud case.
c
      return
      end
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_mixed_ice
c
c  Both ethane and methane nucleation: ethane clouds and methane clouds with
c  ethane growcores
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer_mixed_ice'
      write(*,'(/,a)') '3 Groups - Ethane & Mixed Clouds w/ gc'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
      if( NGROUP .ne. 3 ) then
        write(*,*) 'Reset <NGROUP> in {aerad}'
        stop 1 
      endif
      if( NELEM .ne. 8 ) then
        write(*,*) 'Reset <NELEM> in {globaer}'
        stop 1 
      endif
c
c  Name for each group
c
      groupname(1) = 'Aerosol group'
      groupname(2) = 'Ethane Cloud group'
      groupname(3) = 'Methane Cloud group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 3
      nelemg(3) = 4
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'C2H6 ice crystal'
      elemname(3) = 'tholin core mass'
      elemname(4) = '2nd moment of tholin mass in ice'
      elemname(5) = 'CH4 ice crystal'
      elemname(6) = 'tholin core mass' 
      elemname(7) = '2nd moment of tholin mass in ice'
      elemname(8) = 'C2H6 core mass (growcore)'
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice
c     I_CH4_ICE  is methane ice
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_C2H6_ICE
      icomp(3) = I_CxHyNz
      icomp(4) = I_CxHyNz
      icomp(5) = I_CH4_ICE
      icomp(6) = I_CxHyNz
      icomp(7) = I_CxHyNz
      icomp(8) = I_C2H6_ICE
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
      itype(4) = I_CORE2MOM
      itype(5) = I_VOLATILE
      itype(6) = I_COREMASS
      itype(7) = I_CORE2MOM
      itype(8) = I_GROWCORE
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3e-7
      rmin(2) = 2.08d-6  !1.048e-5
      rmin(3) = 2.62d-6  !1.320e-5

      rmin(1) = 2.08d-6  
      if( rmin(1) .gt. 1.5e-7 ) 
     $   write(*,*) 'SETUPAER: Larger tholin !!!!!!!!!!!'
c
c  Default values: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'ethane'
      gasname(2) = 'methane'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Set to zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1   !map ethane gas to ethane cloud
      igrowgas(5) = 2   !map methane gas to mixed cloud
      igrowgas(8) = 1   !map ethane gas to mixed cloud growcore
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .true.
      is_grp_ice(3) = .true.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
      inucgas(1,1) = 1
      inucgas(1,2) = 2
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c  (<nnuc2elem(ielem)> will be set in {setupnuc})
c
      inuc2elem(1,1) = 3
      inuc2elem(1,2) = 8
      inuc2elem(1,3) = 6
      inuc2elem(1,4) = 7
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c                                           (See top for more options)
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,3) = I_VAPORDEP
      inucproc(2,8) = I_VAPORDEP
      inucproc(3,6) = I_VAPORDEP
      inucproc(4,7) = I_VAPORDEP
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(3) = 1   !evap tholin core in ethane cloud to tholin
      ievp2elem(6) = 3   !evap tholin core in mixed cloud to tholin core in ethane cloud
      ievp2elem(7) = 4   !evap tholin core2m in mixed cloud to tholin core2m in ethane cloud 
      ievp2elem(8) = 2   !evap ethane growcore in mixed cloud to ethane cloud
c
c-------------------------------------------------------------------------------
c
c  Coagulation mapping:
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
        icoag(3,1) = 3
        icoag(3,2) = 3
        icoag(3,3) = 3
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for ethane clouds, methane
c  clouds with ethane growcores model.
c
      return
      end

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_mixed_2phase
c
c
c  Both ethane and methane nucleation, methane clouds can be ice or liquid
c  (liquid clouds include N2)
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer_mixed_2phase'
      write(*,'(/,a)') 'Allow Liquid CH4 in Mixed Clouds'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Name for each group
c
      groupname(1) = 'Aerosol group'
      groupname(2) = 'Ethane Cloud group'
      groupname(3) = 'Methane Ice Cloud group'
      groupname(4) = 'Methane Liq Cloud group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 3
      nelemg(3) = 4
      nelemg(4) = 4
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'C2H6 ice crystal'
      elemname(3) = 'tholin core mass'
      elemname(4) = '2nd moment of tholin mass in ice'
      elemname(5) = 'CH4 ice crystal'
      elemname(6) = 'tholin core mass' 
      elemname(7) = '2nd moment of tholin mass in ice'
      elemname(8) = 'C2H6 core mass (growcore)'
      elemname(9) = 'CH4 droplet'
      elemname(10) = 'tholin core mass' 
      elemname(11) = '2nd moment of tholin mass in drop'
      elemname(12) = 'C2H6 core mass (growcore)'
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice
c     I_CH4_ICE  is methane ice
c     I_CH4_LIQ  is liquid methane 
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_C2H6_ICE
      icomp(3) = I_CxHyNz
      icomp(4) = I_CxHyNz
      icomp(5) = I_CH4_ICE
      icomp(6) = I_CxHyNz
      icomp(7) = I_CxHyNz
      icomp(8) = I_C2H6_ICE
      icomp(9) = I_CH4_LIQ
      icomp(10) = I_CxHyNz
      icomp(11) = I_CxHyNz
      icomp(12) = I_C2H6_ICE
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
      itype(4) = I_CORE2MOM
      itype(5) = I_VOLATILE
      itype(6) = I_COREMASS
      itype(7) = I_CORE2MOM
      itype(8) = I_GROWCORE
      itype(9) = I_VOLATILE
      itype(10) = I_COREMASS
      itype(11) = I_CORE2MOM
      itype(12) = I_GROWCORE
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3e-7
      rmin(2) = 2.08d-6  !1.048e-5
      rmin(3) = 2.62d-6  !1.320e-5
      rmin(4) = 2.62d-6  !1.320e-5
c     rmin(4) = 3.30d-6  !2.62d-6  !1.320e-5
c
c  Defaults: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'ethane'
      gasname(2) = 'methane'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Leave this zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1
      igrowgas(5) = 2
      igrowgas(8) = 1
      igrowgas(9) = 2
      igrowgas(12) = 1
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .true.
      is_grp_ice(3) = .true.
      is_grp_ice(4) = .false.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
      inucgas(1,1) = 1
      inucgas(1,2) = 2
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c  (<nnuc2elem(ielem)> will be set in {setupnuc})
c
      inuc2elem(1,1) = 3
      inuc2elem(1,2) = 8
      inuc2elem(1,3) = 6
      inuc2elem(1,4) = 7
      inuc2elem(2,2) = 12
      inuc2elem(2,3) = 10
      inuc2elem(2,4) = 11
       !include melting
      inuc2elem(1,5) = 9
      inuc2elem(1,6) = 10
      inuc2elem(1,7) = 11
      inuc2elem(1,8) = 12
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c                                           (See top for more options)
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,3) = I_VAPORDEP
      inucproc(2,8) = I_VAPORDEP
      inucproc(3,6) = I_VAPORDEP
      inucproc(4,7) = I_VAPORDEP
      inucproc(2,12) = I_VAPORDEP
      inucproc(3,10) = I_VAPORDEP
      inucproc(4,11) = I_VAPORDEP
     !melting
      inucproc(5,9) = I_ICEMELT
      inucproc(6,10) = I_ICEMELT
      inucproc(7,11) = I_ICEMELT
      inucproc(8,12) = I_ICEMELT
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(3) = 1
      ievp2elem(6) = 3
      ievp2elem(7) = 4
      ievp2elem(8) = 2
      ievp2elem(10) = 3
      ievp2elem(11) = 4
      ievp2elem(12) = 2
c
c-------------------------------------------------------------------------------
c
c  Coagulation mapping:
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
        icoag(3,1) = 3
        icoag(3,2) = 3
        icoag(3,3) = 3
        icoag(4,1) = 4
        icoag(4,2) = 4
        icoag(4,3) = 0  !!no mixing liq/ice mixed clouds
        icoag(4,4) = 4
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for ethane clouds, methane 
c  clouds with ethane growcores - where methane clouds can be ice or liquid.
c
      return
      end
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_methane_ice_mixed_ice
c
c  Both ethane and methane nucleation: ethane clouds and methane clouds with
c  ethane growcores.  Additionally, allow methane nucleation onto tholin
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 
     &    'Enter setupaer_methane_ice_mixed_ice'
      write(*,'(/,a)') '4 Groups - Allow Pure Methane Clouds'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Name for each group
c
      groupname(1) = 'Aerosol group'
      groupname(2) = 'Ethane Cloud group'
      groupname(3) = 'Pure Methane Cloud group'
      groupname(4) = 'Methane Cloud group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 3
      nelemg(3) = 3
      nelemg(4) = 4
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'C2H6 ice crystal'
      elemname(3) = 'tholin core mass'
      elemname(4) = '2nd moment of tholin mass in ice'
      elemname(5) = 'CH4 ice crystal'
      elemname(6) = 'tholin core mass' 
      elemname(7) = '2nd moment of tholin mass in ice'
      elemname(8) = 'CH4 ice crystal'
      elemname(9) = 'tholin core mass' 
      elemname(10) = '2nd moment of tholin mass in ice'
      elemname(11) = 'C2H6 core mass (growcore)'
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice
c     I_CH4_ICE  is methane ice
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_C2H6_ICE
      icomp(3) = I_CxHyNz
      icomp(4) = I_CxHyNz
      icomp(5) = I_CH4_ICE
      icomp(6) = I_CxHyNz
      icomp(7) = I_CxHyNz
      icomp(8) = I_CH4_ICE
      icomp(9) = I_CxHyNz
      icomp(10) = I_CxHyNz
      icomp(11) = I_C2H6_ICE
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
      itype(4) = I_CORE2MOM
      itype(5) = I_VOLATILE
      itype(6) = I_COREMASS
      itype(7) = I_CORE2MOM
      itype(8) = I_VOLATILE
      itype(9) = I_COREMASS
      itype(10) = I_CORE2MOM
      itype(11) = I_GROWCORE
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3d-7
      rmin(2) = 2.08d-6  !1.048e-5
      rmin(3) = 2.08d-6  !1.048e-5
      rmin(4) = 2.62d-6  !1.320e-5
c
c  Default values: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'ethane'
      gasname(2) = 'methane'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Set to zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1   !map ethane gas to ethane cloud
      igrowgas(5) = 2   !map methane gas to methane cloud
      igrowgas(8) = 2   !map methane gas to mixed cloud
      igrowgas(11) = 1  !map ethane gas to ethane growcore
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .true.
      is_grp_ice(3) = .true.
      is_grp_ice(4) = .true.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
c  Expand array to allow for nucleation of more than one gas to a group:
c  inucgas(i,igroup) 
c  Read as: ith gas to nucleate onto igroup is inucgas(i,igroup)
c
      inucgas(1,1) = 1   !Ethane onto tholin
      inucgas(2,1) = 2   !Methane onto tholin 
      inucgas(1,2) = 2   !Methane onto ethane clouds
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c  (<nnuc2elem(ielem)> will be set in {setupnuc})
c
      inuc2elem(1,1) = 3  !from <<tholin>> to <<tholin core in C2H6 cloud>>
      inuc2elem(2,1) = 6  !from <<tholin>> to <<tholin core in CH4 cloud>>
      inuc2elem(1,2) = 11 !from <<C2H6 cloud>> to <<C2H6 growcore>>  
      inuc2elem(1,3) = 9  !from <<th core in C2H6 cl>> to <<th core in CH4/C2H6 cl>>
      inuc2elem(1,4) = 10 !from <<cor2mom in C2H6 cl>> to <<cor2mom in CH4/C2H6 cl>> 
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c                                           (See top for more options)
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,3) = I_VAPORDEP
      inucproc(1,6) = I_DROPACT !!I_VAPORDEP
      inucproc(2,11) = I_VAPORDEP
      inucproc(3,9) = I_VAPORDEP
      inucproc(4,10) = I_VAPORDEP

      if( inucproc(1,6) .eq. I_DROPACT ) 
     $ write(*,*) 
     $ 'setupaer: Methane nucleates to tholin with droplet activation!'
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(3) = 1   !evap tholin core in ethane cloud to tholin
      ievp2elem(6) = 1   !evap tholin core in methane cloud to tholin
      ievp2elem(9) = 3   !evap tholin core in mixed cloud to tholin core in ethane cloud
      ievp2elem(10) = 4  !evap tholin core2m in mixed cloud to tholin core2m in ethane cloud
      ievp2elem(11) = 2  !evap ethane growcore to ethane cloud
c
c-------------------------------------------------------------------------------
c
c  Coagulation mapping
c
c  This <NGROUP> by <NGROUP> array maps aerosol groups for coagulation
c  <icoag(i,j)> is the particle group resulting from a collision
c  between particles in groups <i> and <j>.
c  [This array must be diagonal, so the user need only specify
c  elements (i=1:NGROUP,j=1:i). ]
c
        icoag(1,1) = 1  !tholin + tholin = tholin
        icoag(2,1) = 2  !C2H6 cloud + tholin = C2H6 cloud
        icoag(2,2) = 2  !C2H6 cloud + C2H6 cloud = C2H6 cloud
        icoag(3,1) = 3  !CH4 cloud + tholin = CH4 cloud
        icoag(3,2) = 4  !CH4 cloud + C2H6 cloud = CH4/C2H6 cloud
        icoag(3,3) = 3  !CH4 cloud + CH4 cloud = CH4 cloud
        icoag(4,1) = 4  !CH4/C2H6 cloud + tholin = CH4/C2H6 cloud
        icoag(4,2) = 4  !CH4/C2H6 cloud + C2H6 cloud = CH4/C2H6 cloud
        icoag(4,3) = 4  !CH4/C2H6 cloud + CH4 cloud = CH4/C2H6 cloud 
        icoag(4,4) = 4  !CH4/C2H6 cloud + CH4/C2H6 cloud = CH4/C2H6 cloud
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for purCH4 model.
c
      return
      end
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_methane_volcore
c
c  Methane nucleation only; ice and droplet clouds in same group
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 
     &    'Enter setupaer_methane_volcore'
      write(*,'(/,a)') '2 Groups - Combine ice/liq w/ volcores'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Name for each group
c
      groupname(1) = 'Aerosol group'
      groupname(2) = 'Methane Cloud group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 4
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'CH4 mixed phase'
      elemname(3) = 'CH4 ice core'
      elemname(4) = 'tholin core mass'
      elemname(5) = '2nd moment of tholin mass'
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice
c     I_CH4_ICE  is methane ice
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_CH4_LIQ  !_MIXPHASE
      icomp(3) = I_CH4_ICE
      icomp(4) = I_CxHyNz
      icomp(5) = I_CxHyNz
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
      itype(3) = I_VOLCORE
      itype(4) = I_COREMASS
      itype(5) = I_CORE2MOM
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3d-7
      rmin(2) = 2.08d-6  
c
c  Default values: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'methane'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Set to zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1   !map methane gas to methane cloud
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .false.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
c  Expand array to allow for nucleation of more than one gas to a group:
c  inucgas(i,igroup) 
c  Read as: ith gas to nucleate onto igroup is inucgas(i,igroup)
c
      inucgas(1,1) = 1   !Methane onto tholin
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c  (<nnuc2elem(ielem)> will be set in {setupnuc})
c
      inuc2elem(1,1) = 4  !from <<tholin>> to <<tholin core in CH4 cloud>>
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c                                           (See top for more options)
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,4) = I_VAPORDEP
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(4) = 1   !evap tholin core in methane cloud to tholin
c
c-------------------------------------------------------------------------------
c
c  Coagulation mapping
c
c  This <NGROUP> by <NGROUP> array maps aerosol groups for coagulation
c  <icoag(i,j)> is the particle group resulting from a collision
c  between particles in groups <i> and <j>.
c  [This array must be diagonal, so the user need only specify
c  elements (i=1:NGROUP,j=1:i). ]
c
        icoag(1,1) = 1  !tholin + tholin = tholin
        icoag(2,1) = 2  !CH4 cloud + tholin = CH4 cloud
        icoag(2,2) = 2  !CH4 cloud + CH4 cloud = CH4 cloud
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for purCH4 model.
c
      return
      end
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_mixed_volcore
c
c  Both ethane and methane nucleation: ethane clouds and methane clouds with
c  ethane growcores.  Ice/droplets in same group.
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 
     &    'Enter setupaer_mixed_volcore'
      write(*,'(/,a)') '3 Groups - Combine ice/liq w/ volcores'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Name for each group
c
      groupname(1) = 'Aerosol group'
      groupname(2) = 'Ethane Cloud group'
      groupname(3) = 'Mixed Cloud group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 4
      nelemg(3) = 6
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'C2H6 mixed phase'
      elemname(3) = 'C2H6 ice core'
      elemname(4) = 'tholin core mass'
      elemname(5) = '2nd moment of tholin mass'
      elemname(6) = 'CH4/C2H6 mixed cloud '
      elemname(7) = 'CH4 ice core'
      elemname(8) = 'C2H6 mixed phase growcore'
      elemname(9) = 'C2H6 ice core'
      elemname(10) = 'tholin core mass' 
      elemname(11) = '2nd moment of tholin mass'
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice
c     I_CH4_ICE  is methane ice
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_C2H6_LIQ !_MIXPHASE
      icomp(3) = I_C2H6_ICE
      icomp(4) = I_CxHyNz
      icomp(5) = I_CxHyNz
      icomp(6) = I_CH4_LIQ  !_MIXPHASE
      icomp(7) = I_CH4_ICE
      icomp(8) = I_C2H6_LIQ !_MIXPHASE
      icomp(9) = I_C2H6_ICE
      icomp(10) = I_CxHyNz
      icomp(11) = I_CxHyNz
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
      itype(3) = I_VOLCORE
      itype(4) = I_COREMASS
      itype(5) = I_CORE2MOM
      itype(6) = I_VOLATILE
      itype(7) = I_VOLCORE
      itype(8) = I_GROWCORE
      itype(9) = I_VOLCORE
      itype(10) = I_COREMASS
      itype(11) = I_CORE2MOM
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3d-7
      rmin(2) = 2.08d-6  
      rmin(3) = rmin(2) * 2.**(1./3.)
c
c  Default values: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'ethane'
      gasname(2) = 'methane'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Set to zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1   !map ethane gas to ethane cloud
      igrowgas(6) = 2   !map methane gas to mixed cloud
      igrowgas(8) = 1   !map ethane gas to ethane growcore
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .false.
      is_grp_ice(3) = .false.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
c  Expand array to allow for nucleation of more than one gas to a group:
c  inucgas(i,igroup) 
c  Read as: ith gas to nucleate onto igroup is inucgas(i,igroup)
c
      inucgas(1,1) = 1   !Ethane onto tholin
      inucgas(1,2) = 2   !Methane onto ethane clouds
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c  (<nnuc2elem(ielem)> will be set in {setupnuc})
c
      inuc2elem(1,1) = 4  !from <<tholin>> to <<tholin core in C2H6 cloud>>
      inuc2elem(1,2) = 8  !from <<C2H6 cloud>> to <<C2H6 growcore>>  
      inuc2elem(1,3) = 9  !from <<C2H6 ice core>> to <<C2H6 ice core in mix cl>>
      inuc2elem(1,4) = 10 !from <<th core in C2H6 cl>> to <<th core in CH4/C2H6 cl>>
      inuc2elem(1,5) = 11 !from <<cor2mom in C2H6 cl>> to <<cor2mom in CH4/C2H6 cl>> 
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c                                           (See top for more options)
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,4) = I_VAPORDEP
      inucproc(2,8) = I_VAPORDEP
      inucproc(3,9) = I_VAPORDEP
      inucproc(4,10) = I_VAPORDEP
      inucproc(5,11) = I_VAPORDEP
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(4) = 1   !evap tholin core in ethane cloud to tholin
      ievp2elem(8) = 2   !evap ethane growcore in mixed cloud to ethane cloud
      ievp2elem(9) = 3   !evap ethane ice core in mixed cloud to ethane ice core in ethane cloud
      ievp2elem(10) = 4  !evap tholin core in mixed cloud to tholin core in ethane cloud
      ievp2elem(11) = 5  !evap tholin core2m in mixed cloud to tholin core2m in ethane cloud
c
c-------------------------------------------------------------------------------
c
c  Coagulation mapping
c
c  This <NGROUP> by <NGROUP> array maps aerosol groups for coagulation
c  <icoag(i,j)> is the particle group resulting from a collision
c  between particles in groups <i> and <j>.
c  [This array must be diagonal, so the user need only specify
c  elements (i=1:NGROUP,j=1:i). ]
c
        icoag(1,1) = 1  !tholin + tholin = tholin
        icoag(2,1) = 2  !C2H6 cloud + tholin = C2H6 cloud
        icoag(2,2) = 2  !C2H6 cloud + C2H6 cloud = C2H6 cloud
        icoag(3,1) = 3  !CH4/C2H6 cloud + tholin = CH4/C2H6 cloud
        icoag(3,2) = 3  !CH4/C2H6 cloud + C2H6 cloud = CH4/C2H6 cloud
        icoag(3,3) = 3  !CH4/C2H6 cloud + CH4/C2H6 cloud = CH4/C2H6 cloud
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for purCH4 model.
c
      return
      end
c-------------------------------------------------------------------------------
       subroutine setupaer_methane_2phase
c
c
c  Methane nucleation, methane clouds can be ice or liquid
c  (liquid clouds include N2)
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer_methane_2phase'
      write(*,'(/,a)') 'Ice and liq CH4 clouds, no ethane'
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Name for each group
c
      groupname(1) = 'Aerosol group'
      groupname(2) = 'Methane Ice Cloud group'
      groupname(3) = 'Methane Liq Cloud group'
c
c  Number of elements in each group (elements in a group can include
c  particle number concentration, mass concentrations of cores, and second
c  moments of core mass distributions).
c
      nelemg(1) = 1
      nelemg(2) = 3
      nelemg(3) = 3
c
c  Name for each element
c
      elemname(1) = 'tholin CN'
      elemname(2) = 'CH4 ice crystal'
      elemname(3) = 'tholin core mass' 
      elemname(4) = '2nd moment of tholin mass in ice'
      elemname(5) = 'CH4 droplet'
      elemname(6) = 'tholin core mass' 
      elemname(7) = '2nd moment of tholin mass in drop'
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_C2H6_ICE is ethane ice
c     I_CH4_ICE  is methane ice
c     I_CH4_LIQ  is liquid methane 
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_CH4_ICE
      icomp(3) = I_CxHyNz
      icomp(4) = I_CxHyNz
      icomp(5) = I_CH4_LIQ
      icomp(6) = I_CxHyNz
      icomp(7) = I_CxHyNz
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
      itype(4) = I_CORE2MOM
      itype(5) = I_VOLATILE
      itype(6) = I_COREMASS
      itype(7) = I_CORE2MOM
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3e-7
      rmin(2) = 2.62d-6  !1.320e-5
      rmin(3) = 2.62d-6  !1.320e-5
c
c  Defaults: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
c-------------------------------------------------------------------------------
c
c  Define mapping arrays and parameters for condensational growth, evaporation,
c  and nucleation.
c
c==Set up gas descriptions and mapping arrays used for condensational growth.
c
c  Names of gas species
c
      gasname(1) = 'methane'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Leave this zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1
      igrowgas(5) = 1
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .true.
      is_grp_ice(3) = .false.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
      inucgas(1,1) = 1
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c  (<nnuc2elem(ielem)> will be set in {setupnuc})
c
      inuc2elem(1,1) = 3  !tholin -> core in ch4 ice crystal
      inuc2elem(2,1) = 6  !tholin -> core in ch4 droplet
       !include melting
      inuc2elem(1,2) = 5
      inuc2elem(1,3) = 6
      inuc2elem(1,4) = 7
       !include freezing
      inuc2elem(1,5) = 2
      inuc2elem(1,6) = 3
      inuc2elem(1,7) = 4
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c                                           (See top for more options)
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,3) = I_VAPORDEP
      inucproc(1,6) = I_VAPORDEP
     !melting
      inucproc(2,5) = I_ICEMELT
      inucproc(3,6) = I_ICEMELT
      inucproc(4,7) = I_ICEMELT
     !freezing
      inucproc(5,2) = I_DROPFREEZE
      inucproc(6,3) = I_DROPFREEZE
      inucproc(7,4) = I_DROPFREEZE
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(3) = 1
      ievp2elem(6) = 1
c
c-------------------------------------------------------------------------------
c
c  Coagulation mapping:
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
        icoag(3,1) = 3
        icoag(3,2) = 0  !!no mixing liq/ice mixed clouds
        icoag(3,3) = 3
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for ethane clouds, methane 
c  clouds with ethane growcores - where methane clouds can be ice or liquid.
c
      return
      end
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c-------------------------------------------------------------------------------
       subroutine setupaer_water_ice
c                                                                              
c  Water clouds                                                              
c                                                                           
c  Include global constants and variables                                  
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupaer_water_ice'
      write(*,'(/,a)') '2 Groups - Water Clouds'
c
c-------------------------------------------------------------------------------
c
c==Set up particle types and mapping arrays and structure of size grid.
c
c  Name for each group
c
      groupname(1) = 'Haze Group'
      groupname(2) = 'Water Cloud Group'
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
      elemname(2) = 'H2O ice crystal'
      elemname(3) = 'tholin core mass'
      elemname(4) = '2nd moment of tholin mass in ice'
c
c
c  This array specifies the composition of each element:
c     I_CxHyNz   is tholin                  (See top for more options)
c     I_ICE      is ice water      
c
      icomp(1) = I_CxHyNz
      icomp(2) = I_ICE
      icomp(3) = I_CxHyNz
      icomp(4) = I_CxHyNz
c
c
c  This array specifies the type of each element:
c     I_INVOLATILE is CN (involatile) number concentration [#/cm^3]
c     I_VOLATILE   is water droplet (volatile) number concentration [#/cm^3]
c     I_COREMASS   is core mass concentration [g/cm^3]
c     I_CORE2MOM   is core second moment (mass^2) [g^2/cm^3]
c
      itype(1) = I_INVOLATILE
      itype(2) = I_VOLATILE
      itype(3) = I_COREMASS
      itype(4) = I_CORE2MOM
c
c  Minimum radius for each group (used below to calculate <rmassmin>) [cm]
c
      rmin(1) = 1.3e-7
      rmin(2) = rmin(1)*rmrat(1)**(1./3.)
c
c  Defaults: <rmrat> = 2.0   <ishape> = 1   <eshape> = 1
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
      gasname(1) = 'water'
c
c  Array <igrowgas> maps a particle element to its associated gas for
c  condensational growth/evaporation.
c
c  Leave this zero if there is no growth specific to particle element
c  (i.e., use zero if cores do not grow, even if they are a component
c  of a particle group that does grow; use zero for core second moment).
c
      igrowgas(2) = 1
c
c
c  If <is_grp_ice> = .true. then the particle group is an ice crystal,
c  else the particle group is liquid (or does not grow).  This array
c  is used to select the appropriate ventilation factors in setupgkern.f.
c
      is_grp_ice(1) = .false.
      is_grp_ice(2) = .true.
c
c==Set up mapping arrays for nucleation and total evaporation.
c
c  Array <inucgas> maps a particle group to its associated gas for nucleation:
c  Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
c  Leave this zero if particles are not subject to nucleation.
c
      inucgas(1,1) = 1
c
c
c  Nucleation mapping:
c
c  Nucleation transfers particle mass from element <ielem> to element
c  <inuc2elem(i,ielem)>, where <i> ranges from 0 to <nnuc2elem(ielem)>.
c
      inuc2elem(1,1) = 3
c
c
c  <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
c  particles from element <ielem> to element <ieto>:
c                                           (See top for more options)
c     I_DROPFREEZE: Droplet homogeneous freezing
c     I_VAPORDEP: Heterogeneous nucl by vapor deposition
c  Leave this zero if particles are not subject to nucleation.
c
      inucproc(1,3) = I_DROPFREEZE
c
c
c  Total evaporation mapping: total evaporation transfers particle mass from
c  element <ielem> to element <ievp2elem(ielem)>.
c  Leave this zero if element is not subject to total evaporation.
c
c  This array is not automatically derived from <inuc2elem> because multiple
c  elements can nucleate to a particular element (reverse mapping is not
c  unique).
c
      ievp2elem(3) = 1
c
c
c  Coagulation mapping:
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
c
c  Return to caller with aerosol and cloud microphysics mapping arrays
c  and time-independent parameters defined for the water cloud case.
c
      return
      end
c
c
