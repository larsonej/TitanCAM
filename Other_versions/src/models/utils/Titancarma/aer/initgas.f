       subroutine initgas
c
c
c  @(#) initgas.f  Ackerman  Dec-1995
c  This routine initializes the atmospheric profiles of all gases.
c
c    gc       Gas concentration at layer mid-point [g/cm^3]
c
c  Presently the only vertical coordinate is altitude, zl [cm].
c
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
c  Local declarations
c
      logical do_huygens, do_lellouch, do_fixed_CH4
      logical do_kostiuk, do_fixed_C2H6

c     <qLH> are mixing ratio values for methane
c             from Lellouch et al. Icarus, 79, 328-349 (1989)
c             (values are for 2 km grid spacing)
c
      dimension qLH(51) !NZP1)
c
c     <qHP> are mixing ratio values for methane
c             from Huygens probe data (Niemann et al. 2005,
c             Nature vol. 438, p. 779-784, Fig. 2)

      dimension qHP(15), aHP(15)

c
c  Methane mixing ratio data..................................
      data qHP / 0.0492, 0.0492, 0.0475, 0.0380, 0.0320, 
     $             0.0260, 0.0220, 0.0200, 0.0169, 0.0141, 
     $             0.0161, 0.0162, 0.0150, 0.0144, 0.0141 /
      data aHP / 0.0, 5.0, 8.0, 10.0, 12.0, 16.0, 20.0, 23.0, 
     $            27.0, 33.0, 37.0, 42.0, 85.0, 114.0, 135.0 /

      data qLH / 0.080, 0.080, 0.080, 0.065, 0.055, 0.045,
     $      0.040, 0.034, 0.029, 0.026, 0.023, 0.021, 0.020,
     $      0.018, 0.017, 0.017, 0.016, 0.015, 0.015, 0.015,
     $      0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,
     $      0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,
     $      0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,
     $      0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,
     $      0.015, 0.015, 0.015 /
c.............................................................
c
c   Define formats
c
    1 format(/,'Gas concentrations for ',a,'(initgas)',//,
     $  'Triple point T = ',0p,f6.3,//,a3, 1x, 7(a11,4x), /) 
    2 format(i3,1x,1p,7(e11.3,4x),0p,f11.3)
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter initgas'
c
c-------------------------------------------------------------------------------
c
c***Options for initial profiles************************************
c  Methane
      do_huygens = .true.     ! Use Huygens data
      do_lellouch = .false.   ! Use Lellouch mixing ratio values
      do_fixed_CH4 = .false.  ! Construct methane profile from fixed
                              !  surface humidity <RHs>, which is constant 
        RHs_ch4 = 50.         !  until follows saturation curve and then
                              !  constant again at minimum saturation value

c  Ethane
      do_kostiuk = .true.     ! Construct ethane profile using Kostiuk et al. 
                              !  stratospheric measurement for const mixing ratio
                              !  or follow saturation curve where profile would be
                              !  supersaturated
      do_fixed_C2H6 = .false. ! Fixed surface humidity (see methane explanation) 
        RHs_c2h6 = 20.
c
c********************************************************************
c
c   Calculate vapor pressures
      do ixyz = 1,NXYZ
        call vaporp
      enddo
 
c
c   Parameters for gas profile
c
c   Loop over all spatial dimensions and gases
c
      do igas = 1,NGAS

c   Gas constant for each gas
       rvap = RGAS/gwtmol(igas)

c   Report initialization information, set up alt-independent q
        if( gasname(igas) .eq. 'methane' ) then
         if( do_huygens) then
          write(*,*) 'Initialize methane with huygens probe data'
         elseif( do_lellouch ) then
          write(*,*) 'Initialize methane using Lellouch profile'
         elseif( do_fixed_CH4 ) then
          write(*,*) 'Construct methane profile from 
     $                           fixed surface humidity',RHs_ch4
          q_CH4 = RHs_ch4/100. * pvapl3(1,igas)/ (rhoa3(1)*rvap*t3(1))
         endif
        elseif( gasname(igas) .eq. 'ethane' ) then
         if( do_kostiuk) then
          write(*,*) 'Construct ethane profile from stratospheric
     $                    measurement (Kostiuk et al. 1997)'
         elseif( do_fixed_C2H6 ) then
          write(*,*) 'Construct ethane profile from 
     $                           fixed surface humidity',RHs_c2h6
          q_C2H6 = RHs_c2h6/100. * pvapl3(1,igas)/ (rhoa3(1)*rvap*t3(1))
         endif
        endif

c
       do ix = 1,NX
        do iy = 1,NY
         do iz = 1,NZ

          if( gasname(igas) .eq. 'methane' )then !---------------------

           if( do_huygens ) then
            !Find surrounding z values from input altitude grid
             k = 1
             zc_km = zc3(iz) / 1.d5
             do while( aHP(k) .lt. zc_km )
               k = k + 1
             enddo !while
             q_CH4 = ( qHP(k)-qHP(k-1) ) /
     $                  ( aHP(k) - aHP(k-1) )       
     $                * (zc_km - aHP(k-1)) + qHP(k-1)

           elseif( do_lellouch ) then
             if( dz3(iz) .ne. 2.d5 ) then
               write(*,*) 'Wrong grid for Lellouch profile',dz3(iz)
               stop
             endif
             q_CH4 = qLH(iz)

           endif
c            if( do_combine ) then
c              q_CH4 = 0.5*( q_CH4 + qLH(iz) )  !average
c            endif

             gc(ix,iy,iz,igas) = q_CH4 * rhoa(ix,iy,iz)
     $                            * gwtmol(igas)/WTMOL_AIR

             pp = gc(ix,iy,iz,igas) * Rvap * t(ix,iy,iz)
             if( pp .ge. pvapl3(iz,igas) ) then
                write(*,*) 'Methane supersaturated: ',iz,
     $            gc(ix,iy,iz,igas),pp/pvapl3(iz,igas)

                gc(ix,iy,iz,igas) = pvapl3(iz,igas) 
     $                                  /rvap /t(ix,iy,iz) * 0.99
             endif

             if( do_fixed_CH4 .and. iz.gt.1) then 
                !Follow constant mixing ratio when minimum value is reached
               rmix = gc(ix,iy,iz,igas)/rhoa(ix,iy,iz)
               rmixd = gc(ix,iy,iz-1,igas)/rhoa(ix,iy,iz-1)
               if( rmix .gt. rmixd ) then
               q_CH4 = rmixd
               gc(ix,iy,iz,igas) = q_CH4 * rhoa(ix,iy,iz)
     $                              * gwtmol(igas)/WTMOL_AIR
               endif
             endif
 
          else if( gasname(igas) .eq. 'ethane' )then !--------------------

            if( do_kostiuk ) q_C2H6 = 1.0d-5

            gc(ix,iy,iz,igas) = q_C2H6 * rhoa(ix,iy,iz) 
     $                                       * gwtmol(igas)/WTMOL_AIR

            pp = gc(ix,iy,iz,igas) * Rvap * t(ix,iy,iz)
            if( pp .ge. pvapi3(iz,igas) )
     $               gc(ix,iy,iz,igas) = pvapi3(iz,igas) 
     $                                      /rvap /t(ix,iy,iz) * 0.99

             if( do_fixed_C2H6 .and. iz.gt.1) then 
                !Follow constant mixing ratio when minimum value is reached
               rmix = gc(ix,iy,iz,igas)/rhoa(ix,iy,iz)
               rmixd = gc(ix,iy,iz-1,igas)/rhoa(ix,iy,iz-1)
               if( rmix .gt. rmixd ) then
               q_C2H6 = rmixd
               gc(ix,iy,iz,igas) = q_C2H6 * rhoa(ix,iy,iz)
     $                              * gwtmol(igas)/WTMOL_AIR
               endif
             endif

          else
            write(LUNOPRT,'(/,a)') 'invalid <igas> in initgas.f'
            stop 1
          endif

         enddo
        enddo
       enddo
      enddo
c
c
c  Specify fluxes at top and bottom of model [g/cm^2/s]
c  Define triple-point temperature (K) for pure substance
c
c    Ethane top flux from Yung et al. Ap. J. Suppl. 55, 465-506 (1984)
c    Table 5B

      do igas = 1,NGAS
       !Initialize surface reservoir for liquid methane and ethane
        puddle(igas) = ZERO
        do ixy = 1,NXY
 
          if( gasname(igas) .eq. 'methane' ) then  
 
            ftopgas(ixy,igas) = 0.
 
           !T0(igas) = 90.69 
            T0(igas) = 80.6 !!90.69 
              print *, 'Set methane triple point to mixed value!'
 
          else if( gasname(igas) .eq. 'ethane' ) then  
 
            ftopgas(ixy,igas) = 5.8d9/AVG * gwtmol(igas)  

            ftopgas(ixy,igas) = 1.7d9/AVG * gwtmol(igas)  !Eric's value from AGU2007
            if(ftopgas(ixy,igas) .lt. 5d9) 
     $         write(*,*) 'Ethane using Wilson value'

            T0(igas) = 90.348
          endif
          fbotgas(ixy,igas) = 0.
        enddo
        Tfreez(igas) = T0(igas)  ! Set freezing point temperature
      enddo                      !  to triple point
 
c
c  Scale particle concentrations and boundary fluxes from 
c  cartesian coordinates to coordinate system specified by <igrid>
c
c  Pick indices for top and bottom layers
c  
      if( igridv .eq. I_CART )then
        iztop = NZ
        izbot = 1
      else
        iztop = 1
        izbot = NZ
      endif

      do igas = 1, NGAS

        do ixyz = 1,NXYZ
          gc3(ixyz,igas) = gc3(ixyz,igas) * 
     $                     ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )
        enddo

        do ixy = 1,NXY
          ftopgas(ixy,igas) = ftopgas(ixy,igas) *
     $                        ( xmet2(ixy,iztop)*ymet2(ixy,iztop) )
          fbotgas(ixy,igas) = fbotgas(ixy,igas) *
     $                        ( xmet2(ixy,izbot)*ymet2(ixy,izbot) )
        enddo
      enddo
c
c
c  Initialize <supsati> and <supsatl>
c
      do ixyz = 1,NXYZ
        call supersat
      enddo
c
c
c  Print gas concentrations at 1 horizontal grid point
c
      ix = 1
      iy = 1
      ixy = NX * ( iy - 1 ) + ix 
c
c
c  Set vertical loop index to increment downwards
c
      if( igridv .eq. I_CART )then
        kb  = NZ
        ke  = 1
        idk = -1
      else 
        kb  = 1
        ke  = NZ
        idk = 1
      endif

      do igas = 1,NGAS

        write(LUNOPRT,1) gasname(igas),T0(igas),
     $         'iz','zc','gc [g/cm^3]','pvapl','pvapi',
     $                         'supsatl','supsati','T [K]'

        do iz = kb,ke,idk
          xyzmet = xmet2(ixy,iz)*ymet2(ixy,iz)*zmet2(ixy,iz)
          write(LUNOPRT,2) iz, zc2(ixy,iz), gc2(ixy,iz,igas)/xyzmet,
     $       pvapl2(ixy,iz,igas), pvapi2(ixy,iz,igas),
     $       supsatl2(ixy,iz,igas), supsati2(ixy,iz,igas), t2(ixy,iz)
        enddo

      enddo
c
c
c  Specify the values of <gc> assumed just above(below) the top(bottom)
c  of the model domain.
c
      do igas = 1,NGAS
        do ixy=1,NXY
          gc_topbnd(ixy,igas) = gc2(ixy,NZ,igas)
          gc_botbnd(ixy,igas) = gc2(ixy,1,igas)
        enddo
      enddo
c
c  Return to caller with gas concentrations initialized.
c
      return
      end
