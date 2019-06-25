       subroutine initatm
c
c
c  @(#) initatm.f  McKie  Jul-1997
c  This is the overall grid-related initialization routine.  It is
c  split into subprograms, one for each type of supported grid coord system.
c 
c  A user setting this routine up for a particular modelling problem
c  would typically choose a coordinate system by defining the
c  <igridv> and <igridh> variables to one of the supported symbols
c  (e.g. I_CART, I_SIG for <igridv>, I_LL, I_LC, I_PS, I_ME for <igridh>)
c  in this routine, then do the specific grid setup and atmospheric
c  variables setup for that coordinate system in the appropriate
c  subroutine with name of the form g_*_*.
c
c  This routine initializes various geometrical & atmospheric profiles.
c  The following variables are defined by this routine:
c
c   Variables defined at vertical layer boundaries:
c
c    Geometry:
c      zl       vertical coord
c      zmet     vertical metric scale factor [d(zl)/dz]
c
c    Physical variables:
c      w        vertical velocity
c      u        east-west velocity
c      v        north-south velocity
c      dkz      vertical diffusion coefficient
c
c
c   Variables defined at vertical layer mid-points:
c
c    Geometry:
c      zc       vertical coord
c      xc       x position at center of grid box
c      xu       x position at upper (right) edge of grid box
c      xl       x position at lower (left) edge of grid box
c      dx       x grid-box width
c      xmet     x direction metric scale factor, d(x_distance)/dx 
c      yc       y position at center of grid box
c      yu       y position at upper (right) edge of grid box
c      yl       y position at lower (left) edge of grid box
c      dy       y grid-box width
c      ymet     y direction metric scale factor, d(y_distance)/dy 
c
c    Physical variables:
c      p        Air pressure [dyne/cm^2]
c      rhoa     Scaled air density [g/x_units/y_units/z_units]
c      t        Air temperature [K]
c      ptc      Potential temperature concentration [K g/x_units/y_units/z_units]
c      rmu      Dynamic air viscosity [g/cm/s]
c      thcond   Thermal conductivity of dry air [erg/cm/s/degree_K]
c      dkx      east-west diffusion coefficient
c      dky      north-south diffusion coefficient
c
c  Note on the way the metric scale factors <xmet>, <ymet>, <zmet> are defined:
c
c   These scale factors (the amount of actual measurable Cartesian distance
c   per unit change of the generalized coordinate in the direction of that
c   generalized coordinate) are chosen such that:
c
c     <xmet> * <ymet> * <zmet>  .gt.  0.
c     <xmet> * <ymet>  .gt.  0.
c
c   so that any volume (or area concentrations) are carried as positive
c   quantities in the model.
c
c   Since <xmet> and <ymet> are each naturally positive definite for
c   typical coordinates, and since the mathematical definition of <zmet> could
c   be positive or negative depending on whether the vertical coordinate
c   is positive upward (e.g. cartesian altitude) or positive downward 
c   (e.g. pressure or sigma), then in this model <zmet> is chosen to be
c   the absolute value of its mathematical definition (i.e., always positive).
c   The scaling of vertical wind <w> also follows this convention:
c   in all coordinates <w> is positive for updrafts. Hence, to scale vertical
c   wind in pressure or sigma units, multiply dp/dt or d(sigma)/dt by -1/<zmet>.
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
c
c   Define formats
c
    1 format(i3,1p,6(2x,e11.3))
    2 format(/,'initatm initializations:',/)
    3 format('Error--(initatm) Unknown grid combination: igridv = ',i5,
     $       ', igridh = ',i5)
    4 format(/,a,' at ix,iy=',2(1x,i4))
    5 format(/,a3,6(2x,a11))
    6 format(/,'Coordinate grid type: ',a)
    7 format(a,' = ',1p,e11.4)
    8 format(a,' = ',i5)
    9 format('Warning--(init): ',
     $       'igridv set to I_CART because do_parcel = .true')
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter initatm'
c
c
c  Define pointer to the type of grid coordinate system being used:
c  igridv for vertical coordinate, igridh for horizontal coordinates.
c
c  Possible values for igridv:
c       I_CART    cartesian
c       I_SIG     sigma
c       I_HYB     hybrid
c
c  Possible values for igridh:
c       I_CART   cartesian
c       I_LL     longitude_latitude
c       I_LC     lambert_conformal
c       I_PS     polar_stereographic
c       I_ME     mercator
c
c      igridv = I_CART
      igridv = I_HYB
c     igridv = I_SIG
c
c     igridh = I_CART
      igridh = I_LL
c     igridh = I_LC
c     igridh = I_PS
c     igridh = I_ME
c
c
c  Parcel model in cartesian coordinates only
c  (pressure coordinates would also be fine, but not yet set up)
c  
      if( do_parcel .and. ((igridv .eq. I_SIG) .or. 
     $                     (igridv .eq. I_HYB)))then
        igridv = I_CART
        write(LUNOPRT,9)
      endif

!      print*,igridv, igridh
!      print*,'coordinates',I_CART, I_SIG,I_HYB,I_LL,I_LC,I_PS,I_ME
c
c
c  Do appropriate grid-related initializations depending on grid type
c
!      print*,'testinitatm'
      if( igridh .eq. I_CART .and. igridv .eq. I_CART )then

        call g_cart_cart

      else if( igridh .eq. I_CART .and. igridv .eq. I_SIG )then

        call g_cart_sig

      else if( igridh .eq. I_LL .and. igridv .eq. I_SIG )then

        call g_ll_sig

      else if( igridh .eq. I_LC .and. igridv .eq. I_SIG )then

        call g_lc_sig

      else if( igridh .eq. I_PS .and. igridv .eq. I_SIG )then

        call g_ps_sig

      else if( igridh .eq. I_ME .and. igridv .eq. I_SIG )then

        call g_me_sig
      
      else if( igridh .eq. I_LL .and. igridv .eq. I_HYB )then

        call g_ctm_ll_hyb

      else
        write(LUNOPRT,3) igridv, igridh
        call endcarma
      endif

!      print*, 'testinitatm'
c
c
c  Scale fields defined at midpoint of vertical layers using metric factors.
c   Rule of thumb:  Any quantity with unit of distance factor ds{x,y,z} should
c   have a factor of 1./{x,y,z}met.  
c
      do ixyz=1,NXYZ
        u3(ixyz) = u3(ixyz) / xmet3(ixyz)
        v3(ixyz) = v3(ixyz) / ymet3(ixyz)
        dkx3(ixyz) = dkx3(ixyz) / ( xmet3(ixyz)**2 )
        dky3(ixyz) = dky3(ixyz) / ( ymet3(ixyz)**2 )
        xyzmet = xmet3(ixyz) * ymet3(ixyz) * zmet3(ixyz)
        ptc3(ixyz) = ptc3(ixyz) * xyzmet
        rhoa3(ixyz) = rhoa3(ixyz) * xyzmet
      enddo
c
c
c  Specify the values of <ptc> assumed just above(below) the top(bottom)
c  of the model domain.
c
      do ixy=1,NXY
        ptc_topbnd(ixy) = ptc2(ixy,NZ)
        ptc_botbnd(ixy) = ptc2(ixy,1)
      enddo
c
c
c  Hydrostatically balance initial atmospheric profiles 
c  (consistent with subsequent hydrostat calls)
c
c  NOTE: For I_HYB, use the atmosphere provided by the ctm, which should
c  already be in hydrostatic balance.
      if( ( .not. do_parcel ) .and. ( igridv .ne. I_HYB ) )then
        call hydrostat
      endif
c
c
c  Define appropriate sign factor for <w>.  Non-Cartesian vertical coordinates
c  are assumed to be positive downward, but <w> in this model is always assumed
c  to be positive upward. 
c
      if( igridv .ne. I_CART ) then
        wsign = -1.
      else
        wsign = 1.
      endif
c
c
c  Scale fields defined at boundaries of vertical layers using metric factors.
c   Linearly interpolate/extrapolate nearest 2 vertical midpoint metric factors.
c   Set boundary vertical velocities and diffusion coefficients to zero 
c   (since vertical transport ignores them)
c
      do k=1,NZP1

        if( k .eq. 1 )then
         k1 = 1
         k2 = min( NZ, 2 )
        else if( k .eq. NZP1 )then
         k1 = NZ
         k2 = NZ
        else
         k1 = min( NZ, max( 1, k - 1 ) )
         k2 = k
        endif

        do ixy=1,NXY

          a_mid_k1 = zc2(ixy,k1) 
          a_mid_k2 = zc2(ixy,k2) 

          if(  a_mid_k2 .ne. a_mid_k1 )then
            frac = ( zl2(ixy,k) - a_mid_k1 ) / ( a_mid_k2 - a_mid_k1 ) 
          else
            frac = 0.
          endif

          xmet_k1 = xmet2(ixy,k1)
          xmet_k2 = xmet2(ixy,k2)
          xmet_k = xmet_k1 + frac * ( xmet_k2 - xmet_k1 )
          dkx2(ixy,k) = dkx2(ixy,k) / ( xmet_k**2 )

          ymet_k1 = ymet2(ixy,k1)
          ymet_k2 = ymet2(ixy,k2)
          ymet_k = ymet_k1 + frac * ( ymet_k2 - ymet_k1 )
          dky2(ixy,k) = dky2(ixy,k) / ( ymet_k**2 )

          if( k .eq. 1 .and. ibbnd_pc .eq. I_FLUX_SPEC ) then
            dkz2(ixy,k) = 0.
            w2(ixy,k) = 0.
          endif

          if( k .eq. NZP1 .and. itbnd_pc .eq. I_FLUX_SPEC ) then
            dkz2(ixy,k) = 0.
            w2(ixy,k) = 0.
          endif

          zmet_k1 = zmet2(ixy,k1)
          zmet_k2 = zmet2(ixy,k2)
          zmet_k = zmet_k1 + frac * ( zmet_k2 - zmet_k1 )
          dkz2(ixy,k) = dkz2(ixy,k) / ( zmet_k**2 )

          w2(ixy,k) = w2(ixy,k) / zmet_k
          w2(ixy,k) = wsign * w2(ixy,k)

        enddo
      enddo
c
c
c  Announce this routine's results
c
      if (do_print_setup) then
        write(LUNOPRT,2)
c
c
c  Report type of coordinate grid
c
        write(LUNOPRT,8) 'igridv ', igridv
        write(LUNOPRT,8) 'igridh ', igridh
        write(LUNOPRT,6) gridname
c
c
c  Report coord system/grid/projection selection variables
c
        write(LUNOPRT,7) 'dom_llx', dom_llx
        write(LUNOPRT,7) 'dom_lly', dom_lly
        write(LUNOPRT,7) 'dom_urx', dom_urx
        write(LUNOPRT,7) 'dom_ury', dom_ury
        write(LUNOPRT,7) 'rlon0  ', rlon0
        write(LUNOPRT,7) 'rlat0  ', rlat0
        write(LUNOPRT,7) 'rlat1  ', rlat1
        write(LUNOPRT,7) 'rlat2  ', rlat2
        write(LUNOPRT,7) 'hemisph', hemisph
      endif
c
c
c  Define indices of a horiz grid pt at which to print atm vert structure
c
      ix = 1
      iy = 1
c
c
c  Set vertical loop index to increment downwards
c
c      if( igridv .eq. I_CART )then
c        kb  = NZ
c        ke  = 1
c        idk = -1
c      else if(( igridv .eq. I_SIG ) .or. ( igridv .eq. I_HYB ))then
c        kb  = 1
c        ke  = NZ
c        idk = 1
c      else
c        write(LUNOPRT,8) 'bad igridv = ',igridv
c        call endcarma
c      endif
c
c
c  Print atmospheric structure at horizontal grid point (ix,iy)
c
c      write(LUNOPRT,4) 'Sample atm structure', ix, iy
c      write(LUNOPRT,5) 'k', 'zc', 'p', 'rhoa', 't', 'rmu', 'thcond'
c      do k = kb,ke,idk
c        xyzmet = xmet(ix,iy,k)*ymet(ix,iy,k)*zmet(ix,iy,k)
c        write(LUNOPRT,1)
c     $    k,zc(ix,iy,k),p(ix,iy,k),rhoa(ix,iy,k)/xyzmet,
c     $    t(ix,iy,k),rmu(k),thcond(k)
c      enddo
c
c
c  Print detailed vertical structure at horizontal grid point (ix,iy)
c
c      write(LUNOPRT,4) 'Vertical grid structure at ',ix, iy
c      write(LUNOPRT,5) 'k', 'zc(k)', 'zl(k)','zl(k+1)','w(k)'
c      do k=kb,ke,idk
c       write(LUNOPRT,1)
c     $   k, zc(ix,iy,k),
c     $   zl(ix,iy,k), zl(ix,iy,k+1), w(ix,iy,k)
c      enddo 
c
c
c  Return to caller with atmospheric profiles initialized.
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine g_cart_cart
c
c
c  @(#) g_cart_cart.f  Ackerman  Oct-1995
c  This routine handles atmospheric grid initialization
c  for the horizontal=cartesian and vertical=cartesian grid coordinates.
c
c  For this coord system and grid:
c   Horiz surface is a cartesian x vs y plane,
c   and vertical coordinate is cartesian altitude.
c
c   x,y,z are regular measurable distance in cm.
c   x is positive to the right (east), y is positive north. z is positive up.
c   xmet, ymet metric factors are all unity with no units.
c
c  Modified slightly for use with generalized coordinates.  Jul-1997, McKie.
c  (From v1.11 single geometrical coord system model's initatm.f routine.)
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter g_cart_cart'
c
c
c  Define bogus values for grid selection parameters that have no meaning
c  for this coordinate system.
c
      rlon0 = -999.
      rlat0 = -999.
      rlat1 = -999.
      rlat2 = -999.
      hemisph = -999.
c
c
c  Define descriptive text for type of grid
c
      gridname = 'cartesian horizontal, cartesian vertical'
c
c
c For this example, grid spacing is uniform.
c  Define its east-west, north-south fixed grid box size.
c
      dxfix = 300.e2
      dyfix = 300.e2
c
c
c  Define lower left (southwest) upper right (northest) horizontal
c  domain limits for the grid
c
      dom_llx = 0.
      dom_lly = 0.
      dom_urx = NX * dxfix
      dom_ury = NY * dyfix
c
c
c  Specify horizontal grid coordinates and grid box sizes.
c   <dx> and <dy> are the east-west and north-south grid box thicknesses.
c   <xc>, <yc> are the coord at the center of each grid box.
c   <xl> & <xu> are the lower & upper grid box x coord at y=<yc>.
c   <yl> & <yu> are the lower & upper grid box y coord at x=<xc>.
c   <rlon> & <rlat> are longitude & latitude [deg] at planet surface
c     ( <rlon>, <rlat> usually assumed constant in horiz cartesian coord )
c
      do iy = 1,NY
	      do ix = 1,NX
          rlon(ix,iy) = -90.
          rlat(ix,iy) = 45.
	        do k = 1,NZ
	          dx(ix,iy,k) = dxfix
            xc(ix,iy,k) = dom_llx + ( ix - .5 ) * dx(ix,iy,k)
            xl(ix,iy,k) = dom_llx + ( ix - 1 ) * dx(ix,iy,k)
            xu(ix,iy,k) = dom_llx + ix * dx(ix,iy,k)
	          dy(ix,iy,k) = dyfix
            yc(ix,iy,k) = dom_lly + ( iy - .5 ) * dy(ix,iy,k)
            yl(ix,iy,k) = dom_lly + ( iy - 1 ) * dy(ix,iy,k)
            yu(ix,iy,k) = dom_lly + iy * dy(ix,iy,k)
          enddo
        enddo
      enddo
c
c
c  Define vertical profile of atmospheric variables that are 3-D.
c   In this demo:
c    The vertical profiles do not vary horizontally.
c    Air temp <t> lapse rate is dry adiabatic from a surface value of <t_sfc>,
c    resulting in a constant potential temperature.
c    Pressure <p> is integrated upwards hydrostatically.
c    Air density <rhoa> is a function of <t> & <p>, from dry air equa of state.
c    All units are cgs, deg_K.
c
!      t_sfc = 273.16 + 2.
!      t_sfc = 273.16 + 25. - 6.5625
!      t_sfc = 288.0
       t_sfc = 93.0
c
c
c  Define dry lapse rate 
c
      dlapse = GRAV / CP
!      dlapse = 6.5625e-5
c
c
c  Define (fixed) layer thickness
c
      dzfix = 10.e5
c
c
c  Define surface pressure [dyne/cm^2] 
c
      do iy=1,NY
        do ix=1,NX
          p_surf(ix,iy) = 1496 * RMB2CGS
        enddo
      enddo
c
c
c  Visit each horiz grid point & calculate atmospheric properties.
c
      do ix = 1,NX
       do iy = 1,NY
        t_surf(ix,iy) = t_sfc
c
c
c  Bottom layer is special because it cannot be evaluated from a lower layer.
c
        zbot = 0.
        k = 1
        zl(ix,iy,k) = zbot

        dz(ix,iy,k) = dzfix
        zl(ix,iy,k+1) = zl(ix,iy,k) + dzfix
        zc(ix,iy,k) = zl(ix,iy,k) + dzfix/2.
  
        t(ix,iy,k) = t_sfc - dlapse * ( dzfix / 2. )
        p(ix,iy,k) = p_surf(ix,iy) *
     $               exp( -dzfix/2. * GRAV/(R_AIR*t_sfc) )
        rhoa(ix,iy,k) = p(ix,iy,k)/(R_AIR*t(ix,iy,k))
        ptc(ix,iy,k) = rhoa(ix,iy,k) * 
     $                 t(ix,iy,k) * ( PREF/p(ix,iy,k) )**RKAPPA
c
c
c  Integrate pressure upwards hydrostatically.
c
        do k = 2,NZ

          dz(ix,iy,k) = dzfix
          zl(ix,iy,k+1) = zl(ix,iy,k) + dzfix
          zc(ix,iy,k) = zc(ix,iy,k-1) + dzfix

          t(ix,iy,k) = t(ix,iy,k-1) - dlapse*dzfix
          tbot = 0.5 * ( t(ix,iy,k-1) + t(ix,iy,k) )
          p(ix,iy,k) = p(ix,iy,k-1) *
     $                 exp( -dzfix*GRAV/(R_AIR*tbot) )
          rhoa(ix,iy,k) = p(ix,iy,k)/(R_AIR*t(ix,iy,k))
          ptc(ix,iy,k) = rhoa(ix,iy,k) * 
     $                   t(ix,iy,k) * ( PREF/p(ix,iy,k) )**RKAPPA

        enddo
c
c
c  Pressure at top of model domain [dyne/cm^2]
c
        p_top(ix,iy) = p(ix,iy,NZ) *
     $                 exp( -dzfix/2. * GRAV/(R_AIR*t(ix,iy,NZ)) )

       enddo
      enddo
c
c
c  Define the metric factors that convert coord system diffs to actual distance
c   In cartesian coordinates, the coord are directly distance.
c
      do ixyz=1,NXYZ
        xmet3(ixyz) = 1.
        ymet3(ixyz) = 1.
        zmet3(ixyz) = 1.
      enddo
c
c
c   Compute constants used in viscosity
c
      rmu_0 = 1.8325e-4
      rmu_t0 = 296.16
      rmu_c = 120.
      rmu_const = rmu_0 * (rmu_t0 + rmu_c)
c
c
c  Define vertical profile of atmospheric variables that are 1-D.
c   In this demo:
c    Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
c     Meteorological Tables, in which there is a misprint -- T is deg_K, not deg_C.
c    Thermal conductivity of dry air <thcond> is from Pruppacher and Klett, Eq. 13-16.
c
      do k=1,NZ
        rmu(k) = rmu_const / ( t(1,1,k) + rmu_c ) *
     $              ( t(1,1,k) / rmu_t0 )**1.5d0
        thcond(k) = ( 5.69 + .017*( t(1,1,k) - T0 ) )*4.18e2
      enddo
c
c
c  Define vertical wind speed & diffusion coefficients [in metric cgs units]
c
      do k=1,NZP1
        do ix=1,NX
          do iy=1,NY
            w(ix,iy,k) = 0.e0
            dkx(ix,iy,k) = 0.e5
            dky(ix,iy,k) = 0.e5
            dkz(ix,iy,k) = 0.e1

! From E.Barth - EJL 2-14-13 Commented out diffusion below 
c         Constant diffusion coefficient at/below 90 km
c            if( zl(ix,iy,k) .le. 90.0d5) dkz(ix,iy,k) = 5.d3
c
c            rhoa150 = 4.700000d-06
c            if( zl(ix,iy,k) .eq. 90.0d5)
c     $          y = (dlog10(2.5d13/sqrt(AVG * rhoa150/WTMOL_AIR))
c     $                  - dlog10(5.d3))/(60.d5/dz(ix,iy,k))
c
c            if((zl(ix,iy,k) .gt. 90.0d5) .and.
c     $                          (zl(ix,iy,k) .lt. 150.0d5)) then
c              dkz(ix,iy,k) = 10**(dlog10(5.d3) + y)
c              y = (dlog10(2.5d13/sqrt(AVG * rhoa150/WTMOL_AIR))
c     $                  - dlog10(5.d3))/(60.d5/dzfix) + y
c            endif
c
c            if((zl(ix,iy,k) .ge. 150.0d5) .and.
c     $                          (zl(ix,iy,k) .lt. 400.0d5))
c     $        dkz(ix,iy,k) = 2.5d13/sqrt(AVG * rhoa(ix,iy,k)/WTMOL_AIR)
c
c
c         Constant diffusion coefficient at/above 400 km
c            if( zl(ix,iy,k) .ge. 400.0d5) dkz(ix,iy,k) = 1.d8
c
c         Vertical wind between 320 and 430 km
cc          if( zl(ix,iy,k) .ge. 320.0d5 .and.
cc   $                        zl(ix,iy,k) .lt. 430.0d5)
cc   $        w(ix,iy,k) = 1.25
c          if(k.eq.32) write(*,*) '<initatm> No wind'

cc        Sensitivity Test:  Change diffusion coefficient by constant factor
cc          dkz(ix,iy,k) = 2.*dkz(ix,iy,k)
cc          if( ix.eq.1 .and. iy.eq.1 .and. k.eq.1 )
cc   $       write(*,*) 'Changing eddy diffusion by constant factor !!!'

          enddo
        enddo
      enddo
c
c
c  Define horizontal wind speed [in metric cgs units]
c
      do k=1,NZ
        do ix=1,NX
          do iy=1,NY
            u(ix,iy,k) = 0.
            v(ix,iy,k) = 0.
          enddo
        enddo
      enddo
c
c
c  Return to caller with atm profiles initialized.
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine g_cart_sig
c
c
c  @(#) g_cart_sig.f  Ackerman  Aug-1997
c  This routine handles atmospheric grid initialization
c  for the horizontal=cartesian and vertical=sigma grid coordinates.
c
c  For this coord system and grid:
c   x,y are regular measurable distance in cm.
c   x is positive to the right (east), y is positive north. z is positive up.
c   xmet, ymet metric factors are all unity with no units.
c   z is normalized (unitless) pressure, with fixed pressure (p_top) at top of atm.
c   z=0 at top of atm, z=1 at bottom of atm (positive toward planet surface).
c   zmet metric factor makes z space into vertical distance (altitude)
c
c  Modified slightly for use with generalized coordinates.  Jul-1997, McKie.
c  (From v1.11 single geometrical coord system model's initatm.f routine.)
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter g_cart_sig'
c
c
c  Define bogus values for grid selection parameters that have no meaning
c  for this coordinate system.
c
      rlon0 = -999.
      rlat0 = -999.
      rlat1 = -999.
      rlat2 = -999.
      hemisph = -999.
c
c
c  Define descriptive text for type of grid
c
      gridname = 'cartesian horizontal, sigma vertical'
c
c
c For this example, grid spacing is uniform.
c  Define its east-west, north-south fixed grid box size.
c
      dxfix = 1.e5
      dyfix = 1.e5
c
c
c  Define lower left (southwest) upper right (northest) horizontal
c  domain limits for the grid
c
      dom_llx = 0.
      dom_lly = 0.
      dom_urx = NX * dxfix
      dom_ury = NY * dyfix
c
c
c  Specify horizontal grid coordinates and grid box sizes.
c   <dx> and <dy> are the east-west and north-south grid box thicknesses.
c   <xc>, <yc> are the coord at the center of each grid box.
c   <xl> & <xu> are the lower & upper grid box x coord at y=<yc>.
c   <yl> & <yu> are the lower & upper grid box y coord at x=<xc>.
c   <rlon> & <rlat> are longitude & latitude [deg] at planet surface
c     ( <rlon>, <rlat> usually assumed constant in horiz cartesian coord )
c
      do iy = 1,NY
	      do ix = 1,NX
          rlon(ix,iy) = -90.
          rlat(ix,iy) = 45.
	        do k = 1,NZ
	          dx(ix,iy,k) = dxfix
            xc(ix,iy,k) = dom_llx + ( ix - .5 ) * dx(ix,iy,k)
            xl(ix,iy,k) = dom_llx + ( ix - 1 ) * dx(ix,iy,k)
            xu(ix,iy,k) = dom_llx + ix * dx(ix,iy,k)
	          dy(ix,iy,k) = dyfix
            yc(ix,iy,k) = dom_lly + ( iy - .5 ) * dy(ix,iy,k)
            yl(ix,iy,k) = dom_lly + ( iy - 1 ) * dy(ix,iy,k)
            yu(ix,iy,k) = dom_lly + iy * dy(ix,iy,k)
          enddo
        enddo
      enddo
c
c
c  Define horizontal metric factors that convert coord system diffs to actual distance
c   In cartesian coordinates, the coord are directly distance.
c
      do ixyz=1,NXYZ
        xmet3(ixyz) = 1.
        ymet3(ixyz) = 1.
      enddo
c
c
c  Define vertical sigma values at grid box boundaries.
c   These could be enumerated.
c   This demo uses a simple fractional power function.
c   sigma(k)=((k-1)/NZ)**expon, where 0 < expon < 1.
c   This makes sigma intervals near ground smaller than at top of atm.
c   Use expon closer to 0 to get more dramatically changing sigma layers.
c
      expon = .5
      factor = 1. / float(NZ)
      do k=1,NZP1
        zl2(1,k) = ( float(k-1) * factor )**expon
        do ixy=2,NXY
          zl2(ixy,k) = zl2(1,k)
        enddo
      enddo
c
c
c  Compute sigma values at vertical midpoints of grid boxes & vert box size
c
      do k=1,NZ
        zc2(1,k) = .5 * ( zl2(1,k) + zl2(1,k+1) )
        dz2(1,k) = zl2(1,k+1) - zl2(1,k)
        do ixy=2,NXY
          zc2(ixy,k) = zc2(1,k)
          dz2(ixy,k) = dz2(1,k)
        enddo
      enddo 
c
c
c  Define pressures at bottom and top of model domain [dyne/cm^2]
c
      do iy=1,NY
        do ix=1,NX
          p_surf(ix,iy) = 1496 * RMB2CGS
          p_top(ix,iy) = 900. * RMB2CGS
        enddo
      enddo
c
c
c  Define vertical profile of atmospheric variables that are 3-D.
c   In this demo:
c    The vertical profiles do not vary horizontally.
c    A constant potential temperature <pt> is defined throughout the atm.
c    Air temp <t> is computed as a function of <p> & <pt>.
c    Air density <rhoa> is a function of <t> & <p>, from dry air equa of state.
c    All units are cgs, deg_K.
c
      t_sfc = 94.

      do iy=1,NY
        do ix=1,NX
          t_surf(ix,iy) = t_sfc
          pt_fix = t_sfc * ( PREF / p_surf(ix,iy) )**RKAPPA
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            p(ix,iy,k) = p_top(ix,iy) + zc(ix,iy,k) * pstar
            t(ix,iy,k) = pt_fix * ( p(ix,iy,k) / PREF )**RKAPPA
            rhoa(ix,iy,k) = p(ix,iy,k) / ( R_AIR * t(ix,iy,k) )
            ptc(ix,iy,k) = rhoa(ix,iy,k) * pt_fix
          enddo
        enddo
      enddo
c
c
c  Define <zmet>, the  metric factor that makes <dz> in
c  arbitrary coord system a true distance.
c  See Toon et al., JAS, vol 45, #15, Aug-1988, p 2125.
c
      do iy=1,NY
        do ix=1,NX
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            zmet(ix,iy,k) = pstar / ( GRAV * rhoa(ix,iy,k) )
          enddo
        enddo
      enddo
c
c
c   Compute constants used in viscosity
c
      rmu_0 = 1.8325e-4
      rmu_t0 = 296.16
      rmu_c = 120.
      rmu_const = rmu_0 * (rmu_t0 + rmu_c)
c
c
c  Define vertical profile of atmospheric variables that are 1-D.
c   In this demo:
c    Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
c     Meteorological Tables, in which there is a misprint -- T is deg_K, not deg_C.
c    Thermal conductivity of dry air <thcond> is from Pruppacher and Klett, Eq. 13-16.
c
      do k=1,NZ
        rmu(k) = rmu_const / ( t(1,1,k) + rmu_c ) *
     $              ( t(1,1,k) / rmu_t0 )**1.5d0
        thcond(k) = ( 5.69 + .017*( t(1,1,k) - T0 ) )*4.18e2
      enddo
c
c
c  Define vertical wind speed & diffusion coefficients [in metric cgs units]
c
      do k=1,NZP1
        do ix=1,NX
          do iy=1,NY
            w(ix,iy,k) =  0.
            dkx(ix,iy,k) = 0.e5
            dky(ix,iy,k) = 0.e5
            dkz(ix,iy,k) = 0.e-1
          enddo
        enddo
      enddo
c
c
c  Define horizontal wind speed [in metric cgs units]
c
      do k=1,NZ
        do ix=1,NX
          do iy=1,NY
            u(ix,iy,k) = 0.e4
            v(ix,iy,k) = 0.
          enddo
        enddo
      enddo
c
c
c  Return to caller with atm profiles initialized.
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine g_ll_sig
c
c
c  @(#) g_ll_sig.f  McKie  Jul-1997
c  This routine handles atmospheric grid initialization
c  for horizontal=lon/lat and vertical=sigma grid coordinates.
c
c  For this coord system and grid:
c   x is longitude, y is latitude, both in degrees.
c   x is positive east about Greenwich, y is positive north about equator.
c   x,y Grid node position arrays are in degrees.
c   x,y Grid box size are in degrees.
c   xmet, ymet metric factors make degree space into distance on surface.
c   z is normalized (unitless) pressure, with fixed pressure (p_top) at top of atm.
c   z=0 at top of atm, z=1 at bottom of atm (positive toward planet surface).
c   zmet metric factor makes make z space into vertical distance (altitude)
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter g_ll_sig'
c
c
c  Define bogus values for grid selection parameters that have no meaning
c  for this coordinate system.
c
      rlon0 = -999.
      rlat0 = -999.
      rlat1 = -999.
      rlat2 = -999.
      hemisph = -999.
c
c
c  Define descriptive text for type of grid
c
      gridname = 'longitude & latitude horizontal, sigma vertical'
c
c
c  Define lower left and upper right lon,lat horizontal domain limits
c
      dom_llx = -135.
      dom_urx = -60.
      dom_lly = 15.
      dom_ury = 50.
c
c
c  Compute uniform spacing of grid nodes, independently in lon & lat directions 
c   (This example uses uniform lon & lat spacing, but variable spacing is possible)
c
      dlon = ( dom_urx - dom_llx ) / float(NX)
      dlat = ( dom_ury - dom_lly ) / float(NY)
c
c
c  Compute horizontal grid related things for each grid box:
c   <xc,yc> is lon,lat at grid box center
c   <xl,xu> is lon,lat at left & right edges at yc
c   <yl,yu> is lon,lat at bottom & top edges at xc
c   <dx,dy> is horiz grid box size in lon,lat at xc,yc 
c   <xmet> is metric factor that makes dx in arbitrary coord system a true distance
c   <ymet> is metric factor that makes dy in arbitrary coord system a true distance
c   <rlon> & <rlat> are longitude & latitude [deg] at planet surface
c
      factor = REARTH * DEG2RAD
      do iy=1,NY
        do ix=1,NX
          xc(ix,iy,1) = dom_llx + ( ix - .5 ) * dlon
          yc(ix,iy,1) = dom_lly + ( iy - .5 ) * dlat
          xl(ix,iy,1) = dom_llx + ( ix - 1 ) * dlon
          yl(ix,iy,1) = dom_lly + ( iy - 1 ) * dlat
          xu(ix,iy,1) = dom_llx + ( ix ) * dlon
          yu(ix,iy,1) = dom_lly + ( iy ) * dlat
          dx(ix,iy,1) = dlon
          dy(ix,iy,1) = dlat
          xmet(ix,iy,1) = factor * cos( DEG2RAD * yc(ix,iy,1) )
          ymet(ix,iy,1) = factor
          do k=2,NZ
           xc(ix,iy,k) = xc(ix,iy,1)
           yc(ix,iy,k) = yc(ix,iy,1)
           xl(ix,iy,k) = xl(ix,iy,1)
           yl(ix,iy,k) = yl(ix,iy,1)
           xu(ix,iy,k) = xu(ix,iy,1)
           yu(ix,iy,k) = yu(ix,iy,1)
           dx(ix,iy,k) = dx(ix,iy,1)
           dy(ix,iy,k) = dy(ix,iy,1)
           xmet(ix,iy,k) = xmet(ix,iy,1)
           ymet(ix,iy,k) = ymet(ix,iy,1)
          enddo
          rlon(ix,iy) = xc(ix,iy,NZ)
          rlat(ix,iy) = yc(ix,iy,NZ)
        enddo
      enddo
c
c
c  define vertical sigma values at grid box boundaries.
c   these could be enumerated.
c   this demo uses a simple fractional power function.
c   sigma(k)=((k-1)/NZ)**expon, where 0 < expon < 1.
c   this makes sigma intervals near ground smaller than at top of atm.
c   use expon closer to 0 to get more dramatically changing sigma layers.
c
      expon = .5
      factor = 1. / float(NZ)
      do k=1,NZP1
        zl2(1,k) = ( float(k-1) * factor )**expon
        do ixy=2,NXY
          zl2(ixy,k) = zl2(1,k)
        enddo
      enddo
c
c
c  compute sigma values at vertical midpoints of grid boxes & vert box size
c
      do k=1,NZ
        zc2(1,k) = .5 * ( zl2(1,k) + zl2(1,k+1) )
        dz2(1,k) = zl2(1,k+1) - zl2(1,k)
        do ixy=2,NXY
          zc2(ixy,k) = zc2(1,k)
          dz2(ixy,k) = dz2(1,k)
        enddo
      enddo 
c
c
c  Define pressures at bottom and top of model domain [dyne/cm^2]
c
      do iy=1,NY
        do ix=1,NX
          p_surf(ix,iy) = 1496 * RMB2CGS
          p_top(ix,iy) = 900. * RMB2CGS
        enddo
      enddo
c
c
c  Define vertical profile of atmospheric variables that are 3-D.
c   In this demo:
c    The vertical profiles do not vary horizontally.
c    A constant potential temperature <pt> is defined throughout the atm.
c    Air temp <t> is computed as a function of <p> & <pt>.
c    Air density <rhoa> is a function of <t> & <p>, from dry air equa of state.
c    All units are cgs, deg_K.
c
      t_sfc = 94.

      do iy=1,NY
        do ix=1,NX
          t_surf(ix,iy) = t_sfc
          pt_fix = t_sfc * ( PREF / p_surf(ix,iy) )**RKAPPA
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            p(ix,iy,k) = p_top(ix,iy) + zc(ix,iy,k) * pstar
            t(ix,iy,k) = pt_fix * ( p(ix,iy,k) / PREF )**RKAPPA
            rhoa(ix,iy,k) = p(ix,iy,k) / ( R_AIR * t(ix,iy,k) )
            ptc(ix,iy,k) = rhoa(ix,iy,k) * pt_fix
          enddo
        enddo
      enddo
c
c
c  Define <zmet>, the  metric factor that makes <dz> in
c  arbitrary coord system a true distance.
c  See Toon et al., JAS, vol 45, #15, Aug-1988, p 2125.
c
      do iy=1,NY
        do ix=1,NX
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            zmet(ix,iy,k) = pstar / ( GRAV * rhoa(ix,iy,k) )
          enddo
        enddo
      enddo
c
c
c   Compute constants used in viscosity
c
      rmu_0 = 1.8325e-4
      rmu_t0 = 296.16
      rmu_c = 120.
      rmu_const = rmu_0 * (rmu_t0 + rmu_c)
c
c
c  Define vertical profile of atmospheric variables that are 1-D.
c   In this demo:
c    Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
c     Meteorological Tables, in which there is a misprint -- T is deg_K, not deg_C.
c    Thermal conductivity of dry air <thcond> is from Pruppacher and Klett, Eq. 13-16.
c
      do k=1,NZ
        rmu(k) = rmu_const / ( t(1,1,k) + rmu_c ) *
     $              ( t(1,1,k) / rmu_t0 )**1.5d0
        thcond(k) = ( 5.69 + .017*( t(1,1,k) - T0 ) )*4.18e2
      enddo
c
c
c  Define vertical wind speed & diffusion coefficients [in metric cgs units]
c
      do k=1,NZP1
        do ix=1,NX
          do iy=1,NY
            w(ix,iy,k) = 0.
            dkx(ix,iy,k) = 0.e5
            dky(ix,iy,k) = 0.e5
            dkz(ix,iy,k) = 0.e-1
          enddo
        enddo
      enddo
c
c
c  Define horizontal wind speed [in metric cgs units]
c
      do k=1,NZ
        do ix=1,NX
          do iy=1,NY
            u(ix,iy,k) = 5.e6
            v(ix,iy,k) = 0.
          enddo
        enddo
      enddo
c
c
c  Return to caller with atm profiles initialized
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine g_lc_sig
c
c
c  @(#) g_lc_sig.f  McKie  Jul-1997
c  This routine handles atmospheric grid initialization
c  for horizontal=lambert_conformal and vertical=sigma grid coordinates.
c
c  For this coord system and grid:
c   Planet surface is projected to a cone whose apex is on unit radius planet's
c   polar axis & whose sides pass through 2 latitudes rlat1 & rlat2.  Projection
c   is from center of unit radius planet along straight line through planet
c   surface point to the cone surface.  The cone is then unwrapped along a
c   longitudinal cut line and spread onto the Pu,Pv projection plane, with a
c   central longitude at rlon0.
c
c   x is Pu, y Pv, both unitless.
c   x is positive to the right, y is positive up.
c   x,y Grid node position arrays are mapped back to lon,lat in degrees.
c   Lon is positive east of Greenwich in range [-180,180].
c   Lat is positive north of equator in range [-90,90].
c   x,y Grid box size are in projection plane Pu,Pv units for a unit sphere.
c   xmet, ymet metric factors make Pu,Pv space into distance on planet surface.
c   z is normalized (unitless) pressure, with fixed pressure (p_top) at top of atm.
c   z=0 at top of atm, z=1 at bottom of atm (positive toward planet surface).
c   zmet metric factor makes z space into vertical distance (altitude)
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter g_lc_sig'
c
c
c  Define bogus values for grid selection parameters that have no meaning
c  for this coordinate system.
c
      rlat0 = -999.
c
c
c  Define descriptive text for type of grid
c
      gridname = 'Lambert conformal horizontal, sigma vertical'
c
c
c  Define 2 latitudes at which cone passes through planet's surface
c
      rlat1 = 30.
      rlat2 = 60.
c
c
c  Define central longitude
c
      rlon0 = -100.
c
c
c  Define lon,lat for lower left and upper right lon,lat horizontal domain limits
c   (Note that these are independent of rlon0.  I.e. the domain can be
c    anywhere in the projection plane, not necessarily symmetric with rlon0,
c    or even near rlon0.  In this demo, the domain is chosen symmetrically
c    with rlon0.)
c
      dom_llx = -125.30
      dom_lly = 29.90
      dom_urx = -54.30
      dom_ury = 60.
c
c
c  Compute frequently used forward & reverse projection parameters.
c  Hemisphere constant <hemisph> is chosen automatically based on rlat1 & rlat2.
c   <hemisph> indicates the hemisphere opposite to the projection cone's apex:
c     hemisph = -1.   for northern apex (used for northern hemisph projections)
c     hemisph = +1.   for southern apex (used for southern hemisph projections)
c   <hemisph2> is hemisph/2.
c   <cone> is "cone constant" as a function of rlat1 & rlat2.
c
      call parmlc(rlat1,rlat2, hemisph,hemisph2,cone)
c
c
c  Compute Pu,Pv projections of lon,lat domain limits
c
      call projlc(cone,rlon0,hemisph,hemisph2,
     $            1, dom_llx, dom_lly, Pu_ll, Pv_ll)
      call projlc(cone,rlon0,hemisph,hemisph2,
     $            1, dom_urx, dom_ury, Pu_ur, Pv_ur)
c
c
c  Compute uniform spacing of grid nodes, independently in Pu & Pv directions 
c   (This example uses uniform Pu, Pv spacing, but variable spacing is possible)
c
      dPu = ( Pu_ur - Pu_ll ) / float(NX)
      dPv = ( Pv_ur - Pv_ll ) / float(NY)
c
c
c  Compute preliminary horizontal grid related things for each grid box:
c   <xc,yc> is initially Pu,Pv at grid box center
c   <xl,xu> is initially Pu,Pv at left & right edges at yc
c   <yl,yu> is initially Pu,Pv at bottom & top edges at xc
c   <dx,dy> is horiz grid box size in Pu,Pv at xc,yc 
c
      do iy=1,NY
        do ix=1,NX
          xc(ix,iy,1) = Pu_ll + ( ix - .5 ) * dPu
          yc(ix,iy,1) = Pv_ll + ( iy - .5 ) * dPv
          xl(ix,iy,1) = Pu_ll + ( ix - 1 ) * dPu
          yl(ix,iy,1) = Pv_ll + ( iy - 1 ) * dPv
          xu(ix,iy,1) = Pu_ll + ( ix ) * dPu
          yu(ix,iy,1) = Pv_ll + ( iy ) * dPv
          dx(ix,iy,1) = dPu
          dy(ix,iy,1) = dPv
          do k=2,NZ
           xc(ix,iy,k) = xc(ix,iy,1)
           yc(ix,iy,k) = yc(ix,iy,1)
           xl(ix,iy,k) = xl(ix,iy,1)
           yl(ix,iy,k) = yl(ix,iy,1)
           xu(ix,iy,k) = xu(ix,iy,1)
           yu(ix,iy,k) = yu(ix,iy,1)
           dx(ix,iy,k) = dx(ix,iy,1)
           dy(ix,iy,k) = dy(ix,iy,1)
          enddo
        enddo
      enddo
c
c
c  Compute horizontal metric factors by forming ratio of distance along planet's
c  surface to projected grid box distance in x & y directions. 
c   (u=upper, l=lower)
c
c   Grid box in Pu,Pv space is:
c
c                      rlon_yu, rlat_yu
c                          (xc,yu)
c                        +----*----+
c                        |         |
c                        |         |
c      rlon_xl, rlat_xl  *    +    *  rlon_xu, rlat_xu
c           (xl,yc)      |         |      (xu,yc)
c                        |         |
c                        +----*----+
c                          (xc,yl)
c                      rlon_yl, rlat_yl
c
c   <xmet> is metric factor that makes dx in arbitrary coord system a true distance
c   <ymet> is metric factor that makes dy in arbitrary coord system a true distance
c
c  After computing metric factor, reverse project the position info back from
c  Pu,Pv to lon,lat space.
c  (This is position info intended for history output for post-processing)
c   <xc,yc> becomes lon,lat at grid box center
c   <xl,xu> becomes lon,lat at left & right edges at yc
c   <yl,yu> becomes lon,lat at bottom & top edges at xc
c   <dx,dy> remains horiz grid box size in Pu,Pv at xc,yc 
c
      do ixyz=1,NXYZ

        call invplc(cone,rlon0,hemisph,
     $              1, xl3(ixyz),yc3(ixyz), rlon_xl,rlat_xl)
        call invplc(cone,rlon0,hemisph,
     $              1, xu3(ixyz),yc3(ixyz), rlon_xu,rlat_xu)
        call invplc(cone,rlon0,hemisph,
     $              1, xc3(ixyz),yl3(ixyz), rlon_yl,rlat_yl)
        call invplc(cone,rlon0,hemisph,
     $              1, xc3(ixyz),yu3(ixyz), rlon_yu,rlat_yu)

        ds_x = REARTH * sfirdis(rlon_xl,rlat_xl, rlon_xu,rlat_xu)
        ds_y = REARTH * sfirdis(rlon_yl,rlat_yl, rlon_yu,rlat_yu)

        xmet3(ixyz) = ds_x / dx3(ixyz)
        ymet3(ixyz) = ds_y / dy3(ixyz)

        xl3(ixyz) = rlon_xl
        xu3(ixyz) = rlon_xu
        yl3(ixyz) = rlat_yl
        yu3(ixyz) = rlat_yu

        call invplc(cone,rlon0,hemisph,
     $              1, xc3(ixyz),yc3(ixyz), xc3(ixyz),yc3(ixyz))

      enddo
c
c
c  Define 2-D longitude & latitude [deg] at planet surface, used internally by model
c
      do iy=1,NY
        do ix=1,NX
          rlon(ix,iy) = xc(ix,iy,NZ)
          rlat(ix,iy) = yc(ix,iy,NZ)
        enddo
      enddo
c
c
c  Define vertical sigma values at grid box boundaries.
c   These could be enumerated.
c   This demo uses a simple fractional power function.
c   sigma(k)=((k-1)/NZ)**expon, where 0 < expon < 1.
c   This makes sigma intervals near ground smaller than at top of atm.
c   Use expon closer to 0 to get more dramatically changing sigma layers.
c
      expon = .5
      factor = 1. / float(NZ)
      do k=1,NZP1
        zl2(1,k) = ( float(k-1) * factor )**expon
        do ixy=2,NXY
          zl2(ixy,k) = zl2(1,k)
        enddo
      enddo
c
c
c  Compute sigma values at vertical midpoints of grid boxes & vert box size
c
      do k=1,NZ
        zc2(1,k) = .5 * ( zl2(1,k) + zl2(1,k+1) )
        dz2(1,k) = zl2(1,k+1) - zl2(1,k)
        do ixy=2,NXY
          zc2(ixy,k) = zc2(1,k)
          dz2(ixy,k) = dz2(1,k)
        enddo
      enddo 
c
c
c  Define pressures at bottom and top of model domain [dyne/cm^2]
c
      do iy=1,NY
        do ix=1,NX
          p_surf(ix,iy) = 1496 * RMB2CGS
          p_top(ix,iy) = 900. * RMB2CGS
        enddo
      enddo
c
c
c  Define vertical profile of atmospheric variables that are 3-D.
c   In this demo:
c    The vertical profiles do not vary horizontally.
c    A constant potential temperature <pt> is defined throughout the atm.
c    Air temp <t> is computed as a function of <p> & <pt>.
c    Air density <rhoa> is a function of <t> & <p>, from dry air equa of state.
c    All units are cgs, deg_K.
c
      t_sfc = 94.

      do iy=1,NY
        do ix=1,NX
          t_surf(ix,iy) = t_sfc
          pt_fix = t_sfc * ( PREF / p_surf(ix,iy) )**RKAPPA
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            p(ix,iy,k) = p_top(ix,iy) + zc(ix,iy,k) * pstar
            t(ix,iy,k) = pt_fix * ( p(ix,iy,k) / PREF )**RKAPPA
            rhoa(ix,iy,k) = p(ix,iy,k) / ( R_AIR * t(ix,iy,k) )
            ptc(ix,iy,k) = rhoa(ix,iy,k) * pt_fix
          enddo
        enddo
      enddo
c
c
c  Define <zmet>, the  metric factor that makes <dz> in
c  arbitrary coord system a true distance.
c  See Toon et al., JAS, vol 45, #15, Aug-1988, p 2125.
c
      do iy=1,NY
        do ix=1,NX
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            zmet(ix,iy,k) = pstar / ( GRAV * rhoa(ix,iy,k) )
          enddo
        enddo
      enddo
c
c
c   Compute constants used in viscosity
c
      rmu_0 = 1.8325e-4
      rmu_t0 = 296.16
      rmu_c = 120.
      rmu_const = rmu_0 * (rmu_t0 + rmu_c)
c
c
c  Define vertical profile of atmospheric variables that are 1-D.
c   In this demo:
c    Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
c     Meteorological Tables, in which there is a misprint -- T is deg_K, not deg_C.
c    Thermal conductivity of dry air <thcond> is from Pruppacher and Klett, Eq. 13-16.
c
      do k=1,NZ
        rmu(k) = rmu_const / ( t(1,1,k) + rmu_c ) *
     $              ( t(1,1,k) / rmu_t0 )**1.5d0
        thcond(k) = ( 5.69 + .017*( t(1,1,k) - T0 ) )*4.18e2
      enddo
c
c
c  Define vertical wind speed & diffusion coefficients [in metric cgs units]
c
      do k=1,NZP1
        do ix=1,NX
          do iy=1,NY
            w(ix,iy,k) = 0.
            dkx(ix,iy,k) = 0.e5
            dky(ix,iy,k) = 0.e5
            dkz(ix,iy,k) = 0.e-1
          enddo
        enddo
      enddo
c
c
c  Define horizontal wind speed [in metric cgs units]
c
      do k=1,NZ
        do ix=1,NX
          do iy=1,NY
            u(ix,iy,k) = 0.e4
            v(ix,iy,k) = 0.
          enddo
        enddo
      enddo
c
c
c  Return to caller with atm profiles initialized
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine g_ps_sig
c
c
c  @(#) g_ps_sig.f  McKie  Jul-1997
c  This routine handles atmospheric grid initialization
c  for horizontal=polar_stereographic and vertical=sigma grid coordinates.
c
c  For this coord system and grid:
c   Planet surface is projected to a plane tangent to a unit radius planet'
c   surface at latitude rlat0, usually 90 for north pole point or -90 for
c   south pole point, by straight line from pole point specified by
c   hemisph (hemisph=-1 for south pole, +1 for north) through planet surface
c   point to Pu,Pv projection plane.  The plane is oriented with central
c   longitude rlon0 on negative Pv axis. 
c
c   x is Pu, y Pv, both unitless.
c   x is positive to the right, y is positive up.
c   x,y Grid node position arrays are mapped back to lon,lat in degrees.
c   Lon is positive east of Greenwich in range [-180,180].
c   Lat is positive north of equator in range [-90,90].
c   x,y grid box size are in projection plane Pu,Pv units for a unit sphere.
c   xmet, ymet metric factors make Pu,Pv space into distance on planet surface.
c   z is normalized (unitless) pressure, with fixed pressure (p_top) at top of atm.
c   z=0 at top of atm, z=1 at bottom of atm (positive toward planet surface).
c   zmet metric factor makes z space into vertical distance (altitude)
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter g_lc_sig'
c
c
c  Define bogus values for grid selection parameters that have no meaning
c  for this coordinate system.
c
      rlat1 = -999.
      rlat2 = -999.
c
c
c  Define descriptive text for type of grid
c
      gridname = 'Polar stereographic horizontal, sigma vertical'
c
c
c  Define latitude at which projection plane is tangent to planet
c
      rlat0 = 90.
c
c
c  Define indication of which pole the projection line begins at:
c    <hemisph> = -1.   for south pole point (used for northern hemisph projections)
c    <hemisph> = +1.   for north pole point (used for southern hemisph projections)
c    <rlat0>=90 usually is paired with <hemisph>=-1.
c    <rlat0>=-90 usually is paired with <hemisph>=1.
c
      hemisph = -1.
c
c
c  Define central longitude
c
      rlon0 = -90.
c
c
c  Define lon,lat for lower left and upper right lon,lat horizontal domain limits
c   (Note that these are independent of rlon0.  I.e. the domain can be
c    anywhere in the projection plane, not necessarily symmetric with rlon0,
c    or even near rlon0.  In this demo, the domain is chosen somewhat
c    symmetrically with rlon0.)
c
      dom_llx = -135.00
      dom_lly = 45.
      dom_urx = 45.
      dom_ury = 45.
c
c
c  Compute frequently used forward & reverse projection parameters.
c   <factps> is projection factor computed as a function of rlat0 & hemisph.
c
      call parmps(rlat0,hemisph, factps)
c
c
c  Compute Pu,Pv projections of lon,lat domain limits
c
      call projps(rlon0,hemisph,factps,
     $            1, dom_llx, dom_lly, Pu_ll, Pv_ll)
      call projps(rlon0,hemisph,factps,
     $            1, dom_urx, dom_ury, Pu_ur, Pv_ur)
c
c
c  Compute uniform spacing of grid nodes, independently in Pu & Pv directions 
c   (This example uses uniform Pu, Pv spacing, but variable spacing is possible)
c
      dPu = ( Pu_ur - Pu_ll ) / float(NX)
      dPv = ( Pv_ur - Pv_ll ) / float(NY)
c
c
c  Compute preliminary horizontal grid related things for each grid box:
c   <xc,yc> is initially Pu,Pv at grid box center
c   <xl,xu> is initially Pu,Pv at left & right edges at yc
c   <yl,yu> is initially Pu,Pv at bottom & top edges at xc
c   <dx,dy> is horiz grid box size in Pu,Pv at xc,yc 
c
      do iy=1,NY
        do ix=1,NX
          xc(ix,iy,1) = Pu_ll + ( ix - .5 ) * dPu
          yc(ix,iy,1) = Pv_ll + ( iy - .5 ) * dPv
          xl(ix,iy,1) = Pu_ll + ( ix - 1 ) * dPu
          yl(ix,iy,1) = Pv_ll + ( iy - 1 ) * dPv
          xu(ix,iy,1) = Pu_ll + ( ix ) * dPu
          yu(ix,iy,1) = Pv_ll + ( iy ) * dPv
          dx(ix,iy,1) = dPu
          dy(ix,iy,1) = dPv
          do k=2,NZ
           xc(ix,iy,k) = xc(ix,iy,1)
           yc(ix,iy,k) = yc(ix,iy,1)
           xl(ix,iy,k) = xl(ix,iy,1)
           yl(ix,iy,k) = yl(ix,iy,1)
           xu(ix,iy,k) = xu(ix,iy,1)
           yu(ix,iy,k) = yu(ix,iy,1)
           dx(ix,iy,k) = dx(ix,iy,1)
           dy(ix,iy,k) = dy(ix,iy,1)
          enddo
        enddo
      enddo
c
c
c  Compute horizontal metric factors by forming ratio of distance along planet's
c  surface to projected grid box distance in x & y directions. 
c   (u=upper, l=lower)
c
c   Grid box in Pu,Pv space is:
c
c                      rlon_yu, rlat_yu
c                          (xc,yu)
c                        +----*----+
c                        |         |
c                        |         |
c      rlon_xl, rlat_xl  *    +    *  rlon_xu, rlat_xu
c           (xl,yc)      |         |      (xu,yc)
c                        |         |
c                        +----*----+
c                          (xc,yl)
c                      rlon_yl, rlat_yl
c
c   <xmet> is metric factor that makes dx in arbitrary coord system a true distance
c   <ymet> is metric factor that makes dy in arbitrary coord system a true distance
c
c  After computing metric factor, reverse project the position info back from
c  Pu,Pv to lon,lat space.
c  (This is position info intended for history output for post-processing)
c   <xc,yc> becomes lon,lat at grid box center
c   <xl,xu> becomes lon,lat at left & right edges at yc
c   <yl,yu> becomes lon,lat at bottom & top edges at xc
c   <dx,dy> remains horiz grid box size in Pu,Pv at xc,yc 
c
      do ixyz=1,NXYZ

        call invpps(rlon0,hemisph,factps,
     $              1, xl3(ixyz),yc3(ixyz), rlon_xl,rlat_xl)
        call invpps(rlon0,hemisph,factps,
     $              1, xu3(ixyz),yc3(ixyz), rlon_xu,rlat_xu)
        call invpps(rlon0,hemisph,factps,
     $              1, xc3(ixyz),yl3(ixyz), rlon_yl,rlat_yl)
        call invpps(rlon0,hemisph,factps,
     $              1, xc3(ixyz),yu3(ixyz), rlon_yu,rlat_yu)

        ds_x = REARTH * sfirdis(rlon_xl,rlat_xl, rlon_xu,rlat_xu)
        ds_y = REARTH * sfirdis(rlon_yl,rlat_yl, rlon_yu,rlat_yu)

        xmet3(ixyz) = ds_x / dx3(ixyz)
        ymet3(ixyz) = ds_y / dy3(ixyz)

        xl3(ixyz) = rlon_xl
        xu3(ixyz) = rlon_xu
        yl3(ixyz) = rlat_yl
        yu3(ixyz) = rlat_yu

        call invpps(rlon0,hemisph,factps,
     $              1, xc3(ixyz),yc3(ixyz), xc3(ixyz),yc3(ixyz))

      enddo
c
c
c  Define 2-D longitude & latitude [deg] at planet surface, used internally by model
c
      do iy=1,NY
        do ix=1,NX
          rlon(ix,iy) = xc(ix,iy,NZ)
          rlat(ix,iy) = yc(ix,iy,NZ)
        enddo
      enddo
c
c
c  Define vertical sigma values at grid box boundaries.
c   These could be enumerated.
c   This demo uses a simple fractional power function.
c   sigma(k)=((k-1)/NZ)**expon, where 0 < expon < 1.
c   This makes sigma intervals near ground smaller than at top of atm.
c   Use expon closer to 0 to get more dramatically changing sigma layers.
c
      expon = .5
      factor = 1. / float(NZ)
      do k=1,NZP1
        zl2(1,k) = ( float(k-1) * factor )**expon
        do ixy=2,NXY
          zl2(ixy,k) = zl2(1,k)
        enddo
      enddo
c
c
c  Compute sigma values at vertical midpoints of grid boxes & vert box size
c
      do k=1,NZ
        zc2(1,k) = .5 * ( zl2(1,k) + zl2(1,k+1) )
        dz2(1,k) = zl2(1,k+1) - zl2(1,k)
        do ixy=2,NXY
          zc2(ixy,k) = zc2(1,k)
          dz2(ixy,k) = dz2(1,k)
        enddo
      enddo 
c
c
c  Define pressures at bottom and top of model domain [dyne/cm^2]
c
      do iy=1,NY
        do ix=1,NX
          p_surf(ix,iy) = 1013.5 * RMB2CGS
          p_top(ix,iy) = 900. * RMB2CGS
        enddo
      enddo
c
c
c  Define vertical profile of atmospheric variables that are 3-D.
c   In this demo:
c    The vertical profiles do not vary horizontally.
c    A constant potential temperature <pt> is defined throughout the atm.
c    Air temp <t> is computed as a function of <p> & <pt>.
c    Air density <rhoa> is a function of <t> & <p>, from dry air equa of state.
c    All units are cgs, deg_K.
c
      t_sfc = 288.

      do iy=1,NY
        do ix=1,NX
          t_surf(ix,iy) = t_sfc
          pt_fix = t_sfc * ( PREF / p_surf(ix,iy) )**RKAPPA
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            p(ix,iy,k) = p_top(ix,iy) + zc(ix,iy,k) * pstar
            t(ix,iy,k) = pt_fix * ( p(ix,iy,k) / PREF )**RKAPPA
            rhoa(ix,iy,k) = p(ix,iy,k) / ( R_AIR * t(ix,iy,k) )
            ptc(ix,iy,k) = rhoa(ix,iy,k) * pt_fix
          enddo
        enddo
      enddo
c
c
c  Define <zmet>, the  metric factor that makes <dz> in
c  arbitrary coord system a true distance.
c  See Toon et al., JAS, vol 45, #15, Aug-1988, p 2125.
c
      do iy=1,NY
        do ix=1,NX
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            zmet(ix,iy,k) = pstar / ( GRAV * rhoa(ix,iy,k) )
          enddo
        enddo
      enddo
c
c
c   Compute constants used in viscosity
c
      rmu_0 = 1.8325e-4
      rmu_t0 = 296.16
      rmu_c = 120.
      rmu_const = rmu_0 * (rmu_t0 + rmu_c)
c
c
c  Define vertical profile of atmospheric variables that are 1-D.
c   In this demo:
c    Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
c     Meteorological Tables, in which there is a misprint -- T is deg_K, not deg_C.
c    Thermal conductivity of dry air <thcond> is from Pruppacher and Klett, Eq. 13-16.
c
      do k=1,NZ
        rmu(k) = rmu_const / ( t(1,1,k) + rmu_c ) *
     $              ( t(1,1,k) / rmu_t0 )**1.5d0
        thcond(k) = ( 5.69 + .017*( t(1,1,k) - T0 ) )*4.18e2
      enddo
c
c
c  Define vertical wind speed & diffusion coefficients [in metric cgs units]
c
      do k=1,NZP1
        do ix=1,NX
          do iy=1,NY
            w(ix,iy,k) = 0.
            dkx(ix,iy,k) = 0.e5
            dky(ix,iy,k) = 0.e5
            dkz(ix,iy,k) = 0.e-1
          enddo
        enddo
      enddo
c
c
c  Define horizontal wind speed [in metric cgs units]
c
      do k=1,NZ
        do ix=1,NX
          do iy=1,NY
            u(ix,iy,k) = 0.e4
            v(ix,iy,k) = 0.
          enddo
        enddo
      enddo
c
c
c  Return to caller with atm profiles initialized
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine g_me_sig
c
c
c  @(#) g_me_sig.f  McKie  Jul-1997
c  This routine handles atmospheric grid initialization
c  for horizontal=Mercator and vertical=sigma grid coordinates.
c
c  For this coord system and grid:
c   Planet surface is projected onto a cylinder whose axis is coincident with
c   a unit radius planet's polar axis.  The cylinder's radius may be less than
c   or same as the unit radius of the planet.
c   The northern hemisphere latitude at which the cyclinder cuts the
c   surface of the planet is at latitude rlat0.
c   Projection is from the center of the planet along a straight line
c   through point on surface of planet, and continuing to surface of cylinder.
c   The cycliner is then unrolled into the projection plane Pu,Pv by a
c   cut along the image of a longitude, and the plane shifted so that
c   longitude rlon0 is in the center of the plane.
c
c   x is Pu, y Pv, both unitless.
c   x is positive to the right, y is positive up.
c   x,y grid node position arrays are mapped back to lon,lat in degrees.
c   Lon is positive east of Greenwich in range [-180,180].
c   Lat is positive north of equator in range [-90,90].
c   x,y grid box size are in projection plane Pu,Pv units for a unit sphere.
c   xmet, ymet metric factors make Pu,Pv space into distance on planet surface.
c   z is normalized (unitless) pressure, with fixed pressure (p_top) at top of atm.
c   z=0 at top of atm, z=1 at bottom of atm (positive toward planet surface).
c   zmet metric factor makes z space into vertical distance (altitude)
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter g_lc_sig'
c
c
c  Define bogus values for grid selection parameters that have no meaning
c  for this coordinate system.
c
      rlat1 = -999.
      rlat2 = -999.
      hemisph = -999.
c
c
c  Define descriptive text for type of grid
c
      gridname = 'Mercator horizontal, sigma vertical'
c
c
c  Define latitude at which projection cylinder cuts through surface of planet
c
      rlat0 = 0.
c
c
c  Define central longitude
c
      rlon0 = -90.
c
c
c  Define lon,lat for lower left & upper right lon,lat horizontal domain limits
c   (Note that these are independent of rlon0.  I.e. the domain can be
c    anywhere in the projection plane, not necessarily symmetric with rlon0,
c    or even near rlon0.  In this demo, the domain is chosen symmetrically
c    with rlon0.)
c
      dom_llx = -135.00
      dom_lly = -45.
      dom_urx = -45.
      dom_ury = 45.
c
c
c  Compute frequently used forward & reverse projection parameters.
c   <factme> is projection factor computed as a function of rlat0.
c
      call parmme(rlat0, factme)
c
c
c  Compute Pu,Pv projections of lon,lat domain limits
c
      call projme(rlon0,factme,
     $            1, dom_llx, dom_lly, Pu_ll, Pv_ll)
      call projme(rlon0,factme,
     $            1, dom_urx, dom_ury, Pu_ur, Pv_ur)
c
c
c  Compute uniform spacing of grid nodes, independently in Pu & Pv directions 
c   (This example uses uniform Pu, Pv spacing, but variable spacing is possible)
c
      dPu = ( Pu_ur - Pu_ll ) / float(NX)
      dPv = ( Pv_ur - Pv_ll ) / float(NY)
c
c
c  Compute preliminary horizontal grid related things for each grid box:
c   <xc,yc> is initially Pu,Pv at grid box center
c   <xl,xu> is initially Pu,Pv at left & right edges at yc
c   <yl,yu> is initially Pu,Pv at bottom & top edges at xc
c   <dx,dy> is horiz grid box size in Pu,Pv at xc,yc 
c
      do iy=1,NY
        do ix=1,NX
          xc(ix,iy,1) = Pu_ll + ( ix - .5 ) * dPu
          yc(ix,iy,1) = Pv_ll + ( iy - .5 ) * dPv
          xl(ix,iy,1) = Pu_ll + ( ix - 1 ) * dPu
          yl(ix,iy,1) = Pv_ll + ( iy - 1 ) * dPv
          xu(ix,iy,1) = Pu_ll + ( ix ) * dPu
          yu(ix,iy,1) = Pv_ll + ( iy ) * dPv
          dx(ix,iy,1) = dPu
          dy(ix,iy,1) = dPv
          do k=2,NZ
           xc(ix,iy,k) = xc(ix,iy,1)
           yc(ix,iy,k) = yc(ix,iy,1)
           xl(ix,iy,k) = xl(ix,iy,1)
           yl(ix,iy,k) = yl(ix,iy,1)
           xu(ix,iy,k) = xu(ix,iy,1)
           yu(ix,iy,k) = yu(ix,iy,1)
           dx(ix,iy,k) = dx(ix,iy,1)
           dy(ix,iy,k) = dy(ix,iy,1)
          enddo
        enddo
      enddo
c
c
c  Compute horizontal metric factors by forming ratio of distance along planet's
c  surface to projected grid box distance in x & y directions. 
c   (u=upper, l=lower)
c
c   Grid box in Pu,Pv space is:
c
c                      rlon_yu, rlat_yu
c                          (xc,yu)
c                        +----*----+
c                        |         |
c                        |         |
c      rlon_xl, rlat_xl  *    +    *  rlon_xu, rlat_xu
c           (xl,yc)      |         |      (xu,yc)
c                        |         |
c                        +----*----+
c                          (xc,yl)
c                      rlon_yl, rlat_yl
c
c   <xmet> is metric factor that makes dx in arbitrary coord system a true distance
c   <ymet> is metric factor that makes dy in arbitrary coord system a true distance
c
c  After computing metric factor, reverse project the position info back from
c  Pu,Pv to lon,lat space.
c  (This is position info intended for history output for post-processing)
c   <xc,yc> becomes lon,lat at grid box center
c   <xl,xu> becomes lon,lat at left & right edges at yc
c   <yl,yu> becomes lon,lat at bottom & top edges at xc
c   <dx,dy> remains horiz grid box size in Pu,Pv at xc,yc 
c
      do ixyz=1,NXYZ

        call invpme(rlon0,factme,
     $              1, xl3(ixyz),yc3(ixyz), rlon_xl,rlat_xl)
        call invpme(rlon0,factme,
     $              1, xu3(ixyz),yc3(ixyz), rlon_xu,rlat_xu)
        call invpme(rlon0,factme,
     $              1, xc3(ixyz),yl3(ixyz), rlon_yl,rlat_yl)
        call invpme(rlon0,factme,
     $              1, xc3(ixyz),yu3(ixyz), rlon_yu,rlat_yu)

        ds_x = REARTH * sfirdis(rlon_xl,rlat_xl, rlon_xu,rlat_xu)
        ds_y = REARTH * sfirdis(rlon_yl,rlat_yl, rlon_yu,rlat_yu)

        xmet3(ixyz) = ds_x / dx3(ixyz)
        ymet3(ixyz) = ds_y / dy3(ixyz)

        xl3(ixyz) = rlon_xl
        xu3(ixyz) = rlon_xu
        yl3(ixyz) = rlat_yl
        yu3(ixyz) = rlat_yu

        call invpme(rlon0,factme,
     $              1, xc3(ixyz),yc3(ixyz), xc3(ixyz),yc3(ixyz))

      enddo
c
c
c  Define 2-D longitude & latitude [deg] at planet surface, used internally by model
c
      do iy=1,NY
        do ix=1,NX
          rlon(ix,iy) = xc(ix,iy,NZ)
          rlat(ix,iy) = yc(ix,iy,NZ)
        enddo
      enddo
c
c
c  Define vertical sigma values at grid box boundaries.
c   These could be enumerated.
c   This demo uses a simple fractional power function.
c   sigma(k)=((k-1)/NZ)**expon, where 0 < expon < 1.
c   This makes sigma intervals near ground smaller than at top of atm.
c   Use expon closer to 0 to get more dramatically changing sigma layers.
c
      expon = .5
      factor = 1. / float(NZ)
      do k=1,NZP1
        zl2(1,k) = ( float(k-1) * factor )**expon
        do ixy=2,NXY
          zl2(ixy,k) = zl2(1,k)
        enddo
      enddo
c
c
c  Compute sigma values at vertical midpoints of grid boxes & vert box size
c
      do k=1,NZ
        zc2(1,k) = .5 * ( zl2(1,k) + zl2(1,k+1) )
        dz2(1,k) = zl2(1,k+1) - zl2(1,k)
        do ixy=2,NXY
          zc2(ixy,k) = zc2(1,k)
          dz2(ixy,k) = dz2(1,k)
        enddo
      enddo 
c
c
c  Define pressures at bottom and top of model domain [dyne/cm^2]
c
      do iy=1,NY
        do ix=1,NX
          p_surf(ix,iy) = 1013.5 * RMB2CGS
          p_top(ix,iy) = 900. * RMB2CGS
        enddo
      enddo
c
c
c  Define vertical profile of atmospheric variables that are 3-D.
c   In this demo:
c    The vertical profiles do not vary horizontally.
c    A constant potential temperature <pt> is defined throughout the atm.
c    Air temp <t> is computed as a function of <p> & <pt>.
c    Air density <rhoa> is a function of <t> & <p>, from dry air equa of state.
c    All units are cgs, deg_K.
c
      t_sfc = 94.

      do iy=1,NY
        do ix=1,NX
          t_surf(ix,iy) = t_sfc
          pt_fix = t_sfc * ( PREF / p_surf(ix,iy) )**RKAPPA
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            p(ix,iy,k) = p_top(ix,iy) + zc(ix,iy,k) * pstar
            t(ix,iy,k) = pt_fix * ( p(ix,iy,k) / PREF )**RKAPPA
            rhoa(ix,iy,k) = p(ix,iy,k) / ( R_AIR * t(ix,iy,k) )
            ptc(ix,iy,k) = rhoa(ix,iy,k) * pt_fix
          enddo
        enddo
      enddo
c
c
c  Define <zmet>, the  metric factor that makes <dz> in
c  arbitrary coord system a true distance.
c  See Toon et al., JAS, vol 45, #15, Aug-1988, p 2125.
c
      do iy=1,NY
        do ix=1,NX
          pstar = p_surf(ix,iy) - p_top(ix,iy)
          do k=1,NZ
            zmet(ix,iy,k) = pstar / ( GRAV * rhoa(ix,iy,k) )
          enddo
        enddo
      enddo
c
c
c   Compute constants used in viscosity
c
      rmu_0 = 1.8325e-4
      rmu_t0 = 296.16
      rmu_c = 120.
      rmu_const = rmu_0 * (rmu_t0 + rmu_c)
c
c
c  Define vertical profile of atmospheric variables that are 1-D.
c   In this demo:
c    Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
c     Meteorological Tables, in which there is a misprint -- T is deg_K, not deg_C.
c    Thermal conductivity of dry air <thcond> is from Pruppacher and Klett, Eq. 13-16.
c
      do k=1,NZ
        rmu(k) = rmu_const / ( t(1,1,k) + rmu_c ) *
     $              ( t(1,1,k) / rmu_t0 )**1.5d0
        thcond(k) = ( 5.69 + .017*( t(1,1,k) - T0 ) )*4.18e2
      enddo
c
c
c  Define vertical wind speed & diffusion coefficients [in metric cgs units]
c
      do k=1,NZP1
        do ix=1,NX
          do iy=1,NY
            w(ix,iy,k) = 0.
            dkx(ix,iy,k) = 0.e5
            dky(ix,iy,k) = 0.e5
            dkz(ix,iy,k) = 0.e-1
          enddo
        enddo
      enddo
c
c
c  Define horizontal wind speed [in metric cgs units]
c
      do k=1,NZ
        do ix=1,NX
          do iy=1,NY
            u(ix,iy,k) = 0.e4
            v(ix,iy,k) = 0.
          enddo
        enddo
      enddo
c
c
c  Return to caller with atm profiles initialized
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine g_ctm_ll_hyb
c
c
c  This routine handles atmospheric grid initialization
c  for the grid coordinates:
c     horizontal:  Lat/Lon
c     vertical:    NCAR CTM hybrid grid coordinates.
c
c  For this coord system and grid:
c   x,y grid box size are in degrees.
c   x is positive to the right (east), y is positive north. z is positive up.
c   xmet, ymet metric factors make degree space into distance on surface.
c
c   z is the hybrid coordinate eta
c     z is 0 at infinite distance
c     z is 1 at the surface
c     z is non-zero at the top of atm.
c   zmet metric factor makes z space into vertical distance (altitude)
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter g_ctm_ll_hyb'
c
c
c  Define bogus values for grid selection parameters that have no meaning
c  for this coordinate system.
c
      rlon0 = -999.
      rlat0 = -999.
      rlat1 = -999.
      rlat2 = -999.
      hemisph = -999.
c
c
c  Define descriptive text for type of grid
c
      gridname = 'long & lat horizontal, CTM hybrid vertical'
c
c
c  Compute uniform spacing of grid nodes, independently in lon & lat directions
c   (This example uses uniform lon & lat spacing, but variable spacing is possible)
c
c   NOTE: dom_urx, dom_ury, dom_llx, dom_lly are set in calling routine.
c
c 
      dlon = ( dom_urx - dom_llx ) / float(NX)
      dlat = ( dom_ury - dom_lly ) / float(NY)
c
c
c  Compute horizontal grid related things for each grid box:
c   <xc,yc> is lon,lat at grid box center
c   <xl,xu> is lon,lat at left & right edges at yc
c   <yl,yu> is lon,lat at bottom & top edges at xc
c   <dx,dy> is horiz grid box size in lon,lat at xc,yc
c   <xmet> is metric factor that makes dx in arbitrary coord system a true distance
c   <ymet> is metric factor that makes dy in arbitrary coord system a true distance
c   <rlon> & <rlat> are longitude & latitude [deg] at planet surface
c
c   NOTE: xmet will be reset each timestep as the latitude changes
c
c
      factor = REARTH * DEG2RAD
      do iy=1,NY
        do ix=1,NX
          xc(ix,iy,1) = dom_llx + ( ix - .5 ) * dlon
          yc(ix,iy,1) = dom_lly + ( iy - .5 ) * dlat
          xl(ix,iy,1) = dom_llx + ( ix - 1 ) * dlon
          yl(ix,iy,1) = dom_lly + ( iy - 1 ) * dlat
          xu(ix,iy,1) = dom_llx + ( ix ) * dlon
          yu(ix,iy,1) = dom_lly + ( iy ) * dlat
          dx(ix,iy,1) = dlon
          dy(ix,iy,1) = dlat
          xmet(ix,iy,1) = factor * cos( DEG2RAD * yc(ix,iy,1) )
          ymet(ix,iy,1) = factor
          do k=2,NZ
           xc(ix,iy,k) = xc(ix,iy,1)
           yc(ix,iy,k) = yc(ix,iy,1)
           xl(ix,iy,k) = xl(ix,iy,1)
           yl(ix,iy,k) = yl(ix,iy,1)
           xu(ix,iy,k) = xu(ix,iy,1)
           yu(ix,iy,k) = yu(ix,iy,1)
           dx(ix,iy,k) = dx(ix,iy,1)
           dy(ix,iy,k) = dy(ix,iy,1)
           xmet(ix,iy,k) = xmet(ix,iy,1)
           ymet(ix,iy,k) = ymet(ix,iy,1)
          enddo
          rlon(ix,iy) = xc(ix,iy,NZ)
          rlat(ix,iy) = yc(ix,iy,NZ)
        enddo
      enddo
c
c
c  Compute carma's hybrid values at vertical midpoints of grid boxes & vert box size
c
      do k = 1,NZ
        do ixy = 1,NXY
          dz2(ixy,k) = abs(zl2(ixy, k+1) - zl2(ixy, k))
        enddo
      enddo
c
c
c  Define pressures at bottom (planet surface) and top of model domain [dyne/cm^2]
c   Bottom could be anything .ge. pres_b_ctm(0)
c   Top must be pres_b_ctm(N_VERT_CTM)
c
!      if( pl2(1, NZP1) .lt. 8.e5 .or. pl2(1, NZP1) .gt. 1.6e6 ) then
!        write(LUNOPRT,*) 'initatm: cgs units of pressure, please.'
!        write(LUNOPRT,*) pl2
!        call endcarma
!      endif

      if( pl2(1, NZP1) .lt. pl2(1, 1) )then
        write(LUNOPRT,*) 'initatm: I know this is against convention,'
        write(LUNOPRT,*) '  but please make layer 1 the bottom-most'
        write(LUNOPRT,*) '  layer of the atmosphere.'       
        call endcarma
      endif
!
      do ixy=1,NXY
        p_surf2(ixy) = pl2(ixy, NZP1)
        p_top2(ixy) = pl2(ixy, 1)
      enddo
c
c
c  Initialize 2-D surface temperature [deg K]
c
!      t_sfc = 94.d0  !EJL
!      do ix=1,NX
!        do iy=1,NY
!          t_surf(ix,iy) = t_sfc
!        enddo
!      enddo
c
c
c  Initialize 3-D temperature [deg K]
c   (For this demo, use an idealized vertical temp vs pres profile
c
c  Note: we're not using this because the temperature is passed
c  from the CTM.
c
c
c      write(LUNOPRT,*) 'initatm: k, alt, eta, temperature'
!      sounding = t_sfc
!      scaleheight = 40.e5  !EJL
!      rlapse = 1.e-5  ! from McKay(1997) Icarus
!      oldalt = 0.
!
!      do k = NZ, 1, -1
!        alt = - scaleheight * log(p2(1, k) / pl2(1, NZP1)) 
!        if (p2(1, k) .gt. 2.e5) then
!          sounding = sounding - rlapse * (alt - oldalt)
!        else
!          sounding = sounding + rlapse/2. * (alt - oldalt)
!        endif  
!        
!        do ixy=1,NXY
!          t2(ixy, k) = sounding
!        enddo
!        
!        oldalt = alt
        
c        write(LUNOPRT,*) k, ' ', alt, ' ', zc2(1, k), ' ', t2(1, k)
!      enddo

c
c  Initialize 3-D air density [g/cm^3]
c   (For this demo, use ideal gas law equation of state, p=rho*R*T)
c
!      do k = 1, NZ
!        do ixy=1,NXY
!          rhoa2(ixy,k) = p2(ixy,k) / ( R_AIR * t2(ixy,k) )
!        enddo
!      enddo
c
c
c  Initialize 3-D potential temperature concentration [g*K/cm^3]
c
!      do k=1,NZ
!        do ixy=1,NXY
!          pt_xy = t2(ixy,k) * ( PREF / p2(ixy,k) )**RKAPPA
!          ptc2(ixy,k) = rhoa2(ixy,k) * pt_xy
!        enddo
!      enddo
c
c
c  Define <zmet>, the  metric factor that makes <dz> in
c  an arbitrary vertical coord system a true distance.
c  See Toon et al., JAS, vol 45, #15, Aug-1988, p 2125.
c  
c  For hybid coordinates, eta (n) goes from 1.0 at the surface to 0.0
c  at infinity. The hybrid coordinate is defined in Simmions and Struffing,
c  ECMWF Tech. Report #28, p 32, 1981.
c
c    eta = a(k) + b(k)
c    
c  <zmet> is dh/deta = dh/dp * dp/deta
c  using dp=rho*g*dh, dh/deta = (1/(rho*g)) * (dp/deta)
c      write(LUNOPRT,*) 'initatm: k, p, pl, zl, zmet'
!      do k = 1, NZ
!        do ixy = 1, NXY
!          zmet2(ixy,k) = ((pl2(ixy, k) - pl2(ixy, k+1)) / 
!     $                    (zl2(ixy, k) - zl2(ixy, k+1))) /
!     $                   (GRAV * rhoa2(ixy,k))
!        enddo
c      write(LUNOPRT,*) k, ' ', p2(1, k), ' ',
c     $   pl2(1, k), ' ', zl2(1, k), ' ', zmet2(1, k)
!      enddo
c
c
c  Do vert integrated consistency check of rhoa, dz, zmet for troposphere
c
c      if (do_print_setup) then
c        ixy = 1
c        sum_h = 0.d0
c        sum_pres = 0.d0
c        do k=1,NZ
c         sum_h = sum_h + abs(dz2(ixy,k)) * zmet2(ixy,k)
c         sum_pres = sum_pres +
c     $     GRAV * rhoa2(ixy,k) * abs(dz2(ixy,k)) * zmet2(ixy,k)
c        enddo
c
c        sum_h = sum_h / 1.d+5    ! [km]
c        write(LUNOPRT,*) ' '
c        write(LUNOPRT,*) 'g_ctm_ll_hyb: sum_h = ', sum_h
c        write(LUNOPRT,*) 'g_ctm_ll_hyb: sum_pres = ', sum_pres
c        write(LUNOPRT,*) 'g_ctm_ll_hyb: p_surf = ', p_surf2(ixy)
c        write(LUNOPRT,*) 'g_ctm_ll_hyb: p_top = ', p_top2(ixy)
c        write(LUNOPRT,*) 'g_ctm_ll_hyb: p_top + sum_pres = ',
c     $                                       p_top2(ixy) + sum_pres
c        write(LUNOPRT,*) ' '
c      endif
c
c
c  Compute constants used in viscosity 
c
      rmu_0 = 1.8325e-4
      rmu_t0 = 296.16
      rmu_c = 120.
      rmu_const = rmu_0 * (rmu_t0 + rmu_c)
c
c
c  Define vertical profile of atmospheric variables that are 1-D.
c   In this demo:
c    Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
c     Meteorological Tables, in which there is a misprint -- T is deg_K, not deg_C.
c    Thermal conductivity of dry air <thcond> is from Pruppacher and Klett, Eq. 13-16. - EJL multiplied by 10 to convert to cgs
c
      do k=1,NZ
        rmu(k) = 10.0*rmu_const / ( t(1,1,k) + rmu_c ) *
     $              ( t(1,1,k) / rmu_t0 )**1.5d0
        thcond(k) = ( 5.69 + .017*( t(1,1,k) - T0 ) )*4.18e2
      enddo
c
c
c  Define vertical wind speed & diffusion coefficients [in metric cgs units]
c  This needs to be changed for Titan - EJL
      do k=1,NZP1
        do ix=1,NX
          do iy=1,NY
            w(ix,iy,k) =  0.
            dkx(ix,iy,k) = 0.e5
            dky(ix,iy,k) = 0.e5
            dkz(ix,iy,k) = 0.e-1

! From E.Barth - EJL        
c         Constant diffusion coefficient at/below 90 km
!            if( zl(ix,iy,k) .le. 90.0d5) dkz(ix,iy,k) = 5.d3

!            rhoa150 = 4.700000d-06
!            if( zl(ix,iy,k) .eq. 90.0d5)
!     $          y = (dlog10(2.5d13/sqrt(AVG * rhoa150/WTMOL_AIR))
!     $                  - dlog10(5.d3))/(60.d5/dz(ix,iy,k))

!            if((zl(ix,iy,k) .gt. 90.0d5) .and.
!     $                          (zl(ix,iy,k) .lt. 150.0d5)) then
!              dkz(ix,iy,k) = 10**(dlog10(5.d3) + y)
!              y = (dlog10(2.5d13/sqrt(AVG * rhoa150/WTMOL_AIR))
!     $                  - dlog10(5.d3))/(60.d5/dzfix) + y
!            endif

!            if((zl(ix,iy,k) .ge. 150.0d5) .and.
!     $                          (zl(ix,iy,k) .lt. 400.0d5))
!     $        dkz(ix,iy,k) = 2.5d13/sqrt(AVG * rhoa(ix,iy,k)/WTMOL_AIR)


c         Constant diffusion coefficient at/above 400 km
!            if( zl(ix,iy,k) .ge. 400.0d5) dkz(ix,iy,k) = 1.d8

c         Vertical wind between 320 and 430 km
cc          if( zl(ix,iy,k) .ge. 320.0d5 .and.
cc   $                        zl(ix,iy,k) .lt. 430.0d5)
cc   $        w(ix,iy,k) = 1.25
!          if(k.eq.32) write(*,*) '<initatm> No wind'

cc        Sensitivity Test:  Change diffusion coefficient by constant factor
cc          dkz(ix,iy,k) = 2.*dkz(ix,iy,k)
cc          if( ix.eq.1 .and. iy.eq.1 .and. k.eq.1 )
cc   $       write(*,*) 'Changing eddy diffusion by constant factor !!!'


          enddo
        enddo
      enddo
c
c
c  Define horizontal wind speed [in metric cgs units]
c
      do k=1,NZ
        do ix=1,NX
          do iy=1,NY
            u(ix,iy,k) = 0.e4
            v(ix,iy,k) = 0.
          enddo
        enddo
      enddo
c
c
c  Return to caller with atm profiles initialized
c
      return
      end

