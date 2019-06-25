       subroutine setupvf
c
c
c  @(#) setupvf.f  Ackerman Nov-2000
c
c  This routine evaluates particle fall velocities, vf(k) [cm s^-1]
c  and reynolds' numbers based on fall velocities, re(j,i,k) [dimensionless].
c  indices correspond to vertical level <k>, bin index <i>, and aerosol
c  group <j>.
c
c  Method: first use Stokes flow (with Fuchs' size corrections, 
c  valid only for Stokes flow) to estimate fall velocity, then calculate
c  Reynolds' number (Re) (for spheres, Stokes drag coefficient is 24/Re).
c  Then for Re > 1, correct drag coefficient (Cd) for turbulent boundary
c  layer through standard trick to solving the drag problem: 
c  fit y = log( Re ) as a function of x = log( Cd Re^2 ).  
c  We use the data for rigid spheres taken from Figure 10-6 of
c  Pruppacher and Klett (1978):
c
c   Re     Cd
c  -----  ------
c     1    24
c    10     4.3
c   100     1.1
c  1000     0.45
c
c  Note that we ignore the "drag crisis" at Re > 200,000
c  (as discussed on p. 341 and shown in Fig 10-36 of P&K 1978), where
c  Cd drops dramatically to 0.2 for smooth, rigid spheres, and instead 
c  assume Cd = 0.45 for Re > 1,000
c
c  Note that we also ignore hydrodynamic deformation of liquid droplets
c  as well as any breakup due to Rayleigh-Taylor instability.  
c
c  This routine requires that vertical profiles of temperature <t>,
c  air density <rhoa>, and viscosity <rmu> are defined (i.e., initatm.f
c  must be called before this).  The vertical profile with ix = iy = 1
c  is used.
c
c  We assume spherical particles -- call setupvf_old() to use legacy
c  code from old Toon model for non-spherical effects -- use (better
c  yet, fix) at own risk.
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
c  Define formats
c
    1 format(/,'Non-spherical particles specified for group ',i3,
     $  ' (ishape=',i3,') but spheres assumed in setupvf.f.', 
     $  ' Suggest using/fixing non-spherical code in setupvf_old.f.')
    2 format(/,'Fall velocities and Reynolds'' number in bottom layer',
     $  ' (setupvf): ')
    3 format(/,'Particle group ',i3,/,' bin   r [cm]       dr [cm]',
     $  '        mass [g]       vf [cm/s]         re'/)
    4 format(i3,5(1pe11.3,4x)) 
    5 format(/,'Particle fall velocities for ',a,'(setupvf)',//,
     $  a3, 1x, a11, 4x, a25 /)
    6 format(i4,2x,f6.2,35(2x,1pe13.6))
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupvf'
c 
c-------------------------------------------------------------------------------
c
c  Loop over aerosol groups.
c
      do j = 1,NGROUP
c
c  Warning message for non-spherical particles
c
        if( ishape(j) .ne. 1 )then
          write(*,1) j, ishape
        endif
c
c  Loop over column with ixy = 1
c
        ixy = 1
        do k = 1,NZ
c
c  This is <rhoa> in cartesian coordinates (good old cgs units)
c
          rhoa_cgs = rhoa2(ixy,k) /
     $               (xmet2(ixy,k)*ymet2(ixy,k)*zmet2(ixy,k))
c
c  <vg> is mean thermal velocity of air molecules [cm/s]
c
          vg = sqrt(8./PI * R_AIR*t2(ixy,k))
c
c  <rmfp> is mean free path of air molecules [cm]
c
          rmfp = 2.*rmu(k) / (rhoa_cgs*vg)
c
c  Loop over particle size bins.
c
          do i = 1,NBIN
c
c  <rkn> is knudsen number
c
            rkn = rmfp/r(i,j)
c
c  <bpm> is correction term for non-continuum effects.  Also used to 
c  calculate coagulation kernels and diffusion coefficients.
c
            expon = -.87 / rkn
            expon = max(-POWMAX, expon)
            bpm(k,i,j) = 1. + ( 1.246*rkn + 0.42*rkn*exp(expon) )

c             print *,zl3(k)/1.d5,r(i,j),bpm(k,i,j)
c
c
c  Stokes fall velocity and Reynolds' number
c
            vf(k,i,j) = (ONE*2./9.)*rhop2(ixy,k,i,j)*r(i,j)**2
     $                  *GRAV*bpm(k,i,j)/rmu(k)
 
            re(k,i,j) = 2.*rhoa_cgs*r(i,j)*vf(k,i,j)/rmu(k)
 
            if( re(k,i,j) .ge. 1. )then
c
c   Correct drag coefficient for turbulence 
c
              x = log( re(k,i,j)/bpm(k,i,j) )
              y = x*(0.83 - 0.013*x)

              re(k,i,j) = exp(y)*bpm(k,i,j)

              if( re(k,i,j) .le. 1.e3 )then
c  
c  drag coefficient from quadratic fit y(x) when Re < 1,000
c
                vf(k,i,j) = re(k,i,j) * rmu(k) /
     $                      (2.*r(i,j)*rhoa_cgs)
              else
c
c  drag coefficient = 0.45 independent of Reynolds number when Re > 1,000
c
                cdrag = 0.45 
                vf(k,i,j) = bpm(k,i,j)*
     $                      sqrt( 8.*rhop2(ixy,k,i,j)*r(i,j)*GRAV / 
     $                      (3.*cdrag*rhoa_cgs) )
              endif

            endif
c
c
c  Interpolate <vf> from layer mid-pts to layer boundaries.
c  <vf(k)> is the fall velocity at the lower edge of the layer
c
          if( k .gt. 1 .and. k .lt. NZ )
     $      vf(k,i,j) = sqrt( vf(k-1,i,j)*vf(k,i,j) )

          enddo    ! <i=1,NBIN>

        enddo      ! <k=1,NZ>

        do i = 1,NBIN
          vf(NZ+1,i,j) = vf(NZ,i,j)
          if( NZ .gt. 1 )then
            vf(NZ,i,j) = sqrt( vf(NZ-1,i,j)*vf(NZ,i,j) )
          endif
        enddo
          
      enddo        ! <j=1,NGROUP>
c
c
c  Constant value if <ifall> = 0
c
      if( ifall .eq. 0 )then
        do j = 1,NGROUP
          do i = 1,NBIN
            do k = 1,NZP1
              vf(k,i,j) = vf_const
            enddo
          enddo
        enddo
      endif
c
c
c  Print out fall velocities and reynolds' numbers in lowest model layer.
c
      if (do_print_setup) then
      k = NZP1
      
      write(LUNOPRT,2)

      do j = 1,NGROUP
        
        write(LUNOPRT,3) j

        do i = 1,NBIN
          write(LUNOPRT,4) i,r(i,j),dr(i,j),rmass(i,j),
     $                                  vf(k,i,j),re(k,i,j)
        enddo
      enddo
      endif
c
c
c  Print out all the fall velocities as a table.
c
      open(unit=99, file='carma_vf.txt', status='unknown')
      write(99,*) 'Fall Velocities (cm/s): level, bin'

      do j=1,NGROUP
      
        write(LUNOPRT,3) j

        do k=1,NZP1
          write(99,*) vf(k,:,j)
          write(99,*) ''
        enddo

        write(99,*) ''
      enddo

      close(unit=99)
c
c
c  Scale cartesian fallspeeds to the appropriate vertical coordinate system.
c  Non--cartesion coordinates are assumed to be positive downward, but
c  vertical velocities in this model are always assumed to be positive upward. 
c
      if( igridv .ne. I_CART )then

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

          a_mid_k1 = zc2(ixy,k1)
          a_mid_k2 = zc2(ixy,k2)

          if(  a_mid_k2 .ne. a_mid_k1 )then
            frac = ( zl2(ixy,k) - a_mid_k1 ) / ( a_mid_k2 - a_mid_k1 )
          else
            frac = 0.
          endif

          zmet_k1 = zmet2(1,k1)
          zmet_k2 = zmet2(1,k2)
          zmet_k = zmet_k1 + frac * ( zmet_k2 - zmet_k1 )

          do ig = 1,NGROUP
            do ibin = 1,NBIN
              vf(k,ibin,ig) = vf(k,ibin,ig) / zmet_k
            enddo
          enddo

         enddo

      endif
c
c
c  Write all fall velocites for mass flux calculation
c
      do j = 1,NGROUP
        write(LUNOPRT,5) groupname(j),
     $              'iz','zl [km]','vfall [cm/s], all bins'
        do k = 1,NZ
          write(LUNOPRT,6) k,zl3(k)/1.d5,(vf(k,i,j),i=1,NBIN)
        enddo
      enddo
c
c  Return to caller with particle fall velocities evaluated.
c
      return
      end
