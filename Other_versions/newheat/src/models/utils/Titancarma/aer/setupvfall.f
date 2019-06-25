       subroutine setupvfall
c
c
c  @(#) setupvfall.f  Ackerman Oct-1995
c
c  This routine evaluates particle fall velocities, vf(k) [cm s^-1]
c  and reynolds' numbers based on fall velocities, re(j,i,k) [dimensionless].
c  indices correspond to vertical level <k>, bin index <i>, and aerosol
c  group <j>.
c
c  Non-spherical particles are treated through shape factors <ishape>
c  and <eshape>.
c  
c  General method is to first use Stokes' flow to estimate fall
c  velocity, then calculate reynolds' number, then use "y function" 
c  (defined in Pruppacher and Klett) to reevaluate reynolds' number,
c  from which the fall velocity is finally obtained.
c
c  This routine requires that vertical profiles of temperature <t>,
c  air density <rhoa>, and viscosity <rmu> are defined (i.e., initatm.f
c  must be called before this).  The vertical profile with ix = iy = 1
c  is used.
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
    1 format('In subroutine setupvfall ishape != 1 :',
     $       ' no fall velocity algorithm')
    2 format(/,'Fall velocities and Reynolds'' number in bottom layer',
     $  ' (setupvfall): ')
    3 format(/,'Particle group ',i3,/,' bin   r [cm]       vf [cm/s]',
     $  '         re'/)
    4 format(i3,3(1pe11.3,4x)) 
    6 format(i3,f7.1,i6,1pe13.4)
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupvfall'
c 
c-------------------------------------------------------------------------------
c
c  Loop over aerosol groups.
c
      do j = 1,NGROUP
c
c  First evaluate factors that depend upon particle shape (used in correction
c  factor <bpm> below).  
c
        if( ishape(j) .eq. 1 )then
c
c Spheres
c
          f1 = 1.0
          f2 = 1.0

        else if( ishape(j) .eq. 2 )then
c
c  Hexagons: taken from Turco et al
c  (Planet. Space Sci. Rev. 30, 1147-1181, 1982)
c  with diffuse reflection of air molecules assumed 
c
          f2 = (PI/9./tan(PI/6.))**(ONE/3.)*eshape(j)**(ONE/6.)

        else if( ishape(j) .eq. 3 )then
c
c  Spheroids: taken FROM WHERE?
c
          f2 = (ONE*2./3.)**(ONE/3.)*eshape(j)**(ONE/6.)

        endif
c
c  (following statement yields <f3> = 1.0 for <eshape> = 1)
c
        f3 = 1.39/sqrt((1.14+0.25/eshape(j))*(0.89+eshape(j)/2.))
        f2 = f2*f3

        if( eshape(j) .gt. 1. )then
c
c  For Stokes regime there is no separate data for hexagonal plates or columns,
c  so we use prolate spheroids.  This is from Fuchs' book.
c
          exx = eshape(j)**2 - 1.
          exy = sqrt(exx)
          xcc = 1.333*exx/((2.*eshape(j)**2-1.)*log(eshape(j)+exy)/exy
     $         -eshape(j))
          xa = 2.666*exx/((2.*eshape(j)**2-3.)*log(eshape(j)+exy)/exy
     $         +eshape(j))
          f1 = eshape(j)**(-ONE/3.)*(xcc+2.*xa)/3.

        elseif( eshape(j) .lt. 1. )then
c
c  Use oblate spheroids for disks (eshape < 1.).  Also from Fuchs' book.
c
          bxx = 1./eshape(j)
          exx = bxx**2 - 1.
          exy = sqrt(exx)
          xcc = 1.333*exx/(bxx*(bxx**2-2.)*atan(exy)/exy+bxx)
          xa = 2.666*exx/(bxx*(3.*bxx**2-2.)*atan(exy)/exy-bxx)
          f1 = bxx**(ONE/3.)*(xcc+2.*xa)/3.

        endif
c
c  Loop over column with ixy = 1
c
        ixy = 1
        do k = 1,NZ
c
c  This is <rhoa> in cartesian coordinates.
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
c  <r_shape> is radius of particle used to calculate <re>.
c
            if( ishape(j) .eq. 1 )then

              r_shape = r(i,j)

            else if( ishape(j) .eq. 2 )then

              r_shape = r(i,j)*0.85*eshape(j)**(-ONE/3.)

            else if( ishape(j) .eq. 3 )then

              r_shape = r(i,j)*eshape(j)**(-ONE/3.)

            endif
c
c  <rkn> is knudsen number
c
            rkn = rmfp/r(i,j)
c
c  <bpm> is correction term for particle shape and non-continuum effects.
c  It is also used to calculate coagulation kernels and diffusion coefficients.
c
            expon = -.87 / rkn
            expon = max(-POWMAX, expon)
            bpm(k,i,j) = 1. + f1*f2*( 1.246*rkn + 0.42*rkn*exp(expon) )
c
c
c  These are first guesses for fall velocity and Reynolds' number, 
c  valid for Reynolds' number < 0.01
c
            vf(k,i,j) = (ONE*2./9.)*rhop2(ixy,k,i,j)*(r(i,j)**2)
     $                  *GRAV*bpm(k,i,j)/(f1*rmu(k))
 
            re(k,i,j) = 2. * rhoa_cgs * r_shape *
     $                  vf(k,i,j)/rmu(k)

c           write(LUNOTEMP,6) j,zl3(k)/1.d5,i,re(k,i,j)
c
c
c  <rfix> is used in drag coefficient.
c
            rfix = vol(i,j) * rhop2(ixy,k,i,j) * GRAV *
     $             rhoa_cgs / rmu(k)**2
 
            if( (re(k,i,j) .ge. 0.01) .and. (re(k,i,j) .le. 300.) )then
c
c
c  This is "regime 2" in Pruppacher and Klett (chap. 10).
c
              if( ishape(j) .eq. 1 )then

                x = log(24.*re(k,i,j)/bpm(k,i,j))
                y = -3.18657 + x*(0.992696    - x*(.00153193
     $                       + x*(0.000987059 + x*(.000578878
     $                       - x*(8.55176E-05 - x* 3.27815E-06 )))))
                if( y .lt. -675. ) y = -675.
                if( y .ge.  741. ) y =  741.

                re(k,i,j) = exp(y)*bpm(k,i,j)

              else if( eshape(j) .le. 1. )then

                if( ishape(j) .eq. 2 )then
                  x = log10(16.*rfix/(3.*sqrt(3.*ONE)))
                else if( ishape(j) .eq. 3 )then
                  x = log10(8.*rfix/PI)
                endif

                if( eshape(j) .le. 0.2 )then
                  b0 = -1.33
                  bb1 = 1.0217
                  bb2 = -0.049018
                  bb3 = 0.0
                else if( eshape(j) .le. 0.5 )then
                  ex = (eshape(j)-0.2)/0.3
                  b0 = -1.33+ex*(-1.3247+1.33*ONE)
                  bb1 = 1.0217+ex*(1.0396-1.0217*ONE)
                  bb2 = -0.049018+ex*(-0.047556+0.049018*ONE)
                  bb3 = 0.0+ex*(-0.002327*ONE)
                else
                  ex = (eshape(j)-0.5)/0.5
                  b0 = -1.3247+ex*(-1.310+1.3247*ONE)
                  bb1 = 1.0396+ex*(0.98968-1.0396*ONE)
                  bb2 = -0.047556+ex*(-0.042379+0.047556*ONE)
                  bb3 = -0.002327+ex*(0.002327*ONE)
                endif
                y = b0+x*(bb1+x*(bb2+x*bb3))

                re(k,i,j) = 10.**y*bpm(k,i,j)

              else if( eshape(j) .gt. 1. )then

                x = log10(2.*rfix/eshape(j))

                if( eshape(j) .le. 2 )then
                  ex = eshape(j)-1.
                  b0 = -1.310+ex*(-1.11812+1.310*ONE)
                  bb1 = 0.98968+ex*(0.97084-0.98968*ONE)
                  bb2 = -0.042379+ex*(-0.058810+0.042379*ONE)
                  bb3 = 0.0+ex*(0.002159*ONE)
                else if( eshape(j) .le. 10.)then
                  ex = (eshape(j)-2.)/8.0
                  b0 = -1.11812+ex*(-0.90629+1.11812*ONE)
                  bb1 = 0.97084+ex*(0.90412-0.97084*ONE)
                  bb2 = -0.058810+ex*(-0.059312+0.058810*ONE)
                  bb3 = 0.002159+ex*(0.0029941-0.002159*ONE)
                else
                  ex = 10./eshape(j)
                  b0 = -0.79888+ex*(-0.90629+0.79888*ONE)
                  bb1 = 0.80817+ex*(0.90412-0.80817*ONE)
                  bb2 = -0.030528+ex*(-0.059312+0.030528*ONE)
                  bb3 = 0.0+ex*(0.0029941*ONE)
                endif
                y = b0+x*(bb1+x*(bb2+x*bb3))

                re(k,i,j) = 10.**y*bpm(k,i,j)

              endif
c
c
c  Adjust <vf> for non-sphericicity.
c
              vf(k,i,j) = re(k,i,j) * rmu(k) /
     $                    (2.*r_shape*rhoa_cgs)

            endif
 
            if( re(k,i,j) .gt. 300. )then

              if( ishape(j) .ne. 1 ) write (LUNOPRT,1)

              z  = ((1.e6*rhoa_cgs**2)
     $           /  (GRAV*rhop2(ixy,k,i,j)*rmu(k)**4))**(ONE/6.)
              b0 = (24.*vf(k,i,j)*rmu(k))/100.
              x  = log(z*b0)
              y  = -5.00015 + x*(5.23778   - x*(2.04914 - x*(0.475294
     $                      - x*(0.0542819 - x* 0.00238449 ))))
              if( y .lt. -675. )  y = -675.0
              if( y .ge.  741. )  y =  741.0

              re(k,i,j) = z*exp(y)*bpm(k,i,j)
 
              vf(k,i,j) = re(k,i,j) * rmu(k) /
     $                         ( 2. * r(i,j) * rhoa_cgs )

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
      k = 1
      
      write(LUNOPRT,2)

      do j = 1,NGROUP
        
        write(LUNOPRT,3) j

        do i = 1,NBIN
          write(LUNOPRT,4) i,r(i,j),vf(k,i,j),re(k,i,j)
        enddo
      enddo
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
              vf(k,ibin,ig) = -vf(k,ibin,ig) / zmet_k
            enddo
          enddo

         enddo

      endif
c
c
c  Write all fall velocites for mass flux calculation
c
      do j = 1,NGROUP
        do k = 1,NZ
          do i = 1,NBIN 
 
            write(LUNOTEMP,6) j,zl3(k)/1.d5,i,vf(k,i,j)
 
          enddo
        enddo
      enddo
      close(LUNOTEMP)
c
c  Return to caller with particle fall velocities evaluated.
c
      return
      end
