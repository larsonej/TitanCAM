      subroutine charge
c
c
c  @(#) charge.f  Barth  Jul-1999
c  This routine calculates charge to radius ratio for determining
c  sticking coefficient
c
c  Argument list input:
c    k (altitude)
c    prrat (charge to radius ratio return value)
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
c  Local declarations
c
c    <emass> are ion masses [amu]
c    <erate> are ionization rates due to galactic cosmic rays
c    <erate2> are ionization rates due to precipitating electrons
c		   ionization rate units are [ion pairs/cm**3/s]

       dimension emass(61), erate(61), erate2(61)
c    $		 em(NZ),    q(NZ),     q2(NZ)

c    Ion masses and ion rates from Borucki et al. (Icarus 72, 604-622, 1987)
c    Values are defined for 10km layers starting at 10km.
c    (Surface value for erate was interpolated)

      data emass/ 1.11e+02,
     $ 1.11e+02, 1.11e+02, 1.11e+02, 1.11e+02, 1.08e+02, 1.06e+02,
     $ 1.01e+02, 9.60e+01, 8.90e+01, 8.30e+01, 7.80e+01, 7.30e+01,
     $ 6.80e+01, 6.30e+01, 5.80e+01, 5.40e+01, 5.80e+01, 6.30e+01,
     $ 6.90e+01, 7.80e+01, 8.80e+01, 9.40e+01, 1.01e+02, 1.09e+02,
     $ 1.09e+02, 1.09e+02, 1.09e+02, 1.09e+02, 1.09e+02, 1.09e+02,
     $ 1.09e+02, 1.08e+02, 1.06e+02, 1.04e+02, 1.02e+02, 1.01e+02,
     $ 9.90e+01, 9.70e+01, 9.50e+01, 9.30e+01, 9.10e+02, 8.90e+01,
     $ 8.70e+01, 8.50e+01, 8.30e+01, 8.00e+01, 8.00e+01, 8.00e+01,
     $ 8.00e+01, 8.00e+01, 8.00e+01, 8.00e+01, 8.00e+01, 8.00e+01,
     $ 8.00e+01, 8.00e+01, 8.00e+01, 8.00e+01, 8.00e+01, 8.00e+01 /

      data erate/ 2.48e-01,
     $ 3.14e-01, 3.80e-01, 4.47e-01, 5.30e-01, 6.70e-01, 1.38e+00, 
     $ 3.76e+00, 8.94e+00, 1.16e+01, 1.22e+01, 1.08e+01, 8.47e+00, 
     $ 6.27e+00, 4.62e+00, 3.45e+00, 2.59e+00, 1.93e+00, 1.49e+00, 
     $ 1.14e+00, 8.94e-01, 6.90e-01, 5.33e-01, 4.15e-01, 3.21e-01, 
     $ 2.59e-01, 1.96e-01, 1.15e-01, 1.12e-01, 9.35e-02, 7.29e-02,
     $ 5.68e-02, 4.43e-02, 3.45e-02, 2.69e-02, 2.10e-02, 1.64e-02,
     $ 1.28e-02, 9.94e-03, 7.75e-03, 6.00e-03, 5.00e-03, 4.00e-03,
     $ 3.00e-03, 2.00e-03, 1.80e-03, 1.30e-03, 1.00e-03, 8.00e-04, 
     $ 6.00e-04, 5.00e-04, 4.00e-04, 3.00e-04, 2.50e-04, 1.80e-04, 
     $ 1.30e-04, 1.00e-04, 8.00e-05, 6.00e-05, 5.00e-05, 3.50e-05 /

      data erate2/ 0.00e+00,
     $ 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     $ 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     $ 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     $ 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     $ 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     $ 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 1.00e-03,
     $ 3.00e-03, 5.00e-03, 7.00e-03, 9.00e-03, 1.00e-02, 2.00e-02,
     $ 4.00e-02, 5.00e-02, 7.00e-02, 9.00e-02, 1.00e-01, 2.00e-01,
     $ 4.00e-01, 5.00e-01, 7.00e-01, 9.00e-01, 1.00e+00, 2.00e+00, 
     $ 3.00e+00, 4.00e+00, 6.00e+00, 8.00e+00, 9.00e+00, 1.00e+01 /
c
c
c  Define formats
c
    1 format(i4,4e13.6)
    2 format(i4,'pc:  ',3e13.6)
    3 format('nr2 sum: ',2e13.6)
    4 format(i4,' pr ',e13.6)
    5 format(i4,' diff ',2e13.6)
    6 format('<charge> 'i4,e13.6)
    7 format(i4,4e14.6)
c
c-----------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter charge'
c
c-----------------------------------------------------------------------------
c
c  electron charge [esu]
      e = -4.8d-10

c
      onex = -ONE*exp(ONE)

c  Charging of particles below 300km

      do k = 1,NZ

      if(zl(1,1,k) .gt. 300.d5) return

c  Interpolate for layer thicknesses other than 10km
c
       modval = mod( zl(1,1,k), 10.d5 )
       i = int( zl(1,1,k) / 10.d5 ) + 1

       if( modval .eq. 0 ) then
         em = emass(i)
         q  = erate(i)
         q2 = erate2(i)
       else
         deltay = dble(modval)
         em = ( emass(i+1)-emass(i) )/10.d5 * deltay + emass(i)
         q  = ( erate(i+1)-erate(i) )/10.d5 * deltay + erate(i)
         q2 = ( erate2(i+1)-erate2(i) )/10.d5 * deltay + erate2(i)
       endif

c 	thermal velocity
        vth = sqrt(8.0 * BK * t(1,1,k)/ em * AVG)

             Sum_nc = 0.
             Sum_nr2 = 0.
             do jbin = 1,NBIN

                Sum_nc = Sum_nc + pc3(k,jbin,1) /
     $                                (e*e)*rm(jbin,1,k)*BK*t(1,1,k)
                Sum_nr2 = Sum_nr2 + pc3(k,jbin,1) * rm(jbin,1,k)**2
             enddo

        xemx = -(q + q2)/(PI * vth * Sum_nc * Sum_nr2)
c       write(*,1) k,xemx,vth,Sum_nc,Sum_nr2

c  Solve for x

        x = -10.
        do while (xemx .lt. onex) 		!x < -1
          if(xemx .lt. x*exp(-x)) then
             x = x - 10.
          else
            do while(xemx .gt. x*exp(-x)) 
                x = x + 1.
            end do
                 xinc = 1.0d-1
               iflag = -1
         goto 1000
           endif
          end do

c  For (-1. <= x < 0.) Taylor expand x*exp(-x) and use the quadratic formula
c  to solve for x

          if(xemx .ge. onex) then
             x = 0.5* (1. - sqrt(1. - 4.*xemx) )
             xinc = 1.0d-1
             iflag = 1 
           endif

 1000         dif = x*exp(-x) / xemx
              if(dif .lt. 0.999) then
                    if(iflag .eq. 1) xinc = xinc/10.
                     x = x - xinc
                    iflag = -1 
                    else if(dif .gt. 1.001) then
                     if(iflag .eq. -1) xinc = xinc/10.
                        x = x + xinc
                      iflag = 1	
                else
                  goto 9000
                endif
                 goto 1000

c  Calculate charge to radius ratio

 9000    prrat(k) = 3.0 * x * BK*t(1,1,k)/ (e*e) 
c       write(*,7) k,Sum_nc,Sum_nr2,x,prrat(k)       
        enddo

!        open(unit=99, file='carma_charge.txt', status='unknown')
!        write(99,*) 'k, t(1,1,k), e, Sum_nc, Sum_nr2, x, prrat(k)'
!        do k=1,NZ
!          write(99,*) k, t(1,1,k), e, Sum_nc, Sum_nr2,x,prrat(k)
!          write(99,*) ''
!         enddo
!        close(unit=99)
c
c
c  Return to caller with ratio calculated
c
      return
      end