      subroutine freezaerl
c
c
c  @(#) freezaerl.f  Jensen  Mar-1995
c  This routine evaluates particle loss rates due to nucleation <rnuclg>:
c  aerosol freezing only.
c  The parameterization described by Tabazadeh et al. [GRL, 27, 1111, 2000.] is
c  used.
c  
c  The loss rates for all particle elements in a particle group are equal.
c
c  To avoid nucleation into an evaporating bin, this subroutine must
c  be called after growp, which evaluates evaporation loss rates <evaplg>.
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c
c  Argument list input:
c
c  Argument list output:
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c
      logical evapfrom_nucto
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter freezaerl'
c
c
c  Loop over particle groups.
c
      do igroup = 1,NGROUP
       do igs = 1,NGAS

       igas = inucgas(igs,igroup)                ! condensing gas
       iepart = ienconc( igroup )            ! particle number density element
 
       if( igas .ne. 0 )then
c
c
c  Calculate nucleation loss rates.  Do not allow nucleation into
c  an evaporating bin.
c
c  <ienucto> is index of target nucleation element;
c  <ignucto> is index of target nucleation group.
c
c        if( nnuc2elem(iepart) .gt. 1 )then
        do inuc = 1,nnuc2elem(iepart)

         ienucto = inuc2elem(inuc,iepart)
         if( ienucto .ne. 0 )then
           ignucto = igelem( ienucto )
         else
           ignucto = 0
         endif
c
c
c  Only compute nucleation rate for aerosol freezing
c
         if( inucproc(iepart,ienucto) .eq. I_AERFREEZE ) then
c
c
c  Loop over particle bins.  Loop from largest to smallest for 
c  evaluation of index of smallest bin nucleated during time step <inucstep>.
c
          do ibin =NBIN,1,-1
c
c
c  <inucto> is index of target nucleation bin.
c
           if( ignucto .ne. 0 )then
             inucto = inuc2bin(ibin,igroup,ignucto)
           else
             inucto = 0
           endif
c
c
c  Bypass calculation if few particles are present
c
           if( pconmax(ixyz,igroup) .gt. FEW_PC )then
c
c
c  Set <evapfrom_nucto> to .true. when target droplets are evaporating
c  
            if( inucto .ne. 0 )then
              evapfrom_nucto = evaplg(inucto,ignucto) .gt. 0.
            else
              evapfrom_nucto = .false.
            endif
c
c
c  Calculate approximate critical saturation needed for homogeneous freezing
c  of sulfate aerosols (see Jensen and Toon, GRL, 1994).
c
            sifreeze = 0.3
c
c
c  Homogeneous freezing of sulfate aerosols should only occur of SL < Scrit
c  and SI > <sifreeze>.
c            
            if( supsati3(ixyz,igas) .gt. sifreeze .and.
     $          .not. evapfrom_nucto )then

             rlogs = log( supsatl3(ixyz,igas) + 1. )
c
c  Calculate mean ice density and latent heat of freezing over temperature
c  interval [T0,T]
c
             rhoibar = ( 0.916 * (t3(ixyz)-T0(igas)) -
     $         1.75e-4/2. * ((t3(ixyz)-T0(igas))**2) -
     $         5.e-7 * ((t3(ixyz)-T0(igas))**3)/3. ) / 
     $                                (t3(ixyz)-T0(igas))

             rlhbar = ( 79.7 * (t3(ixyz)-T0(igas)) +
     $                0.485/2. * (t3(ixyz)-T0(igas))**2 -
     $                2.5e-3/3. * (t3(ixyz)-T0(igas))**3 )
     $                / (t3(ixyz)-T0(igas)) * 4.186e7*18.
c
c
c  Equilibrium H2SO4 weight percent for fixed water activity
c
             act = supsatl3(ixyz,igas) + 1.
             IF(act .LT. 0.05) THEN
               CONTL = 12.37208932 * (act**(-0.16125516114)) -
     $                 30.490657554 * act  - 2.1133114241
               CONTH = 13.455394705 * (act**(-0.1921312255))  -
     $                 34.285174604 * act - 1.7620073078
             END IF
             IF(act .GE. 0.05 .and. act .LE. 0.85) THEN
               CONTL = 11.820654354 * (act**(-0.20786404244)) -
     $                 4.807306373 * act  - 5.1727540348
               CONTH = 12.891938068 * (act**(-0.23233847708)) -
     $                 6.4261237757 * act - 4.9005471319
             END IF
             IF(act .GT. 0.85) THEN
               CONTL = -180.06541028 * (act**(-0.38601102592)) -
     $                 93.317846778 * act + 273.88132245
               CONTH = -176.95814097 * (act**(-0.36257048154)) -
     $                 90.469744201 * act + 267.45509988
             END IF
             H2SO4m = CONTL + ((CONTH - CONTL) * (t3(ixyz) -190.)/70.)
             WT = (98.0 * H2SO4m)/(1000. + 98. * H2SO4m)
             WT = 100. * WT
c
c
c  Volume ratio of wet/dry aerosols.
c
             vrat = rhosol(1)/RHO_W * ((100.-wt)/wt) + 1.
c
c
c  Calculation sulfate solution density from Myhre et al. (1998).
c
             wtfrac = WT/100.
             C1      = t3(ixyz) - 273.15
             C2      = C1**2
             C3      = C1**3
             C4      = C1**4
             A0 = 999.8426 + 334.5402e-4*C1 - 569.1304e-5*C2
             A1 = 547.2659 - 530.0445e-2*C1 + 118.7671e-4*C2
     $            + 599.0008e-6*C3
             A2 = 526.295e+1 + 372.0445e-1*C1 + 120.1909e-3*C2
     $            - 414.8594e-5*C3 + 119.7973e-7*C4
             A3 = -621.3958e+2 - 287.7670*C1 - 406.4638e-3*C2
     $            + 111.9488e-4*C3 + 360.7768e-7*C4
             A4 = 409.0293e+3 + 127.0854e+1*C1 + 326.9710e-3*C2
     $            - 137.7435e-4*C3 - 263.3585e-7*C4
             A5 = -159.6989e+4 - 306.2836e+1*C1 + 136.6499e-3*C2
     $            + 637.3031e-5*C3
             A6 = 385.7411e+4 + 408.3717e+1*C1 - 192.7785e-3*C2
             A7 = -580.8064e+4 - 284.4401e+1*C1
             A8 = 530.1976e+4 + 809.1053*C1
             A9 = -268.2616e+4
             A10 = 576.4288e+3
             den = A0 + wtfrac*A1 + wtfrac**2 * A2 +
     $             wtfrac**3 * A3 + wtfrac**4 * A4
             den = den + wtfrac**5 * A5 + wtfrac**6 * A6 +
     $             wtfrac**7 * A7
             den = den + wtfrac**8 * A8 + wtfrac**9 * A9 +
     $             wtfrac**10 * A10
c
c
c  Activation energy is based on Koop's lab data.
c
             IF(t3(ixyz) .GT. 220) then
             A0      = 104525.93058
             A1      = -1103.7644651
             A2      = 1.070332702
             A3      = 0.017386254322
             A4      = -1.5506854268e-06
             A5      = -3.2661912497e-07
             A6      = 6.467954459e-10
             ELSE
             A0      = -17459.516183
             A1      = 458.45827551
             A2      = -4.8492831317
             A3      = 0.026003658878
             A4      = -7.1991577798e-05
             A5      = 8.9049094618e-08
             A6      = -2.4932257419e-11
             END IF
             diffact = ( A0 + A1*t3(ixyz) + A2*t3(ixyz)**2 +
     $           A3*t3(ixyz)**3 + A4*t3(ixyz)**4 +
     $           A5*t3(ixyz)**5 + A6*t3(ixyz)**6 ) * 1.0e-13
c
c
c  Surface energy
c
c  Weight percent function for T = 260 K
c
             c0      = 77.40682664
             c1      = -0.006963123274
             c2      = -0.009682499074
             c3      = 0.00088797988
             c4      = -2.384669516e-05
             c5      = 2.095358048e-07
             S260 = c0 + c1*wt + c2*wt**2 + c3*wt**3 +
     $              c4*wt**4 + c5*wt**5
c
c  Weight percent function for T = 220 K
c
             d0      = 82.01197792
             d1      = 0.5312072092
             d2      = -0.1050692123
             d3      = 0.005415260617
             d4      = -0.0001145573827
             d5      = 8.969257061e-07
             S220 = d0 + d1*wt + d2*wt**2 + d3*wt**3 +
     $              d4*wt**4 + d5*wt**5
c
c  Weight percent function for T = 180K
c
             e0      = 85.75507114
             e1      = 0.09541966318
             e2      = -0.1103647657
             e3      = 0.007485866933
             e4      = -0.0001912224154
             e5      = 1.736789787e-06
             S180 = e0 + e1*wt + e2*wt**2 + e3*wt**3 +
     $              e4*wt**4 + e5*wt**5

             if( t3(ixyz) .GE. 220. ) then
               sigma = S260 + ((260.-t3(ixyz))*(S220-S260))/40.
             else
               sigma = S220 + ((220.-t3(ixyz))*(S180-S220))/40.
             endif

             sigsula = sigma
             sigicea = 105.
             sigsulice = abs( sigsula - sigicea )
c
c
c  Critical ice germ radius formed in the sulfate solution
c
             ag = 2.*gwtmol(igas)*sigsulice /
     $            ( rlhbar * rhoibar * log(T0(igas)/t3(ixyz)) +
     $            rhoibar * rgas * 0.5 * (T0(igas)+t3(ixyz)) *
     $            log(supsatl3(ixyz,igas)+1.) )
             if( ag .lt. 0. ) ag = 1.e10
C
C     Gibbs free energy of ice germ formation in the ice/sulfate solution

             delfg = 4./3.*PI * Sigsulice * (ag**2)
C
C    Ice nucleation rate in a 0.2 micron aerosol (/sec)
C
             expon = ( -diffact - delfg ) / BK / t3(ixyz)
             expon = max( -100.*ONE, expon )
             rnuclg(ibin,igroup,ignucto) = prenuc *
     $              sqrt(Sigsulice*t3(ixyz)) *
     $              vrat*vol(ibin,igroup) * exp( expon )
c
c             xh = 0.1 * r(ibin,igroup) / ag
c             phih = sqrt( 1. - 2.*rmiv*xh + xh**2 )
c             rath = (xh-rmiv) / phih
c             fv3h = xh**3 * ( 2.*ONE - 3.*rath + rath**3 )
c             fv4h = 3. * rmiv * xh**2 * (rath-1.)
c             if( abs(rath) .gt. 1.e0-1.e-8 ) fv3h = 0.
c             if( abs(rath) .gt. 1.e0-1.e-10 ) fv4h = 0.
c 
c             fh = 0.5 * ( ONE + ((ONE-rmiv*xh)/phih)**3 +
c     $                              fv3h + fv4h )
c
c             expon = ( -delfwat2ice - delfg ) / BK / t3(ixyz)
c             expon = max( -POWMAX, expon )
c
              if(ibin .lt. inucstep(igroup)) then
               inucstep(igroup) = ibin
              endif

            endif

           endif   ! pconmax(ixyz,igroup) .gt. FEW_PC 
          enddo      ! ibin = 1,NBIN
         endif       ! inucproc(iepart,ienucto) .eq. I_DROPACT
c        endif       ! (nnuc2elem(iepart) .gt. 1)
        enddo        ! inuc = 1,nnuc2elem(iepart)
       endif        ! (igas = inucgas(igroup) .ne. 0)
       enddo        ! igs = 1,NGAS
      enddo         ! igroup = 1,NGROUP
c
c
c  Return to caller with particle loss rates due to nucleation evaluated.
c
      return
      end
