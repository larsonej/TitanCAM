      SUBROUTINE RADTRAN
C
C     **************************************************************
C     Purpose:    Driver routine for radiative transfer model.  
C
C     Input:      Temperature, vapor, and aerosol profiles and solar
C                 zenith angle are taken from interface common block.
C
C     Output:     Profiles of radiative fluxes, heating
C                 rates for air and particles; vertically 
C                 integrated optical depths; and albedos (which
C                 are all loaded into interface common block).
C     **************************************************************
C
      include 'globrad.h'
c
c  Reset flag for computation of solar fluxes
c
      ISL        = isl_aerad
c
c     Get atmospheric profiles from interface common block
c
      do k = 1,nvert
        t(k) = t_aerad(k)
        q(k) = qv_aerad(k)
      enddo
C
C     INTERPOLATE TEMPERATURE FROM LAYER CENTER (T) TO LAYER EDGE (TT)
C
      TT(1) = T(1)
      DO 12 J = 2, NVERT
         TT(J) = T(J-1) * (PRESS(J)/P(J-1)) **
     1             (log(T(J)/T(J-1))/log(P(J)/P(J-1)))
   12 CONTINUE
C
C     WATER VAPOR (G / CM**2)
C
      DO 10 J = 2, NLAYER
         RDH2O(J)   = Q(J-1) * DPG(J-1)
   10 CONTINUE
C
C     AEROSOL CONCENTRATIONS (# / CM**2)
C
      do ig = 1, NGROUP
        DO 15 J = 2, NVERT
           DO 15 I = 1, NRAD
              CAER(I,J,ig)  = pc_aerad(J-1,I,ig) 
   15   CONTINUE
      enddo
c
c     Solar zenith angle 
c
      u0 = u0_aerad
C
C     SURFACE REFLECTIVITY AND EMISSIVITY
C
      DO 20 L =  1,NSOLP
         RSFX(L) =  ALBEDO_SFC
         EMIS(L) =  0.0
 20   CONTINUE
C
      DO 30 L =  NSOLP+1,NTOTAL
         EMIS(L) =  EMISIR
         RSFX(L) = 1.0 - EMIS(L)
 30   CONTINUE
C
C     SET WAVELENGTH LIMITS LLA AND LLS BASED ON VALUES OF ISL AND IR
C
        LLA                   =  NTOTAL
        LLS                   =  1
        IF(ISL  .EQ. 0) THEN
          LLS   =  NSOLP+1
        ENDIF
C
        IF(IR   .EQ. 0) THEN
          LLA  =  NSOLP
        ENDIF
C
C     CALCULATE THE OPTICAL PROPERTIES
C
        CALL OPPR
C
C     IF INFRARED CALCULATIONS ARE REQUIRED THEN CALCULATE
C     THE PLANK FUNCTION
C
        IF(IR .NE. 0) THEN
           CALL OPPR1
        ENDIF
C
C     IF NO INFRARED SCATTERING THEN SET INDEX TO NUMBER OF
C     SOLAR INTERVALS
C
        IF(IRS .EQ. 0) THEN
          LLA  =  NSOLP
        ENDIF
C
C     IF EITHER SOLAR OR INFRARED SCATTERING CALCULATIONS ARE REQUIRED
C     CALL THE TWO STREAM CODE AND FIND THE SOLUTION
C
        IF(ISL .NE. 0 .OR. IRS. NE. 0 ) THEN
          CALL TWOSTR
          CALL ADD
        ENDIF
C
C     IF INFRARED CALCULATIONS ARE REQUIRED THEN CALL NEWFLUX1 FOR
C     A MORE ACCURATE SOLUTION
C
        IF(IR .NE. 0) THEN
          CALL NEWFLUX1
        ENDIF

C     CALCULATE INFRAFRED AND SOLAR HEATING RATES (DEG/DAY),
C
        DO 500 J      =  1,NVERT
           HEATS(J)   =  0.0
           HEATI(J)   =  0.0
           TERM1      =  FDEGDAY/(DPG(J)*G)
C
           IF(ISL .NE. 0) THEN
             DO 480 L     =  1,NSOLP
               HEATS(J)   =  HEATS(J)+(FNET(L,J+1)-FNET(L,J))*TERM1
 480         CONTINUE
           ENDIF
C
           IF(IR .NE. 0) THEN
             DO 490 L     =  NSOLP+1,NTOTAL
                HEATI(J)  =  HEATI(J)+( DIRECTU(L,J+1)-DIREC(L,J+1)
     1                       -(DIRECTU(L,J)-DIREC(L,J)) )*TERM1
 490         CONTINUE
           ENDIF
           HEAT(J)        =  HEATS(J)+HEATI(J)
c
c     Load heating rates [deg_K/s] into interface common block
c
           heats_aerad(j) =  heats(j)/scday
           heati_aerad(j) =  heati(j)/scday
C
 500    CONTINUE
c
c     Here we Calculate (4 * pi * mean_intensity) for the IR.
c
      IF (IR .NE. 0) THEN
        DO J = 1, NVERT
          DO L = NSOLP+1, NTOTAL
            TMI(L,J) = TMIU(L,J)+TMID(L,J)
          end do
        end do
      ENDIF
C
c     Here we compute the heating rates for droplets
c     (C11 converts W/m^2 to erg/cm^2)
c
      C11 = 1000.
      do ig = 1, NGROUP
       DO I = 1, NRAD
        DO J = 1, NVERT
C 
        QRAD(j) = 0.
C
        IF (IR .NE. 0) THEN
          DO L = NSOLP+1, NTOTAL
            X = TMI(L,J)-4.0*PI*PTEMP(L,J)
            if( abs(X/TMI(L,J)) .lt. EPSILON ) x = 0.
            QRAD(j) = QRAD(j) + X*C11*XSECTA(I,ig) *
     1     (RDQEXT(I,ig,NPROB(L))-QSCAT(I,ig,NPROB(L)))
          end do
        ENDIF
C        
        IF (ISL .NE. 0) THEN
          DO L = 1, NSOLP
            QRAD(j) = QRAD(j) + TMI(L,J)*C11*XSECTA(I,ig) *
     1     (RDQEXT(I,ig,NPROB(L))-QSCAT(I,ig,NPROB(L)))
          end do
        ENDIF
c
c     Load layer averages of droplet heating rates into interface common block
c
        if (j.eq.nvert) then
          qrad_aerad(i,j,ig) = qrad(j)
        else if (j.gt.1) then
          qrad_aerad(i,j-1,ig) = 0.5 * ( qrad(j) + qrad(j-1) )
        endif
C
        end do     ! j=1,nvert
       end do      ! i=1,nrad
      end do       ! ig=1,ngroup
c
c
c   Calculate some diagnostic quantities (formerly done in radout.f) and
c   load them into the interface common block.  None of these presently 
c   influence any microphysical processes -- hence, the following code only
c   needs to be executed before the aerosol model writes its output.  
c   Not all of the calculated quantities are presently being
c   loaded into the interface common block.  
c
c   Load optical depths into interface common block
c
        do i = 1, nwave
          opd_aerad(i) = uopd(i,nlayer)
        enddo
c
c     <tsLu> and <tsLd> are total upwelling and downwelling solar
c     fluxes at top-of-atmosphere
c
        tsLu = 0.
        tsLd = 0.
c
c     <fupbs>, <fdownbs>, and <fnetbs> are total upwelling, downwelling, and net
c     solar fluxes at grid boundaries
c
        do 507 j = 1, nlayer
          fupbs(j) = 0.
          fdownbs(j) = 0.
          fnetbs(j) = 0.
 507    continue
c
c     <fsLu> and <fsLd> are upwelling, downwelling, and net
c     solar fluxes at top-of-atmosphere (spectrally-resolved)
c
c     <alb_toa> is albedo at top-of-atmosphere (spectrally-resolved)
c
        do 509 i = 1, nsoL
          fsLu(i) = 0.
          fsLd(i) = 0.
          alb_toa(i) = 0.
 509    continue
c
c      <alb_tomi> and <alb_toai> are total solar albedos at top-of-model 
c      and top-of-atmosphere
c
        alb_tomi = 0.
        alb_toai = 0.
C
C     CALCULATE SOLAR ABSORBED BY GROUND, SOLNET, AND UPWARD AND DOWNWARD
C     LONGWAVE FLUXES AT SURFACE, XIRUP AND XIRDOWN (WATTS/M**2)
C
        SOLNET  = 0.0

        IF (ISL .NE. 0) THEN
          DO 510 L       =  1,NSOLP
            SOLNET = SOLNET - FNET(L,NLAYER)
            fp = ck1(L,1)*eL2(L,1) - ck2(L,1)*em2(L,1) + cp(L,1)
            fsLu( nprob(L) ) = fsLu( nprob(L) ) + fp
            do 510 j = 1, nlayer
              fp = ck1(L,j)*eL1(L,j) + ck2(L,j)*em1(L,j) + cpb(L,j)
              fupbs(j) = fupbs(j) + fp
              fnetbs(j) = fnetbs(j) + fnet(L,j)
              if (L.eq.nsolp) fdownbs(J) = fupbs(j) - fnetbs(j)
 510      CONTINUE

          do 508 i = 1, nsoL
            fsLd(i) = u0*solfx(i)
            alb_toa(i) = fsLu(i)/fsLd(i)
            tsLu = tsLu + fsLu(i)
            tsLd = tsLd + fsLd(i)
 508      continue
c
          alb_tomi = fupbs(1)/fdownbs(1)
          alb_toai = tsLu/tsLd
c
c      Load albedos into interface common block
c
          alb_toai_aerad = alb_toai
          alb_tomi_aerad = alb_tomi
          do i = 1, NSOL
            alb_toa_aerad(i) = alb_toa(i)
          enddo
c
c      Load fluxes into interface common block
c
          do j = 1, nlayer
            fsl_up_aerad(j) = fupbs(j)
            fsl_dn_aerad(j) = fdownbs(j)
          enddo

        ENDIF
c
c     <tiru> is total upwelling infrared flux at top-of-atmosphere;
c     <fupbi>, <fdownbi>, and <fnetbi> are total upwelling, downwelling, and net
c     infrared fluxes at grid boundaries
c
        tiru = 0.
        do 606 j = 1, nlayer
          fupbi(j)   =  0.0
          fdownbi(j)   =  0.0
          fnetbi(j)   =  0.0
 606    continue
c
c     <firu> is upwelling infrared flux at top-of-atmosphere (spectrally-resolved)
c
        do 609 i = 1, nir
          firu(i) = 0.
 609    continue

        XIRDOWN = 0.0
        XIRUP   = 0.0

        IF (IR .NE. 0) THEN
          DO 520 L        =  NSOLP+1,NTOTAL
             XIRDOWN = XIRDOWN + DIREC  (L,NLAYER)
             XIRUP   = XIRUP   + DIRECTU(L,NLAYER)
             firu( nprob(L)-nsol ) = firu( nprob(L)-nsol ) +
     1                               directu(L,1)
             do 520 j = 1, nlayer
               fupbi(j) = fupbi(j) + directu(L,j)
               fdownbi(j) = fdownbi(j) + direc  (L,j)
               fnetbi(j) = fnetbi(j) + directu(L,j) - direc(L,j)
 520      CONTINUE

          do 529 i = 1, nir
            tiru = tiru + firu(i)
 529      continue
c
c      Load fluxes into interface common block
c
          do j = 1, nlayer
            fir_up_aerad(j) = fupbi(j)
            fir_dn_aerad(j) = fdownbi(j)
          enddo


        ENDIF

      return
      END
