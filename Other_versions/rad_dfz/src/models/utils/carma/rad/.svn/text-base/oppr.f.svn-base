      SUBROUTINE OPPR
C
C     **************************************************************
C     *  Purpose             :  CaLculates optical properties      *
C     *                         such as single scattering albedo,  *
C     *                         asymmetry parameter, etc.          *
C     *                         This routine is case dependent and *
C     *                         wiLL have to be repLaced by the    *
C     *                         user.                              *
C     *  Subroutines Called  :  None                               *
C     *  Input               :  PAH2O, RDH2O, CO2, O3, ETC         *
C     *  Output              :  TAUL, W0, G0, OPD, Y3              *
C     * ************************************************************
C
      include 'globrad.h'
C
C     W0(NWAVE,NLAYER) : SINGLE SCATTERING ALBEDO *** delta scaled ***
C     G0(NWAVE,NLAYER) : ASYMMETRY PARAMETER *** delta scaled ***
C     OPD(NWAVE,NLAYER): cumulative OPTICAL DEPTH *** delta scaled ***
C     SFL(NWAVE)       : SOLAR FLUX
C    uW0(NWAVE,NLAYER)  : unscaled SINGLE SCATTERING ALBEDO 
C    uG0(NWAVE,NLAYER)  : unscaled ASYMMETRY PARAMETER 
C    uTAUL(NWAVE,NLAYER): unscaled OPTICAL DEPTH of layer
C
C     ASSUME THAT P IS SAME ON ALL SIGMA LEVELS. IF PSURFACE
C     VARIES A LOT, THEN WE MAY NEED TO CALCULATE TAUO3,
C     TAUCO2, TAUO2 FOR EACH POINT.
C
C     NOTE : THE TOP LAYER IS A DUMMY. IT CONTAINS A DEFINED GAS
C            AMOUNT. DIFFERENT MODELS WILL REQUIRE DIFFERENT
C            TREATMENT OF THIS.
C     CALCULATE TOTAL OPTICAL DEPTH INCLUDING GASES. THEN
C     GIVEN THE AEROSOL OPTICAL DEPTHS AND CLOUD OPTICAL DEPTHS,
C     CALCULATE FINAL OPTICAL PROPERTIES. WE USE A DELTA
C     TWO STREAM APPROACH TO FIND W0, SINGLE SCATTERING ALBEDO,
C     G0, ASYMMMETRY PARAMETER, TAUL, LAYER OPTICAL DEPTH,
C     OPD, CUMULATIVE OPTICAL DEPTH TO BASE OF LAYER.
C
      DO 200  J          =   1,NLAYER
C
       DO 110 L           =   1,NWAVE
          TAUA(L,J)       =   0.0
          TAUS(L,J)       =   0.0
          G01(L,J)        =   0.0
 110   CONTINUE
C
       do ig            =   1,NGROUP
         DO 120 I       =   1,NRAD
             R2Z        =   XSECTA(I,ig)*CAER(I,J,ig)
          DO 115 L      =   1,NWAVE
             TAUA(L,J)  =   TAUA(L,J)+RDQEXT(I,ig,L)*R2Z
             TAUS(L,J)  =   TAUS(L,J)+QSCAT(I,ig,L)*R2Z
             G01(L,J)   =   G01(L,J) +QBRQS(I,ig,L)*R2Z
 115      CONTINUE
 120     CONTINUE
       enddo
C
       DO 130 L          = 1,NTOTAL
       TAUAER(L,J)       = MAX(TAUA(NPROB(L),J),EPSILON)
       WOL(L,J)          = TAUS(NPROB(L),J)/TAUAER(L,J)
       TAUAER(L,J)       = TAUA(NPROB(L),J)
       if( WOL(L,J) .ne. 0. ) then
         TTAS = TAUS(NPROB(L),J)
       else
         TTAS = 1.
       endif
       GOL(L,J)          = G01(NPROB(L),J)/TTAS
 130   CONTINUE
C
 200  CONTINUE
c
c     iradgas = 0: no gas in radiative xfer
c
      iradgas = 1     

      DO 500 J           = 1,NLAYER
          j1             = max( 1, j-1 )
c
c     Bergstrom water vapor continuum fix:
c
c      <qcorr> is layer average water vapor mixing ratio
c      <pcorr> is layer average pressure
c
c     For layer 0, calculate mixing ratio [g/g] from vapor column [g/cm^2]
c     and average pressure [dyne/cm^2]
c
          if( j .eq. 1 )then
            qcorr = rdh2o(1) * g / ptop
            pcorr = p(1) / 2.
          else
            qcorr = q(j1)
            pcorr = p(j1)
          endif
c
          cco = exp(1800./t(j1))*(qcorr*pcorr/2.87 + pcorr/4610.)
c
          DO 400 L       = LLS,LLA
             TAUH2O(L,J) = RDH2O(J)*PAH2O(L,J)
c
c     Bergstrom water vapor continuum fix (next two statements)
c
             if( L .GT. NSOLP+30 .AND. L .LE. NSOLP+36 ) then
               TAUH2O(L,J) = RDH2O(J)*PAH2O(L,J)*cco
             else
               TAUH2O(L,J) = TAUH2O(L,J)
             endif
             if (L.gt.nsolp+36) 
     1       tauh2o(L,j) = tauh2o(L,j) + cco*rdh2o(j)*contnm(L-nsolp)
 
             TAUL(L,J)   = TAUH2O(L,J)+TAUGAS(L,J)+
     1                        PARAY(L,J)+TAUAER(L,J)+TAUCLD(L,J)
             if (iradgas.eq.0) tauL(L,j) = tauaer(L,j)
             if( TAUL(L,J) .lt. EPSILON ) TAUL(L,J) = EPSILON
             utauL(L,j)  = TAUL(L,J)
             WOT         = (PARAY(L,J)+TAUAER(L,J)*WOL(L,J)+
     1                          TAUCLD(L,J)*WCLD(L,J))/TAUL(L,J)
             if (iradgas.eq.0) wot = woL(L,j)

             WOT         = min(1.-EPSILON,WOT)
             uw0(L,j)    = WOT
             DENOM       = (PARAY(L,J)+TAUCLD(L,J)*WCLD(L,J)+
     1                         TAUAER(L,J)*WOL(L,J))
             if( denom .le. EPSILON ) denom = EPSILON
             if( DENOM .GT. EPSILON ) then
               got = ( WCLD(L,J)*GCLD(L,J)*TAUCLD(L,J) + GOL(L,J)*
     $               WOL(L,J)*TAUAER(L,J) ) / DENOM
             else
               got = 0.
             endif
             if (iradgas.eq.0) got = goL(L,j)
             ug0(L,j)    = GOT
             FO          = GOT**2
             DEN         = 1.-WOT*FO
             TAUL(L,J)   = TAUL(L,J) * DEN
             W0(L,J)     = (1.-FO)*WOT/DEN
             G0(L,J)     = GOT/(1.+GOT)
             OPD(L,J)    = 0.0
             OPD(L,J)    = OPD(L,J1)+TAUL(L,J)
             uOPD(L,J)   = 0.0
             uOPD(L,J)   = uOPD(L,J1)+uTAUL(L,J)
 400        CONTINUE
C
          IF(IR .EQ. 1) THEN
             DO 450 I        =   1,NGAUSS
                DO 425 L     =   LLS,LLA
                   Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
 425            CONTINUE
 450         CONTINUE
          ENDIF
 500  CONTINUE

      RETURN
      END
