      SUBROUTINE ADD
C
C     ***************************************************************
C     *  Purpose             :  Defines source terms, form matrix   *
C     *                         for multiple layers and solve tri-  *
C     *                         diagnol equations to obtain mean    *
C     *                         intensity and net flux.             *
C     *  Subroutines Called  :  None                                *
C     *  Input               :  G0, U0, RSFX, TAUL, OPD             *
C     *  Output              :  B3, EE3, DIRECT, SFCS, CPB, CMB, DS *
C     * *************************************************************
C
      include 'globrad.h'
C
C     THIS SUBROUTINE FORMS THE MATRIX FOR THE MULTIPLE LAYERS AND
C     USES A TRIDIAGONAL ROUTINE TO FIND RADIATION IN THE ENTIRE
C     ATMOSPHERE.
C
C     ******************************
C     *   CALCULATIONS FOR SOLAR   *
C     ******************************
      IF(ISL .NE. 0)  THEN
        DU0                =  1./U0
        DO 10 J            =  1,NLAYER
            j1 = max( 1, j-1 )
            DO 10 L        =  1,NSOLP
               B3(L,J)     =  0.5*(1.-SQ3*G0(L,J)*U0)
               B4          =  1. - B3(L,J)
               X2          =  TAUL(L,J)*DU0
               EE3(L,J)    =  EXP(-X2)
               X3          =  OPD(L,J)*DU0
               EL3(L,J)    =  EXP(-X3)*SOL(L)
               DIRECT(L,J) =  U0*EL3(L,J)
               C1          =  B1(L,J) - DU0
               if( ABS(C1).lt.EPSILON ) c1 = SIGN(EPSILON,C1)
               C2          =  AK(L,J)*AK(L,J) - DU0*DU0
               if( ABS(C2) .le. EPSILON ) c2 = EPSILON
               CP1         =  W0(L,J)*(B3(L,J)*C1+B4*B2(L,J))/C2
               CPB(L,J)    =  CP1 * EL3(L,J)
               if( j .ne. 1 ) then
                 x4 = eL3(L,j1)
               else
                 x4 = soL(L)
               endif
               CP(L,J)     =  CP1 * X4
               CM1         =  ( CP1*B2(L,J) + W0(L,J)*B4 )/C1
               CMB(L,J)    =  CM1 * EL3(L,J)
               CM(L,J)     =  CM1 * X4
  10  CONTINUE
C
C       CALCULATE SFCS, THE SOURCE AT THE BOTTOM.
C
        DO 20 L            =  1,NSOLP
  20       SFCS(L)         =  DIRECT(L,NLAYER) * RSFX(L)
C
      END IF
C
C     ******************************
C     * CALCULATIONS FOR INFRARED. *
C     ******************************
C
      IF(IRS .NE. 0)  THEN
        DO 30 J           =   1,NLAYER
           KINDEX         = max(1,J-1)
           DO 30 L        = NSOLP+1,NTOTAL
              B3(L,J)     = 1.0/(B1(L,J)+B2(L,J))
              CP(L,J)     = (PTEMP(L,KINDEX)+SLOPE(L,J)*B3(L,J))*U1S(L)
              CPB(L,J)    = CP(L,J) + SLOPE(L,J)*TAUL(L,J)*U1S(L)
              CM(L,J)     = (PTEMP(L,KINDEX)-SLOPE(L,J)*B3(L,J))*U1S(L)
              CMB(L,J)    = CM(L,J) + SLOPE(L,J)*TAUL(L,J)*U1S(L)
              EL3(L,J)    = 0.0
              DIRECT(L,J) = 0.0
              EE3(L,J)    = 0.0
  30  CONTINUE
C
      DO 40 L             = NSOLP+1,NTOTAL
  40     SFCS(L)          = EMIS(L)*PTEMPG(L)*PI
C
      END IF
C
      J                =  0
      DO 42 JD         =  2,JN,2
         J             =  J + 1
         DO 42 L       =  LLS,LLA
C           HERE ARE THE EVEN MATRIX ELEMENTS
            DF(L,JD) = (CP(L,J+1) - CPB(L,J))*EM1(L,J+1) -
     1                 (CM(L,J+1) - CMB(L,J))*EM2(L,J+1)
C           HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
            DF(L,JD+1) =  EL2(L,J) * (CP(L,J+1)-CPB(L,J)) +
     2                    EL1(L,J) * (CMB(L,J) - CM(L,J+1))
  42  CONTINUE
C
C     HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
C     BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME NO
C     DIFFUSE RADIATION IS INCIDENT AT THE TOP.
C
      DO 44 L        = LLS,LLA
         DF(L,1)     = -CM(L,1)
         DF(L,JDBLE) = SFCS(L)+RSFX(L)*CMB(L,NLAYER)-CPB(L,NLAYER)
         DS(L,JDBLE) = DF(L,JDBLE)/BF(L,JDBLE)
  44     AS(L,JDBLE) = AF(L,JDBLE)/BF(L,JDBLE)
C
C     ********************************************
C     *     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
C     ********************************************
C
      DO 46 J               = 2, JDBLE
         DO 46 L            = LLS,LLA
            X               = 1./(BF(L,JDBLE+1-J) -
     1                        EF(L,JDBLE+1-J)*AS(L,JDBLE+2-J))
            AS(L,JDBLE+1-J) = AF(L,JDBLE+1-J)*X
            DS(L,JDBLE+1-J) = (DF(L,JDBLE+1-J) - EF(L,JDBLE+1-J)
     2                        *DS(L,JDBLE+2-J))*X
  46  CONTINUE
C
      DO 48 L       = LLS,LLA
  48     XK(L,1)    = DS(L,1)
C
      DO 50 J       = 2, JDBLE
         DO 50 L    = LLS,LLA
            XK(L,J) = DS(L,J) - AS(L,J)*XK(L,J-1)
  50  CONTINUE
C
C  ***************************************************************
C     CALCULATE LAYER COEFFICIENTS, NET FLUX AND MEAN INTENSITY
C  ***************************************************************
C
      do J = 1,NLAYER
        do L = LLS,LLA
          CK1(L,J)   = XK(L,2*J-1)                                         
          CK2(L,J)   = XK(L,2*J)                                          
          FNET(L,J)  = CK1(L,J)  *( EL1(L,J) -EL2(L,J))   +                
     3                 CK2(L,J) *( EM1(L,J)-EM2(L,J) ) + CPB(L,J) -  
     4                 CMB(L,J) - DIRECT(L,J)                        
C                                                                      
          TMI(L,J)   =  EL3(L,J) + U1I(L) *( CK1(L,J)  *                   
     5                ( EL1(L,J) + EL2(L,J))   +                         
     6                  CK2(L,J) *( EM1(L,J)+EM2(L,J) ) +
     7                  CPB(L,J) + CMB(L,J) )                             
        enddo
      enddo

      RETURN
      END
