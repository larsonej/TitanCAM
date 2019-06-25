      SUBROUTINE DMIESS( RO, RFR, RFI, THETD, JX, QEXT, QSCAT, CTBRQS,
     1                   R, RE2, TMAG2, WVNO  )
C
C **********************************************************************
C    THIS SUBROUTINE COMPUTES MIE SCATTERING BY A STRATIFIED SPHERE,
C    I.E. A PARTICLE CONSISTING OF A SPHERICAL CORE SURROUNDED BY A
C    SPHERICAL SHELL.  THE BASIC CODE USED WAS THAT DESCRIBED IN THE
C    REPORT: " SUBROUTINES FOR COMPUTING THE PARAMETERS OF THE
C    ELECTROMAGNETIC RADIATION SCATTERED BY A SPHERE " J.V. DAVE,
C    I B M SCIENTIFIC CENTER, PALO ALTO , CALIFORNIA.
C    REPORT NO. 320 - 3236 .. MAY 1968 .
C
C    THE MODIFICATIONS FOR STRATIFIED SPHERES ARE DESCRIBED IN
C        TOON AND ACKERMAN, APPL. OPTICS, IN PRESS, 1981
C
C    THE PARAMETERS IN THE CALLING STATEMENT ARE DEFINED AS FOLLOWS :
C      RO IS THE OUTER (SHELL) RADIUS;
C      R  IS THE CORE RADIUS;
C      RFR, RFI  ARE THE REAL AND IMAGINARY PARTS OF THE SHELL INDEX
C          OF REFRACTION IN THE FORM (RFR - I* RFI);
C      RE2, TMAG2  ARE THE INDEX PARTS FOR THE CORE;
C          ( WE ASSUME SPACE HAS UNIT INDEX. )
C      THETD(J): ANGLE IN DEGREES BETWEEN THE DIRECTIONS OF THE INCIDENT
C          AND THE SCATTERED RADIATION.  THETD(J) IS< OR= 90.0
C          IF THETD(J) SHOULD HAPPEN TO BE GREATER THAN 90.0, ENTER WITH
C          SUPPLEMENTARY VALUE, SEE COMMENTS BELOW ON ELTRMX;
C      JX: TOTAL NUMBER OF THETD FOR WHICH THE COMPUTATIONS ARE
C          REQUIRED.  JX SHOULD NOT EXCEED IT UNLESS THE DIMENSIONS
C          STATEMENTS ARE APPROPRIATEDLY MODIFIED;
C
C      THE DEFINITIONS FOR THE FOLLOWING SYMBOLS CAN BE FOUND IN"LIGHT
C          SCATTERING BY SMALL PARTICLES,H.C.VAN DE HULST, JOHN WILEY '
C          SONS, INC., NEW YORK, 1957" .
C      QEXT: EFFICIENCY FACTOR FOR EXTINCTION,VAN DE HULST,P.14 ' 127.
C      QSCAT: EFFICIENCY FACTOR FOR SCATTERING,V.D. HULST,P.14 ' 127.
C      CTBRQS: AVERAGE(COSINE THETA) * QSCAT,VAN DE HULST,P.128
C      ELTRMX(I,J,K): ELEMENTS OF THE TRANSFORMATION MATRIX F,V.D.HULST
C          ,P.34,45 ' 125. I=1: ELEMENT M SUB 2..I=2: ELEMENT M SUB 1..
C          I = 3: ELEMENT S SUB 21.. I = 4: ELEMENT D SUB 21..
C      ELTRMX(I,J,1) REPRESENTS THE ITH ELEMENT OF THE MATRIX FOR
C          THE ANGLE THETD(J).. ELTRMX(I,J,2) REPRESENTS THE ITH ELEMENT
C          OF THE MATRIX FOR THE ANGLE 180.0 - THETD(J) ..
C      QBS IS THE BACK SCATTER CROSS SECTION.
C
C      IT: IS THE DIMENSION OF THETD, ELTRMX, CSTHT, PI, TAU, SI2THT,
C          IT MUST CORRESPOND EXACTLY TO THE SECOND DIMENSION OF ELTRMX.  
C      IACAP IS THE DIMENSION OF ACAP 
C          IN THE ORIGINAL PROGRAM THE DIMENSION OF ACAP WAS 7000.
C          FOR CONSERVING SPACE THIS SHOULD BE NOT MUCH HIGHER THAN
C          THE VALUE, N=1.1*(NREAL**2 + NIMAG**2)**.5 * X + 1
C      WVNO: 2*PI / WAVELENGTH
C
C    ALSO THE SUBROUTINE COMPUTES THE CAPITAL A FUNCTION BY MAKING USE O
C    DOWNWARD RECURRENCE RELATIONSHIP.
C
C      TA(1): REAL PART OF WFN(1).  TA(2): IMAGINARY PART OF WFN(1).
C      TA(3): REAL PART OF WFN(2).  TA(4): IMAGINARY PART OF WFN(2).
C      TB(1): REAL PART OF FNA.     TB(2): IMAGINARY PART OF FNA.
C      TC(1): REAL PART OF FNB.     TC(2): IMAGINARY PART OF FNB.
C      TD(1): REAL PART OF FNAP.    TD(2): IMAGINARY PART OF FNAP.
C      TE(1): REAL PART OF FNBP.    TE(2): IMAGINARY PART OF FNBP.
C      FNAP, FNBP  ARE THE PRECEDING VALUES OF FNA, FNB RESPECTIVELY.
C **********************************************************************
c
c   Include implicit declarations
c
c   New idea: implicit double precision, explicit single precision for
c   arguments
c
c     include 'precision.h'
c
c  Define implicit type for floating point variables
c
      implicit double precision ( a-h, o-z )
c
c
c  Define implicit type for integer variables
c
      implicit integer ( i-n )
c
c
      parameter( EPSILON_MIE = 1.d-14 )
c
c
c  Explicitly declare floating point arguments
c
      double precision RO, RFR, RFI, THETD, QEXT, QSCAT, CTBRQS,
     1                   R, RE2, TMAG2, WVNO
c
c   Define various dimensions of local arrays
c
      parameter( IACAP = 200000, IT = 1 )
c
c   Declare local variables
c
c     COMPLEX FNAP,   FNBP,   ACAP(IACAP), W,
      double COMPLEX FNAP,   FNBP,   ACAP(IACAP), W,
     2        FNA,    FNB,    RF,       RRF,
     3        RRFX,   WM1,    FN1,      FN2,
     4        TC1,    TC2,    WFN(2),   Z(4),
     5        K1,     K2,     K3,
     6        RC,     U(8),   DH1,
     7        DH2,    DH4,    P24H24,   P24H21,
     8        PSTORE, HSTORE, DUMMY,    DUMSQ
C
      COMMON / WARRAY / W(3,IACAP)
C
      DIMENSION     T(5),      TA(4),     TB(2),        TC(2),
     2              TD(2),     TE(2),     PI( 3,IT ),   TAU( 3,IT ),
     3              CSTHT(IT), THETD(IT), SI2THT(IT),   ELTRMX( 4,IT,2 )
C
      EQUIVALENCE   (FNA,TB(1)),(FNB,TC(1)),(FNAP,TD(1)),(FNBP,TE(1))
c
c
c  Some compilers (e.g. absoft) don't the support imag(z) generic intrinsic function.
c  In that situation, uncomment the following 2 lines to define stmt function to
c  redefine imag() function to the alternative function.  I.e. change "alt_function"
c  in following stmt function to the appropriate function that will return the
c  imaginary part of a double complex for the compiler being used.
c  For compilers that support imag(z), simply leave following 2 statments commented.
c  The "c--alt_imag" comment prefix is designed to be recognized by the
c  automated editing script create_dmiess, so please do not modify this prefix
c  in the dmiess.f.template file.
c  -bm  Nov-1999
c
c--alt_imag      double complex z_dum_arg
c--alt_imag      imag(z_dum_arg) = alt_function(z_dum_arg)
C
C
C   IF THE CORE IS SMALL SCATTERING IS COMPUTED FOR THE SHELL ONLY
C
      IFLAG = 1
      IF ( R/RO .LT. 1.d-6 )   IFLAG = 2
      IF ( JX .LE. IT )   GO TO 20
         WRITE( *,7 )
         WRITE( *,6 )
         STOP 30
   20 RF =  CMPLX( RFR,  -RFI )
      RC =  CMPLX( RE2, -TMAG2 )
      X  =  RO * WVNO
      K1 =  RC * WVNO
      K2 =  RF * WVNO
      K3 =  CMPLX( WVNO, 0.0 )
      Z(1) =  K2 * RO
      Z(2) =  K3 * RO
      Z(3) =  K1 * R
      Z(4) =  K2 * R
      X1   =  REAL( Z(1) )
      X4   =  REAL( Z(4) )
c aimag for AbSoft
c      Y1   =  aimag( Z(1) )
c      Y4   =  aimag( Z(4) )
      Y1   =  imag( Z(1) )
      Y4   =  imag( Z(4) )
      RRF  =  1.0 / RF
      RX   =  1.0 / X
      RRFX =  RRF * RX
      T(1) =  ( X**2 ) * ( RFR**2 + RFI**2 )
      T(1) =  SQRT( T(1) )
      NMX1 =  1.10 * T(1)
C
      IF ( NMX1 .LE. IACAP-1 )   GO TO 21
         WRITE(*,8)
         STOP 32
   21 NMX2 = T(1)
      IF ( NMX1 .GT.  150 )   GO TO 22
         NMX1 = 150
         NMX2 = 135
C
   22 ACAP( NMX1+1 )  =  ( 0.0,0.0 )
      IF ( IFLAG .EQ. 2 )   GO TO 26
         DO N = 1,3
           W( N,NMX1+1 )  =  ( 0.0,0.0 )
         enddo
   26 CONTINUE
      DO N = 1,NMX1
         NN = NMX1 - N + 1
         ACAP(NN) = (NN+1) * RRFX - 1.0 / ( (NN+1) * RRFX + ACAP(NN+1) )
         IF ( IFLAG .EQ. 2 )   GO TO 23
            DO M = 1,3
              W( M,NN ) = (NN+1) / Z(M+1)  -
     1                   1.0 / (  (NN+1) / Z(M+1)  +  W( M,NN+1 )  )
            enddo
   23 CONTINUE
      enddo
C
      DO 30   J = 1,JX
      IF ( THETD(J) .LT. 0.0 )  THETD(J) =  ABS( THETD(J) )
      IF ( THETD(J) .GT. 0.0 )  GO TO 24
      CSTHT(J)  = 1.0
      SI2THT(J) = 0.0
      GO TO 30
   24 IF ( THETD(J) .GE. 90.0 )  GO TO 25
      T(1)      =  ( 3.14159265359 * THETD(J) ) / 180.0
      CSTHT(J)  =  COS( T(1) )
      SI2THT(J) =  1.0 - CSTHT(J)**2
      GO TO 30
   25 IF ( THETD(J) .GT. 90.0 )  GO TO 28
      CSTHT(J)  =  0.0
      SI2THT(J) =  1.0
      GO TO 30
   28 WRITE( *,5 )  THETD(J)
      WRITE( *,6 )
      STOP 34
   30 CONTINUE
C
      DO 35  J = 1,JX
      PI(1,J)  =  0.0
      PI(2,J)  =  1.0
      TAU(1,J) =  0.0
      TAU(2,J) =  CSTHT(J)
   35 CONTINUE
C
C INITIALIZATION OF HOMOGENEOUS SPHERE
C
      T(1)   =  COS(X)
      T(2)   =  SIN(X)
      WM1    =  CMPLX( T(1),-T(2) )
      WFN(1) =  CMPLX( T(2), T(1) )
      TA(1)  =  T(2)
      TA(2)  =  T(1)
      WFN(2) =  RX * WFN(1) - WM1
      TA(3)  =  REAL(WFN(2))
c aimag for AbSoft
c     TA(4)  =  aimag(WFN(2))
      TA(4)  =  imag(WFN(2))
C
      IF ( IFLAG .EQ. 2 )   GO TO 560
      N = 1
C
C INITIALIZATION PROCEDURE FOR STRATIFIED SPHERE BEGINS HERE
C
      SINX1   =  SIN( X1 )
      SINX4   =  SIN( X4 )
      COSX1   =  COS( X1 )
      COSX4   =  COS( X4 )
      EY1     =  EXP( Y1 )
      E2Y1    =  EY1 * EY1
      EY4     =  EXP( Y4 )
      EY1MY4  =  EXP( Y1 - Y4 )
      EY1PY4  =  EY1 * EY4
      EY1MY4  =  EXP( Y1 - Y4 )
      AA  =  SINX4 * ( EY1PY4 + EY1MY4 )
      BB  =  COSX4 * ( EY1PY4 - EY1MY4 )
      CC  =  SINX1 * ( E2Y1 + 1.0 )
      DD  =  COSX1 * ( E2Y1 - 1.0 )
      DENOM   =  1.0  +  E2Y1 * ( 4.0 * SINX1 * SINX1 - 2.0 + E2Y1 )
      REALP   =  ( AA * CC  +  BB * DD ) / DENOM
      AMAGP   =  ( BB * CC  -  AA * DD ) / DENOM
      DUMMY   =  CMPLX( REALP, AMAGP )
      AA  =  SINX4 * SINX4 - 0.5
      BB  =  COSX4 * SINX4
      P24H24  =  0.5 + CMPLX( AA,BB ) * EY4 * EY4
      AA  =  SINX1 * SINX4  -  COSX1 * COSX4
      BB  =  SINX1 * COSX4  +  COSX1 * SINX4
      CC  =  SINX1 * SINX4  +  COSX1 * COSX4
      DD  = -SINX1 * COSX4  +  COSX1 * SINX4
      P24H21  =  0.5 * CMPLX( AA,BB ) * EY1 * EY4  +
     2           0.5 * CMPLX( CC,DD ) * EY1MY4
      DH4  =  Z(4) / ( 1.0 + ( 0.0,1.0 ) * Z(4) )  -  1.0 / Z(4)
      DH1  =  Z(1) / ( 1.0 + ( 0.0,1.0 ) * Z(1) )  -  1.0 / Z(1)
      DH2  =  Z(2) / ( 1.0 + ( 0.0,1.0 ) * Z(2) )  -  1.0 / Z(2)
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
C
C NOTE:  THE DEFINITIONS OF U(I) IN THIS PROGRAM ARE NOT THE SAME AS
C        THE USUBI DEFINED IN THE ARTICLE BY TOON AND ACKERMAN.  THE
C        CORRESPONDING TERMS ARE:
C          USUB1 = U(1)                       USUB2 = U(5)
C          USUB3 = U(7)                       USUB4 = DUMSQ
C          USUB5 = U(2)                       USUB6 = U(3)
C          USUB7 = U(6)                       USUB8 = U(4)
C          RATIO OF SPHERICAL BESSEL FTN TO SPHERICAL HENKAL FTN = U(8)
C
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0,-1.0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
C
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /
     2               ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /
     2               ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
      GO TO 561
  560 TC1  =  ACAP(1) * RRF  +  RX
      TC2  =  ACAP(1) * RF   +  RX
      FNA  =  ( TC1 * TA(3)  -  TA(1) ) / ( TC1 * WFN(2)  -  WFN(1) )
      FNB  =  ( TC2 * TA(3)  -  TA(1) ) / ( TC2 * WFN(2)  -  WFN(1) )
C
  561 CONTINUE
      FNAP = FNA
      FNBP = FNB
      T(1) = 1.50
C
C    FROM HERE TO THE STATMENT NUMBER 90, ELTRMX(I,J,K) HAS
C    FOLLOWING MEANING:
C    ELTRMX(1,J,K): REAL PART OF THE FIRST COMPLEX AMPLITUDE.
C    ELTRMX(2,J,K): IMAGINARY PART OF THE FIRST COMPLEX AMPLITUDE.
C    ELTRMX(3,J,K): REAL PART OF THE SECOND COMPLEX AMPLITUDE.
C    ELTRMX(4,J,K): IMAGINARY PART OF THE SECOND COMPLEX AMPLITUDE.
C    K = 1 : FOR THETD(J) AND K = 2 : FOR 180.0 - THETD(J)
C    DEFINITION OF THE COMPLEX AMPLITUDE: VAN DE HULST,P.125.
C
      TB(1) = T(1) * TB(1)
      TB(2) = T(1) * TB(2)
      TC(1) = T(1) * TC(1)
      TC(2) = T(1) * TC(2)
      DO 60 J = 1,JX
          ELTRMX(1,J,1) = TB(1) * PI(2,J) + TC(1) * TAU(2,J)
          ELTRMX(2,J,1) = TB(2) * PI(2,J) + TC(2) * TAU(2,J)
          ELTRMX(3,J,1) = TC(1) * PI(2,J) + TB(1) * TAU(2,J)
          ELTRMX(4,J,1) = TC(2) * PI(2,J) + TB(2) * TAU(2,J)
          ELTRMX(1,J,2) = TB(1) * PI(2,J) - TC(1) * TAU(2,J)
          ELTRMX(2,J,2) = TB(2) * PI(2,J) - TC(2) * TAU(2,J)
          ELTRMX(3,J,2) = TC(1) * PI(2,J) - TB(1) * TAU(2,J)
          ELTRMX(4,J,2) = TC(2) * PI(2,J) - TB(2) * TAU(2,J)
   60 CONTINUE
C
      QEXT   = 2.0 * ( TB(1) + TC(1))
      QSCAT  = ( TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2 ) / 0.75
      CTBRQS = 0.0
      QBSR   = -2.0*(TC(1) - TB(1))
      QBSI   = -2.0*(TC(2) - TB(2))
      RMM    = -1.0
      N = 2
   65 T(1) = 2*N - 1
      T(2) =   N - 1
      T(3) = 2*N + 1
      DO 70  J = 1,JX
          PI(3,J)  = ( T(1) * PI(2,J) * CSTHT(J) - N * PI(1,J) ) / T(2)
          TAU(3,J) = CSTHT(J) * ( PI(3,J) - PI(1,J) )  -
     1                          T(1) * SI2THT(J) * PI(2,J)  +  TAU(1,J)
   70 CONTINUE
C
C HERE SET UP HOMOGENEOUS SPHERE
C
      WM1    =  WFN(1)
      WFN(1) =  WFN(2)
      TA(1)  =  REAL(WFN(1))
c aimag for AbSoft
c      TA(2)  =  aimag(WFN(1))
c      TA(4)  =  aimag(WFN(2))
      TA(2)  =  imag(WFN(1))
      TA(4)  =  imag(WFN(2))
      WFN(2) =  T(1) * RX * WFN(1)  -  WM1
      TA(3)  =  REAL(WFN(2))
C
      IF ( IFLAG .EQ. 2 )   GO TO 1000
C
C HERE SET UP STRATIFIED SPHERE
C
      DH2  =  - N / Z(2)  +  1.0 / ( N / Z(2) - DH2 )
      DH4  =  - N / Z(4)  +  1.0 / ( N / Z(4) - DH4 )
      DH1  =  - N / Z(1)  +  1.0 / ( N / Z(1) - DH1 )
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
C
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0,-1.0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
C
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /
     2               ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /
     2               ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
C
 1000 CONTINUE
      TC1  =  ACAP(N) * RRF  +  N * RX
      TC2  =  ACAP(N) * RF   +  N * RX
      FN1  =  ( TC1 * TA(3)  -  TA(1) ) /  ( TC1 * WFN(2) - WFN(1) )
      FN2  =  ( TC2 * TA(3)  -  TA(1) ) /  ( TC2 * WFN(2) - WFN(1) )
      M    =  WVNO * R
      IF ( N .LT. M )   GO TO 1002
      IF ( IFLAG .EQ. 2 )   GO TO 1001
      IF ( abs(  ( FN1-FNA ) / FN1  )  .LT.  EPSILON_MIE   .AND.
     1     abs(  ( FN2-FNB ) / FN2  )  .LT . EPSILON_MIE  ) IFLAG = 2
      IF ( IFLAG .EQ. 1 )   GO TO 1002
 1001 FNA  =  FN1
      FNB  =  FN2
C
 1002 CONTINUE
      T(5)  =  N
      T(4)  =  T(1) / ( T(5) * T(2) )
      T(2)  =  (  T(2) * ( T(5) + 1.0 )  ) / T(5)
C
      CTBRQS  =  CTBRQS  +  T(2) * ( TD(1) * TB(1)  +  TD(2) * TB(2)
     1                   +           TE(1) * TC(1)  +  TE(2) * TC(2) )
     2                   +  T(4) * ( TD(1) * TE(1)  +  TD(2) * TE(2) )
      QEXT    =   QEXT  +  T(3) * ( TB(1) + TC(1) )
c     $        T(3), TB(1), TC(1), QEXT
      T(4)    =  TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2
      QSCAT   =  QSCAT  +  T(3) * T(4)
      RMM     =  -RMM
      QBSR    =  QBSR + T(3)*RMM*(TC(1) - TB(1))
      QBSI    =  QBSI + T(3)*RMM*(TC(2) - TB(2))
C
      T(2)    =  N * (N+1)
      T(1)    =  T(3) / T(2)
      K = (N/2)*2
      DO 80 J = 1,JX
       ELTRMX(1,J,1) = ELTRMX(1,J,1)+T(1)*(TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,1) = ELTRMX(2,J,1)+T(1)*(TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,1) = ELTRMX(3,J,1)+T(1)*(TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,1) = ELTRMX(4,J,1)+T(1)*(TC(2)*PI(3,J)+TB(2)*TAU(3,J))
      IF ( K .EQ. N )  THEN
       ELTRMX(1,J,2) =ELTRMX(1,J,2)+T(1)*(-TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,2) =ELTRMX(2,J,2)+T(1)*(-TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,2) =ELTRMX(3,J,2)+T(1)*(-TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,2) =ELTRMX(4,J,2)+T(1)*(-TC(2)*PI(3,J)+TB(2)*TAU(3,J))
      ELSE
       ELTRMX(1,J,2) = ELTRMX(1,J,2)+T(1)*(TB(1)*PI(3,J)-TC(1)*TAU(3,J))
       ELTRMX(2,J,2) = ELTRMX(2,J,2)+T(1)*(TB(2)*PI(3,J)-TC(2)*TAU(3,J))
       ELTRMX(3,J,2) = ELTRMX(3,J,2)+T(1)*(TC(1)*PI(3,J)-TB(1)*TAU(3,J))
       ELTRMX(4,J,2) = ELTRMX(4,J,2)+T(1)*(TC(2)*PI(3,J)-TB(2)*TAU(3,J))
      END IF
   80 CONTINUE
C
      IF ( T(4) .LT. EPSILON_MIE )   GO TO 100
      N = N + 1
      DO 90 J = 1,JX
         PI(1,J)   =   PI(2,J)
         PI(2,J)   =   PI(3,J)
         TAU(1,J)  =  TAU(2,J)
         TAU(2,J)  =  TAU(3,J)
   90 CONTINUE
      FNAP  =  FNA
      FNBP  =  FNB
      IF ( N .LE. NMX2 )   GO TO 65
         WRITE( *,8 )
         STOP 36
  100 continue
      DO J = 1,JX
      DO K = 1,2
         DO  115  I= 1,4
         T(I)  =  ELTRMX(I,J,K)
  115    CONTINUE
         ELTRMX(2,J,K)  =      T(1)**2  +  T(2)**2
         ELTRMX(1,J,K)  =      T(3)**2  +  T(4)**2
         ELTRMX(3,J,K)  =  T(1) * T(3)  +  T(2) * T(4)
         ELTRMX(4,J,K)  =  T(2) * T(3)  -  T(4) * T(1)
      enddo
      enddo
      T(1)    =    2.0 * RX**2
      QEXT    =   QEXT * T(1)
      QSCAT   =  QSCAT * T(1)
      CTBRQS  =  2.0 * CTBRQS * T(1)
C
C QBS IS THE BACK SCATTER CROSS SECTION
C
c      PIG   = ACOS(-1.0)
c      RXP4  = RX*RX/(4.0*PIG)
c      QBS   = RXP4*(QBSR**2 + QBSI**2)
C
    5 FORMAT( 10X,' THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN
     1 90.0 DEGREES. IT IS ', E15.4 )
    6 FORMAT( // 10X, 'PLEASE READ COMMENTS.' // )
    7 FORMAT( // 10X, 'THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT')
    8 FORMAT( // 10X, 'THE UPPER LIMIT FOR ACAP IS NOT ENOUGH. SUGGEST
     1 GET DETAILED OUTPUT AND MODIFY SUBROUTINE' // )
C
      RETURN
      END