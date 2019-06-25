      SUBROUTINE PLNK(E,T1,D)
C
C     ******************************************************
C     *  Purpose             :  Calculate Planck Function  *
C     *  Subroutines Called  :  None                       *
C     *  Input               :  WAVE, NCOUNT               *
C     *  Output              :  PLANK                      *
C     * ****************************************************
C
C  THIS SUBROUTINE COMPUTES THE INTEGRAL OF THE PLANCK FUNCTION BETWEEN
C  ZERO AND THE SPECIFIED VALUE OF LAMBDA.  THUS (USING XL AS LAMBDA)
C  WE WANT TO INTEGRATE
C  R = INTEGRAL(XL=0 TO XL=XLSPEC) ( C1*XL**-5* / (EXP(C2/XL*T)-1) )*DXL
C  SUBSTITUTING U=C2/(XL*T), THE INTEGRAL BECOMES
C  R = A CONSTANT TIMES INTEGRAL (USPEC TO INFINITY) OF
C            ( U**3 / (EXP(U) - 1) )*DU
C  THE APPROXIMATIONS SHOWN HERE ARE ON PAGE 998 OF ABRAMOWITZ AND SEGUN
C  UNDER THE HEADING OF DEBYE FUNCTIONS.  C2 IS THE PRODUCT OF PLANCK'S
C  CONSTANT AND THE SPEED OF LIGHT DIVIDED BY BOLTZMANN'S CONSTANT.
C  C2 = 14390 WHEN LAMBDA IS IN MICRONS.
C  THE FACTOR 0.15399 IS THE RECIPROCAL OF SIX TIMES
C  THE SUM OF (1/N**2) FOR ALL N FROM ONE TO INFINITY.  IT IS CHOSEN TO
C  NORMALIZE THE INTEGRAL TO A MAXIMUM VALUE OF UNITY.
C  RADIATION IN REAL UNITS IS OBTAINED BY MULTIPLYING THE INTEGRAL BY
C  THE STEFAN-BOLTZMANN CONSTANT TIMES T**4.
c
c
c   Include implicit declarations
c
      include 'precision.h'
c
c   Local declarations
C
      DIMENSION AM(5)
C
      D            =   0.0
      V1           =   E/T1
C
      IF (V1 .LE. 1.) THEN
         D         =  1.0 - 0.15399*V1**3 *
     1                (1./3.-V1/8. + V1**2/60. - V1**4/5040. +
     2                V1**6/272160. - V1**8/13305600         )
      ENDIF
C
      IF ( V1 .GT. 1. .AND. V1 .LE. 50.) THEN
         DO 100 M   =  1,5
            A       =  FLOAT(M)*V1
            AM(M)   =  0.15399 * EXP(-A)/M**4 *
     1                 (((A+3.)*A+6.)*A+6.)
 100     CONTINUE
C
         D          =  AM(1)+AM(2)+AM(3)+AM(4)+AM(5)
      ENDIF
C
      D             =  D*T1**4
C
      RETURN
      END
