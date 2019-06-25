      SUBROUTINE TWOSTR
C
C    ******************************************************************
C    *  Purpose             :  Defines matrix properties and sets up  *
C    *                         matrix coefficients that do not depend *
C    *                         on zenith angle or temperature.        *
C    *  Subroutines Called  :  None                                   *
C    *  Input               :  W0, G0                                 *
C    *  Output              :  B1, B2, GAMI, ACON, EL1, AF, ETC       *
C    * ****************************************************************
C
      include 'globrad.h'
      parameter (TWO = 2.d0)
C
       DO 10 L    =  LLS,LLA
          if( L .LE. NSOLP ) then
            U1I(L) = SQ3
          else
            U1I(L) = TWO
          endif
  10      U1S(L)  =  TPI/U1I(L)
      
C
C      HERE WE DEFINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
C      OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
C      NEEDED FOR MATRIX.
C
       DO 14 J          =  1,NLAYER
          DO 14 L       =  LLS,LLA
C            THESE ARE FOR TWO STREAM AND HEMISPHERIC MEANS
             B1(L,J)    =  0.5*U1I(L)*(2. - W0(L,J)*(1. + G0(L,J)))
             B2(L,J)    =  0.5*U1I(L)*W0(L,J)*(1. - G0(L,J))
             AK(L,J)    =  SQRT(ABS(B1(L,J)**2 - B2(L,J)**2))
             GAMI(L,J)  =  B2(L,J)/(B1(L,J) + AK(L,J))
             EE1(L,J)   =  EXP(-AK(L,J)*TAUL(L,J))
             EL1(L,J)   =  1.0 + GAMI(L,J) *EE1(L,J)                           
             EM1(L,J)   =  1.0 - GAMI(L,J) * EE1(L,J)                       
             EL2(L,J)   =  GAMI(L,J) + EE1(L,J)                               
             EM2(L,J)   =  GAMI(L,J) - EE1(L,J)                               
  14  CONTINUE
C
C     WE SEEK TO SOLVE AX(L-1)+BX(L)+EX(L+1) = D.
C     L=2N FOR EVEN L, L=N+1 FOR ODD L. THE MEAN INTENSITY (TMI/4PI)
C     AND THE NET FLUX (FNET) ARE RELATED TO X'S AS NOTED IN ADD.
C     FIRST WE SET UP THE COEFFICIENTS THAT ARE INDEPENDENT OF SOLAR
C     ANGLE OR TEMPARATURE: A(I),B(I),E(I). D(I) IS DEFINED IN ADD.
C
      J                 =  0
      DO 18 JD          =  2,JN,2
         J              =  J + 1
         DO 18 L        =  LLS,LLA
C          HERE ARE THE EVEN MATRIX ELEMENTS
             AF(L,JD)   =  EM1(L,J+1)*EL1(L,J)-EM2(L,J+1)*EL2(L,J)
             BF(L,JD)   =  EM1(L,J+1)* EM1(L,J)-EM2(L,J+1)*EM2(L,J)
             EF(L,JD)   =  EL1(L,J+1)*EM2(L,J+1) - EL2(L,J+1)*EM1(L,J+1)
C          HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP. 
             AF(L,JD+1) =  EM1(L,J)*EL2(L,J)-EL1(L,J)*EM2(L,J)
             BF(L,JD+1) =  EL1(L,J+1)*EL1(L,J) - EL2(L,J+1)*EL2(L,J)      
             EF(L,JD+1) =  EL2(L,J)*EM2(L,J+1)-EL1(L,J)*EM1(L,J+1)
  18  CONTINUE
C
C     HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
C     BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
C     NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
C
      DO 20 L        = LLS,LLA
         AF(L,1)     = 0.0
         BF(L,1)     = EL1(L,1)
         EF(L,1)     = -EM1(L,1)
         AF(L,JDBLE) = EL1(L,NLAYER)-RSFX(L)*EL2(L,NLAYER)
         BF(L,JDBLE) = EM1(L,NLAYER)-RSFX(L)*EM2(L,NLAYER)
         EF(L,JDBLE) = 0.0
  20  CONTINUE
      RETURN
      END
