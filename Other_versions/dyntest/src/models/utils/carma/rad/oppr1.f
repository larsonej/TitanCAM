      SUBROUTINE OPPR1
C
C     **********************************************************
C     *  Purpose             :  Calculate Planck Function and  *
C     *                         and its derivative at ground   *
C     *                         and at all altitudes.          *
C     *  Subroutines Called  :  None                           *
C     *  Input               :  TGRND, NLOW, WEIGHT            *
C     *  Output              :  PTEMP, PTEMPG, SLOPE           *
C     * ********************************************************
C
      include 'globrad.h'
C
      DIMENSION  PTEMP2(NTOTAL), PLTEMP1(NTOTAL)
C
C     **************************************
C     * CALCULATE PTEMP AND SLOPE          *
C     **************************************
C
C     CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE GROUND.
      ITG                 = ANINT(100.*TGRND) - NLOW
C
       CALL GATHER(NIRP,PLTEMP1,PLANK(1,ITG),LTEMP)
       DO 100 L            =   NSOLP+1,NTOTAL
          PTEMPG(L)        =   PLTEMP1(L-NSOLP)*WEIGHT(L)
 100   CONTINUE
C
      if( iblackbody_above .ne. 0 )then
C
C       CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE TOP 
C       OF THE MODEL.
        ITP                 = ANINT(100.*t_above) - NLOW
C
         CALL GATHER(NIRP,PLTEMP1,PLANK(1,ITP),LTEMP)
         DO 400 L            =   NSOLP+1,NTOTAL
            PTEMPT(L)        =   PLTEMP1(L-NSOLP)*WEIGHT(L)
 400     CONTINUE
C
      endif
C
        DO 300 J            =   1,NLAYER
            kindex          = max( 1, j-1 )
            IT1             = ANINT(100.*TT(J)) - NLOW
            CALL GATHER(NIRP,PTEMP2,PLANK(1,IT1),LTEMP)
C
C           KINDEX MAKES THE TOP LAYER ISOTHERMAL. USING KINDEX, FIND
C           PLANK FUNCTION AT BOTTOM OF EACH LAYER.
C           NOTE: IF YOU FORCE SLOPE=0, THEN YOU HAVE ISOTHERMAL
C           LAYERS WITH TT(J) CORRESPONDING TO AVERAGE TEMPERATURE
C           OF LAYER AND TT(NLAYER) SHOULD BE SET TO TGRND.
            DO 200 L        = NSOLP+1,NTOTAL
               PTEMP(L,J)   = PTEMP2(L-NSOLP)*WEIGHT(L)
               SLOPE(L,J)   = (PTEMP(L,J)-PTEMP(L,KINDEX))/TAUL(L,J)
               if( TAUL(L,J) .le. 1.0E-6 ) SLOPE(L,J) = 0.
 200        CONTINUE
 300     CONTINUE
C
      RETURN
      END
