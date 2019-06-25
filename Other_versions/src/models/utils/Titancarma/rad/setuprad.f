      SUBROUTINE SETUPRAD
C
C     *********************************************************
C     *  Purpose            :  Defines all constants, and     *
C     *                        calculates pressure averaged   *
C     *                        absorption coefficients.       *
C     * *******************************************************
C
      include 'globrad.h'
C
C **********************************************************************
C
C           LOCAL DECLARATIONS
C
C **********************************************************************
C
      DIMENSION AKO3(4,6), AKCO2(6,6), AKH2O(54,6), PJ(6)
c
      logical all_ok
C
C **********************************************************************
C
C            DEFINE CONSTANTS 
C
C **********************************************************************
C
C     NVERT   - NUMBER OF VERTICAL LAYERS
C     NLAYER  - NUMBER OF EDGES OF VERTICAL LAYERS
C
      DATA NVERT   / IVERT        /
      DATA NLAYER  / ILAYER       /
C
C     NRAD    - NUMBER OF RADIUS BINS
C
      DATA NRAD    / IRAD         /
C
C     REAL REFRACTIVE INDEX FOR LIQUID WATER
C
      DATA (TREAL(1,I),I=1,IWAVE) /
c     VISIBLE
     1      1.38    , 1.37    , 1.365   , 1.36    , 1.357   ,
     2      1.353   , 1.342   , 1.336   , 1.333   , 1.332   ,
     3      1.332   , 1.332   , 1.331   , 1.330   , 1.329   ,
     4      1.328   , 1.327   , 1.323   , 1.319   , 1.313   ,
     5      1.305   , 1.295   , 1.252   , 1.455   , 1.362   , 
     6      1.334   ,
C     INFRARED
     7      1.326   , 1.320   , 1.308   , 1.283   , 1.278   ,
     8      1.313   , 1.326   , 1.310   , 1.293   , 1.270   ,
     9      1.227   , 1.164   , 1.173   , 1.287   , 1.415   ,
     1      1.508   , 1.541   , 1.669   /
C
C     IMAGINARY REFRACTIVE INDEX FOR LIQUID WATER
C
      DATA (TTMAG(1,I),I=1,IWAVE) /
C     VISIBLE
     1      1.e-15  , 1.e-15  , 1.e-15  , 1.e-15  , 1.e-15  ,
     2      1.e-15  , 1.e-15  , 1.e-15  , 1.e-15  , 2.70E-08,
     3      1.68E-07, 6.66E-08, 1.24E-07, 2.94E-07, 8.77E-07,
     4      3.06E-06, 2.05E-06, 2.39E-05, 1.20E-04, 1.18E-04,
     5      6.79E-04, 3.51E-04, 2.39E-03, 0.0442  , 0.00339 ,
     6      0.00833 ,
C     INFRARED
     7                0.0139  , 0.0125  , 0.011   , 0.015   ,
     8      0.075   , 0.086   , 0.039   , 0.035   , 0.035   ,
     9      0.038   , 0.051   , 0.161   , 0.308   , 0.39    ,
     1      0.42    , 0.395   , 0.373   , 0.5     /
c
c     Refractive index for ice taken from Warren,
c     Applied Optics, 23, 1206-1225, 1984, Table 1.
c
c     For the solar bins, Warren's values were first weighted by the
c     solar flux per micron which was taken from Table 131 of the Smithsonian
c     Tables, interpolated to Warren's grid. These values were integrated
c     over the model grid and then normalized by the total solar flux.
c     Weighting the refractive indicies by S probably had only a small effect
c     since S changes only slightly in each bin. Averaging over the entire
c     bin yielded values that in some cases were 10 or 100 times the value
c     obtained by interpolating between the nearest values in Warren's table.
c
c     The refractive indicies for the infrared bins change more slowly within
c     a bin so we have simply interpolated between the nearest values in
c     Warren's table. Of course, no solar weighting is applied.
c
      DATA (TREAL(2,I),I=1,IWAVE) /
c     VISIBLE
     1  1.3480E+00,  1.3407E+00,  1.3353E+00,  1.3307E+00,  1.3275E+00, 
     1  1.3231E+00,  1.3177E+00,  1.3140E+00,  1.3098E+00,  1.3071E+00, 
     1  1.3065E+00,  1.3056E+00,  1.3047E+00,  1.3038E+00,  1.3028E+00, 
     1  1.3014E+00,  1.2997E+00,  1.2955E+00,  1.2906E+00,  1.2837E+00, 
     1  1.2717E+00,  1.2629E+00,  1.1815E+00,  1.4310E+00,  1.3874E+00, 
     1  1.3473E+00,
C     INFRARED
     1               1.3427E+00,  1.3288E+00,  1.3140E+00,  1.2972E+00, 
     1  1.2908E+00,  1.3155E+00,  1.3186E+00,  1.3189E+00,  1.3125E+00, 
     1  1.2848E+00,  1.2268E+00,  1.1921E+00,  1.4906E+00,  1.5588E+00, 
     1  1.5143E+00,  1.4427E+00,  1.3047E+00,  1.4859E+00 /
c
      DATA (TTMAG(2,I),I=1,IWAVE) /
c     VISIBLE
     1  8.0822E-09,  6.7516E-09,  5.7565E-09,  4.8780E-09,  4.2693E-09, 
     1  3.4207E-09,  2.2612E-09,  1.7429E-09,  7.5115E-09,  2.4450E-08, 
     1  3.8345E-08,  7.5497E-08,  1.3765E-07,  2.3306E-07,  5.2670E-07, 
     1  1.7241E-06,  2.2288E-06,  7.3358E-05,  4.8415E-04,  2.6270E-04, 
     1  1.2123E-03,  3.0644E-04,  2.7738E-02,  2.7197E-01,  7.5589E-03, 
     1  1.6397E-02,
C     INFRARED
     1               1.9980E-02,  1.2469E-02,  1.4622E-02,  2.5275E-02, 
     1  5.4226E-02,  6.5501E-02,  5.6727E-02,  5.4499E-02,  4.6806E-02, 
     1  4.0056E-02,  4.9961E-02,  3.0257E-01,  3.6528E-01,  1.7601E-01, 
     1  8.6618E-02,  3.4573E-02,  8.5167E-02,  4.4443E-01 /
C
C     NPROB IS THE NUMBER OF PROBABILITY INTERVALS IN EACH WAVELENGTH
C     INTERVAL. NOTE THAT WAVE BINS 11 AND 12 ARE REVERSED IN ORDER.
C     THIS IS DONE FOR HISTORICAL REASONS.  CROSS SECTIONS, WEIGHTS,
C     REFRACTIVE INDICIES ETC. FOR BINS 11 AND 12 MUST BE REVERSED ALSO.
C
      DATA NPROB   / 1,2,3,4,5,6,7,8,9,10,4*11,3*12,3*13,14,4*15,
     1               16,4*17,4*18,3*19,4*20,4*21,22,20*23,
     2               4*24,25,6*26,3*27,4*28,3*29,3*30,3*31,
     3               4*32,4*33,3*34,3*35,36,4*37,38,12*39,
     4               9*40,3*41,4*42,4*43,3*44               /
C
C     WAVE REFERS TO THE WAVE LENGTHS OF THE CENTERS OF THE SOLAR FLUX
C     BINS AND THE EDGES OF THE INFRARED BINS.
C     FOR THE CURRENT MODEL SETUP, WAVE BINS 11 AND 12 ARE REVERSED
C     IN ORDER, THAT IS 12 PRECEEDS 11. THEREFORE, CROSS SECTIONS,
C     WEIGHTS, REFRACTIVE INDICIES ETC FOR BINS 11 AND 12 IN DATA
C     STATEMENTS MUST BE REVERSED ALSO. THIS IS DONE FOR HISTORICAL
C     REASONS.
C
      DATA WAVE / 0.256, 0.280, 0.296,0.319,0.335,0.365,0.420,0.482,
     1     0.598, 0.690, 0.762, 0.719, 0.813, 0.862, 0.926, 1.005,
     2     1.111, 1.333, 1.562, 1.770, 2.051, 2.210, 2.584, 3.284,
     3     3.809, 4.292,
     4     4.546, 4.878, 5.128, 5.405, 5.714, 6.061, 6.452, 6.897,
     5     7.407, 8.333, 9.009, 10.309, 12.500, 13.889, 16.667,
     6     20.000, 26.316, 35.714, 62.50                         /
C
C     SOLAR FLUXES ( W/M**2)
      DATA SOLFX  /  4.1712E0, 2.5074E0, 1.2024E1, 1.7296E1, 1.2299E1,
     1               5.6975E1, 1.0551E2, 1.3250E2, 2.7804E2, 2.8636E1,
     2               5.9268E1, 5.0747E1, 5.7410E1, 4.3280E1, 7.4598E1,
     3               5.2732E1, 8.6900E1, 1.2126E2, 2.5731E1, 6.0107E1,
     4               1.8400E1, 9.5952E0, 3.5664E1, 1.2764E1, 4.0354E0,
     5               4.5364E0                                         /
C
C     HERE ARE THE CROSS SECTIONS FOR VARIOUS GASES. THERE ARE NWAVE
C     OF THESE. A IS CROSS SECTION, W IS WEIGHT, PS IS PRESSURE SCALE.
C
C     ***********************
C     *  DATA FOR SOLAR     *
C     ***********************
C
C     UNITS ARE (CM**2)/GM
      DATA (AH2O(I),I=1,77)  /      14*0.0, 0.0000, 0.1965, 9.2460,
     1      0.0000, 0.1765, 9.2460, 0.0000, 0.0000, 0.2939, 0.5311,
     2      113.40, 0.0000, 0.0000, 0.3055, 5.2180, 113.00, 0.0000,
     3      0.3355, 5.5090, 124.90, 3*0.00, 0.0000, 0.3420, 7.1190,
     4      95.740, 4*0.00,   0.00, 4*0.00, 4*.3012,4*5.021,4*63.17,
     5      4*699.1,0.0000, 6.3321, 5.336,  123.4,  7*0.0    /
C
C      UNITS ARE (CM AMAGAT)
       DATA (ACO2(I),I=1,77)  /
     1       34*0.0,  0.0,    0.0035, 0.1849, 4*0.,   0.0,
     2       0.0032,  0.0705, 1.2980, 0.0,    0.0,    0.0077,
     3       0.2027,  3.8160, 0.0000, 0.0077, 0.2027, 3.8160,
     4       0.0000,  0.0077, 0.2027, 3.8160, 0.0000, 0.0077,
     5       0.2027,  3.8160, 0.0000, 0.0077, 0.2027, 3.8160,
     6       5*0.00,  0.0000, 0.0053, 0.0705, 0.6732, 6.2880,
     7       83.900             /
C
C     UNITS ARE (CM AMAGAT)
      DATA (AO2(I),I=1,77)   /
     1      10*0.0,  0.00, 0.00, 0.0001, 0.0022, 63*0.0  /
C
C     UNITS ARE (CM AMAGAT)
      DATA (AO3(I),I=1,77)   /
     1      260.0, 100.9, 11.93, 0.7370, 0.0872, 0.0, 0.0,
     2      0.0, 0.1180, 0.0, 67*0.0   /
C
      DATA (PSO2(I),I=1,77)  /  10*0.0, 4*0.75, 63*0.0      /
C
      DATA (PSO3(I),I=1,77)  /  77*0.0                      /
C
      DATA (PSH2O(I),I=1,77) /
     1      14*0.0, 3*0.54, 3*0.54, 0.0, 4*0.54, 0.0, 4*0.52,
     2      4*0.44, 3*0.00, 4*0.62, 5*0.0, 20*0.60,   4*0.60,
     3      7*0.0                     /
C
      DATA (PSCO2(I),I=1,77) /
     1       41*0.0, 4*0.82, 0.0, 20*0.88, 5*0.0, 6*0.93    /
C
      DATA (WEIGHT(I),I=1,77)  /
     1     10*1.00, 0.8471, 0.1316, 0.0158, 0.0055, 0.7825,
     2     0.19240, 0.0251, 0.8091, 0.1689, 0.0221, 1.0000,
     3     0.41070, 0.4885, 0.0751, 0.0258, 1.0000, 0.6385,
     4     0.29660, 0.0484, 0.0164, 0.2964, 0.3854, 0.2583,
     5     0.05990, 0.9486, 0.0455, 0.0059, 0.2962, 0.3234,
     6     0.25880, 0.1216, 0.5866, 0.3524, 0.0449, 0.0160,
     7     1.00000, 2.2823E-2, 3.8272E-3, 3.1557E-3, 2.7938E-3,
     8     2.3677E-1, 3.9705E-2, 3.2738E-2, 2.8984E-2,
     9     1.8182E-1, 3.0489E-2, 2.5139E-2, 2.2256E-2,
     1     1.7783E-1, 2.9820E-2, 2.4587E-2, 2.1768E-2,
     2     8.0862E-2, 1.3560E-2, 1.1180E-2, 9.8984E-3,
     3     0.3807, 0.3416, 0.2256, 0.0522, 1.0000, 0.3908,
     4     0.1383, 0.1007, 0.0967, 0.0984, 0.1751        /
C
C     ***********************
C     *  DATA FOR INFRARED  *
C     ***********************
C
      DATA (WEIGHT(I),I=78,148) /   0.3710, 0.4400, 0.1890,
     1      0.3070, 0.3980, 0.2050, 0.0900, 0.4030, 0.4150,
     2      0.1820, 0.3770, 0.4330, 0.1900, 0.3100, 0.4580,
     3      0.2320, 0.3160, 0.3980, 0.1850, 0.1010, 0.2030,
     4      0.4010, 0.2720, 0.1240, 0.3820, 0.4240, 0.1940,
     5      0.4120, 0.1800, 0.4080, 1.0000, 0.2220, 0.3400,
     6      0.1400, 0.2980, 1.0000, 0.0117, 0.0311, 0.0615,
     7      0.0458, 0.0176, 0.0468, 0.0927, 0.0689, 0.0487,
     8      0.1292, 0.2558, 0.1903, 0.0242, 0.0633, 0.0665,
     9      0.0589, 0.1541, 0.1620, 0.0740, 0.1936, 0.2035,
     1      0.1570, 0.4100, 0.4330, 0.0810, 0.2080, 0.4090,
     2      0.3020, 0.0990, 0.2200, 0.4040, 0.2770, 0.3670,
     3      0.5870, 0.0460                                  /
C
      DATA (PSCO2(I),I=78,148)  /   71*0.0                  /
C
      DATA (PSH2O(I),I=78,148)  /   71*0.0                  /
C
      DATA (PSO3(I),I=78,148)   /   71*0.0                  /
C
      DATA (PSO2(I),I=78,148)   /   71*0.0                  /
C
C     GAUSS ANGLES AND GAUSS WEIGHTS FOR GAUSSIAN INTEGRATION
C     MOMENTS (USE FIRST MOMENT VALUES) N=3
C
      DATA GANGLE  /  0.2123405382, 0.5905331356,
     1                0.9114120405                       /
      DATA GRATIO/0.4679139346, 0.3607615730, 0.1713244924/
      DATA GWEIGHT /  0.0698269799, 0.2292411064,
     1                0.2009319137                       /
C
C     GAUSS ANGLES AND WEIGHTS FOR GAUSSIAN INTEGRATION MOMENTS
C     (USE FIRST MOMENT ONLY)  N=8
C
C      DATA GANGLE  /  0.0446339553, 0.1443662570,
C                      0.2868247571, 0.4548133152, 0.6280678354,
C                      0.7856915206, 0.9086763921, 0.9822200849  /
C
C      DATA GWEIGHT /  0.0032951914, 0.0178429027,
C                      0.0454393195, 0.0791995995, 0.1060473594,
C                      0.1125057995, 0.0911190236, 0.0445508044  /
C
C     ALOS   - LOCSHMIDT'S NUMBER (#/CM**3)
C     AM     - MOLECULAR WEIGHT OF AIR (G/MOL)
C     AVG    - AVAGODROS' NUMBER (#/MOL)
C     G      - GRAVITY (CM/S**2)
C     PI     - PI
C     RGAS   - UNIVERSAL GAS CONSTANT (ERG / MOL K)
C     SCDAY  - NUMBER OF SECONDS IN ONE DAY (S)
C
      DATA AM     / 28.966       /
      DATA ALOS   / 2.68719E19   /
      DATA AVG    / 6.02252E+23  /
      DATA G      / 980.6        /
      DATA PI     / 3.1415926536 /
      DATA RGAS   / 8.31430E+07  /
      DATA SBK    / 5.6697E-8    /
      DATA SCDAY  / 86400.0      /
C
C     EPSILON - roundoff error precision
C
      DATA EPSILON / ALMOST_ZERO  /
C
C     EXPMAX - LARGEST (NEGATIVE) EXP ARGUMENT
C
      DATA EXPMAX  / POWMAX       /
C
C     THIS ROUTINE ASSUMES THAT PRESSURE DOES NOT VARY WITH
C     LOCATION OR TIME.
C
C     CO2MOL - MOLECULAR WEIGHT OF CO2 (G/MOL)
C     O3MOL  - MOLECULAR WEIGHT OF O3 (G/MOL)
C     O2MOL  - MOLECULAR WEIGHT OF O2 (G/MOL)
C
      DATA CO2MOL  /  44.         /
      DATA O3MOL   /  48.         /
      DATA O2MOL   /  32.         /
C
      DATA (PJ(J),J=1,6)  /  1.0, 0.75, 0.50, 0.25, 0.10, 0.01 /
C     
C     OZONE MASS MIXING RATIO (G/G) AT PRESSURE LEVELS
C     DEFINED BY PJ VECTOR. CURRENT VALUES ARE TAKEN FROM
C     U.S. STANDARD ATMOSPHERE MID-LATITUDE PROFILE
C
      DATA O3MIXP / 5.3E-08, 5.5E-08, 6.5E-08,
     1              2.2E-07, 1.4E-06, 1.14E-05 /
C
      DATA (AKO3(1,I),I=1,6)   / 12103.15537, 12402.63363, 12815.75276,
     1                           13293.58735, 12792.67732, 5263.93060 /
      DATA (AKO3(2,I),I=1,6)   / 4073.95522, 3824.69235, 3430.8699,
     1                           2683.45418, 1705.58659, 294.00564    /
      DATA (AKO3(3,I),I=1,6)   / 256.02950, 246.52991, 229.39518,
     1                           189.37447, 123.25784, 19.16440       /
      DATA (AKO3(4,I),I=1,6)   / 0.02424, 0.02417, 0.02402, 0.02357,
     1                           0.02223, 0.00373                     /
C
      DATA (AKCO2(1,I),I=1,6)  / 21.55629, 19.04724, 15.82831,
     1                           11.24202,  6.87763, 1.93565          /
      DATA (AKCO2(2,I),I=1,6)  / 1.41477, 1.19470, 0.94231, 0.63257,
     1                           0.38015, 0.06712                     /
      DATA (AKCO2(3,I),I=1,6)  / 0.02091, 0.01822, 0.01471, 0.00950,
     1                           0.00443, 0.00784                     /
      DATA (AKCO2(4,I),I=1,6)  / 4285.18491, 4190.41876, 3956.55409,
     1                           3319.57528, 2275.73601, 527.26813    /
      DATA (AKCO2(5,I),I=1,6)  / 490.02297, 420.10700, 331.10330,
     1                           211.02566, 109.50739, 18.39552       /
      DATA (AKCO2(6,I),I=1,6)  / 52.24044, 42.14958, 31.06085,
     1                           18.17948, 8.80769, 1.40437           /
C
      DATA (AKH2O(1,I),I=1,6)  / 0.02080, 0.01550, 0.01040, 0.00510,
     1                           0.00200, 0.00062       /
      DATA (AKH2O(2,I),I=1,6)  / 0.17070, 0.12910, 0.08630, 0.04300,
     1                           0.01710, 0.00043       /
      DATA (AKH2O(3,I),I=1,6)  / 3.69350, 3.03920, 2.18270, 1.18620,
     1                           0.48190, 0.00012       /
      DATA (AKH2O(4,I),I=1,6)  / 0.09550, 0.07170, 0.04775, 0.02390,
     1                           0.00970, 0.00096       /
      DATA (AKH2O(5,I),I=1,6)  / 0.56670, 0.42490, 0.28390, 0.14190,
     1                           0.05620, 0.00565       /
      DATA (AKH2O(6,I),I=1,6)  / 4.42380, 3.37860, 2.27150, 1.14360,
     1                           0.45980, 0.04582       /
      DATA (AKH2O(7,I),I=1,6)  / 52.2782, 46.7368, 37.4543, 22.1443,
     1                           9.48870, 0.96239       /
      DATA (AKH2O(8,I),I=1,6)  / 0.74700, 0.55860, 0.37410, 0.18640,
     1                           0.07480, 0.00757      /
      DATA (AKH2O(9,I),I=1,6)  / 6.21910, 4.71170, 3.17460, 1.57760,
     1                           0.63070, 0.06400      /
      DATA (AKH2O(10,I),I=1,6) / 108.8498, 93.64720, 70.7601, 40.4311,
     1                           16.7743, 1.64448      /
      DATA (AKH2O(11,I),I=1,6) / 5.32480, 4.02140, 2.69230, 1.33950,
     1                           0.53460, 0.05812      /
      DATA (AKH2O(12,I),I=1,6) / 37.4911, 28.7735, 19.3551, 9.8670,
     1                           3.93660, 0.48536      /
      DATA (AKH2O(13,I),I=1,6) / 399.4563, 372.4726, 317.0208, 199.0702,
     1                           91.07190, 9.51138     /
      DATA (AKH2O(14,I),I=1,6) / 19.5458, 14.6823, 9.89610, 4.88830,
     1                           1.96260, 0.19703      /
      DATA (AKH2O(15,I),I=1,6) / 129.8270, 100.9301, 68.2980, 34.9858,
     1                           13.9141, 1.39263      /
      DATA (AKH2O(16,I),I=1,6) / 1289.8795, 1225.3730, 1077.8933,
     1                           713.0905, 338.1374, 34.83768    /
      DATA (AKH2O(17,I),I=1,6) / 10.1830, 7.64440, 5.09240, 2.54900,
     1                           1.04120, 0.10272      /
      DATA (AKH2O(18,I),I=1,6) / 61.2038, 46.4646, 31.2440, 15.7152,
     1                           6.24550, 0.62696      /
C
      DATA (AKH2O(19,I),I=1,6) / 363.8096, 297.8458, 214.1972, 112.6198,
     1                           45.6958, 4.55650      /
      DATA (AKH2O(20,I),I=1,6) / 2005.95, 2053.91, 1991.06, 1561.74,
     1                           825.730, 89.7681      /
      DATA (AKH2O(21,I),I=1,6) / 15.3380, 11.5228, 8.5296, 3.8592,
     1                           1.54910, 0.15449      /
      DATA (AKH2O(22,I),I=1,6) / 72.4480, 54.6480, 32.2028, 18.2840,
     1                           7.31640, 0.73264      /
      DATA (AKH2O(23,I),I=1,6) / 440.2140, 352.1620, 512.833, 128.293,
     1                           51.6640, 5.14812      /
      DATA (AKH2O(24,I),I=1,6) / 3343.51, 3347.37, 512.833, 233.720,
     1                           1182.99, 124.7161     /
      DATA (AKH2O(25,I),I=1,6) / 3.76110, 2.82540, 1.88200, 0.943200,
     1                           0.37750, 0.037776     /
      DATA (AKH2O(26,I),I=1,6) / 29.5500, 22.4360, 15.1239, 7.57360,
     1                           3.02760, 0.30409      /
      DATA (AKH2O(27,I),I=1,6) / 424.620, 379.480, 304.175, 180.650,
     1                           78.0440, 7.94202      /
      DATA (AKH2O(28,I),I=1,6) / 1.25230, 0.93840, 0.62750, 0.31300,
     1                           0.12450, 0.01251      /
      DATA (AKH2O(29,I),I=1,6) / 27.4300, 22.1200, 15.7100, 8.2700,
     1                           3.37560, 0.341119     /
      DATA (AKH2O(30,I),I=1,6) / 0.14650, 0.11010, 0.07330, 0.03670,
     1                           0.01480, 0.014761     /
C
      DATA (AKH2O(31,I),I=1,6) / 6*4.95e-8                 /
      DATA (AKH2O(32,I),I=1,6) / 6*6.06e-8                 /
      DATA (AKH2O(33,I),I=1,6) / 6*1.12e-7                 /
C
      DATA (AKH2O(34,I),I=1,6) / 1.32700, 1.38470, 0.97515, 0.50750,
     1                           0.20570, 0.02090      /
      DATA (AKH2O(35,I),I=1,6) / 0.07038, 0.05380, 0.03595, 0.01790,
     1                           0.00720, 0.00068      /
      DATA (AKH2O(36,I),I=1,6) / 0.00782, 0.00590, 0.00393, 0.00190,
     1                           0.00078, 0.00002      /
      DATA (AKH2O(37,I),I=1,6) / 0.00128, 0.00096, 0.00642, 0.00032,
     1                           0.00013, 0.00010      /
      DATA (AKH2O(38,I),I=1,6) / 3.00300, 2.91330, 1.97760, 1.00080,
     1                           0.40170, 0.0402       /
      DATA (AKH2O(39,I),I=1,6) / 0.13090, 0.09870, 0.06580, 0.03290,
     1                           0.00131, 0.00131      /
      DATA (AKH2O(40,I),I=1,6) / 0.01379, 0.01350, 0.00690, 0.00350,
     1                           0.00014, 0.00014      /
C
      DATA (AKH2O(41,I),I=1,6) / 29.2600, 27.7800, 19.0800, 9.7400,
     1                           3.92000, 0.39310      /
      DATA (AKH2O(42,I),I=1,6) / 1.27900, 0.96860, 0.64580, 0.32290,
     1                           0.12910, 0.01290      /
      DATA (AKH2O(43,I),I=1,6) / 0.13480, 0.10150, 0.06770, 0.03380,
     1                           0.01350, 0.00135      /
      DATA (AKH2O(44,I),I=1,6) / 716.500, 666.670, 491.250, 266.300,
     1                           109.670, 11.2050      /
      DATA (AKH2O(45,I),I=1,6) / 39.7800, 29.7300, 19.8800, 9.9500,
     1                           3.98000, 0.37610      /
      DATA (AKH2O(46,I),I=1,6) / 4.47600, 3.35300, 2.23640, 1.11900,
     1                           0.44740, 0.00996      /
      DATA (AKH2O(47,I),I=1,6) / 0.73960, 0.55410, 0.36940, 0.18460,
     1                           0.07390, 0.06010      /
      DATA (AKH2O(48,I),I=1,6) / 5631.00, 6090.99, 4475.54, 2422.79,
     1                           998.340, 100.540       /
      DATA (AKH2O(49,I),I=1,6) / 414.800, 309.620, 206.840, 103.450,
     1                           41.3600, 4.13620       /
      DATA (AKH2O(50,I),I=1,6) / 53.4600, 40.2600, 26.8600, 13.4300,
     1                           5.37520, 0.53760       /
      DATA (AKH2O(51,I),I=1,6) / 9.56300, 7.17130, 4.78020, 2.38970,
     1                           0.95580, 0.09560       /
      DATA (AKH2O(52,I),I=1,6) / 2779.00, 2844.55, 1937.23, 982.710,
     1                           394.850, 39.5000       /
      DATA (AKH2O(53,I),I=1,6) / 134.600, 108.410, 72.2400, 36.1100,
     1                           14.4000, 1.44600       /
      DATA (AKH2O(54,I),I=1,6) / 0.008995, 8.62480, 5.75790, 2.8800,
     1                           1.15300, 0.10640       /
C     
C     CORERAD - RADIUS OF CORE OF AEROSOL PARTICLES
C     COREREAL- REAL PART OF REFRACTIVE INDEX OF CORES
C     COREIMAG- IMAGINARY PART OF REFRACTIVE INDEX OF CORES
C
      DATA CORERAD  / 0.0        /
      DATA COREREAL / 1.25       /
      DATA COREIMAG / 0.5        /
c
c     Bergstrom's water vapor continuum fix
c 
      data (contnm(i),i=1,71) / 36*0., 12*2.e-7,
     1      9*4.266e-7, 3*9.65e-7, 4*2.423e-6,
     2      4*5.765e-6, 3*1.482e-5 /
C
C     NSOL    - NUMBER OF SOLAR WAVELENGTHS
C     NIR     - NUMBER OF INFRARED WAVELENGTHS
C     NWAVE   - TOTAL NUMBER OF SOLAR+INFRARED WAVELENGTHS
C     NSOLP   - SOLAR WAVELENGTHS ON FINER GRID
C     NIRP    - INFRARED WAVELENGTHS ON FINER GRID
C     NTOTAL  - NUMBER OF PROBABILITY INTERVALS
C     NGAUSS  - NUMBER OF GAUSS QUADRATURE POINTS
C
      data NSOL     /   ISOL     /
      data NIR      /   IIR      /
      data NWAVE    /   IWAVE    /
      data NSOLP    /   ISOLP    /
      data NIRP     /   IIRP     /
      data NTOTAL   /   ITOTAL   /
      data NGAUSS   /   IGAUSS   /
C
C     DERIVED PARAMETERS
C
      SQ3     =   SQRT(3.)
      JDBLE   =   2*NLAYER
      JN      =   JDBLE-1
      TPI     =   2.*PI
      CPCON   =   1.006
      FDEGDAY =   1.0E-4*G*SCDAY/CPCON
c
c     Initialize the aerosol concentrations
c
      do ig = 1, NGROUP
        do j = 1, NLAYER
          do i = 1, NRAD
            caer(i,j,ig) = 0.
          enddo
        enddo
      enddo
c
c     Inverse of nprob matrix:  NPROBI(i,1) is first probability interval
c     for wavelength interval i, and NPROBI(i,2) is number of probability
c     intervals in the wavelength interval
c
      do i = 1, nwave
        kount = 0
        do j = 1, ntotal
          if (nprob(j).eq.i) then
            if (kount.eq.0) nprobi(i,1) = j
            kount = kount + 1
          endif
        enddo
        nprobi(i,2) = kount
      enddo
c
c     Load wavelengths into interface common block
c
      do i = 1, NWAVE+1
        wave_aerad(i) = wave(i)
      enddo
c
c     Get scalars from interface common block:
c
c     ISL        - do solar calculations when = 1
c     IR         - do infrared calculations when = 1
C     UO         - SOLAR ZENITH ANGLE
C     ALBEDO_SFC - SURFACE ALBEDO
C     EMISIR     - SURFACE IR EMISSIVITY
C     PTOP       - PRESSURE AT TOP OF MODEL (DYNES/CM**2)
C     PBOT       - PRESSURE AT BOTTOM OF MODEL (DYNES/CM**2)
C     TGRND      - GROUND TEMPERATURE, DEG K
c
      ISL        = isl_aerad
      IR         = ir_aerad
      U0         = u0_aerad
      ALBEDO_SFC = sfc_alb_aerad
      EMISIR     = emisir_aerad
      PTOP       = ptop_aerad
      PBOT       = pbot_aerad
      TGRND      = tsfc_aerad
c
c     Get atmospheric temperature and pressure profiles
c     from interface common block
C
C     P - PRESSURE IN DYNES/CM**2
C     T - TEMPERATURE IN DEG K
C
      do k = 1,nvert
        p(k) = p_aerad(k)
        t(k) = t_aerad(k)
      enddo
C
C     PBAR  - LAYER AVERAGE PRESSURE (BARS)
C     (NOTE - THE TOP LAYER IS FROM PTOP TO 0, SO AVERAGE = PTOP/2)
C     PRESS - PRESSURE AT EDGE OF LAYER (dyne/cm^2)
C     DPG   - MASS OF LAYER (G / CM**2)
C
      PBAR(1)  = PTOP/2.0E6
      PRESS(1) = PTOP
      DO 45 J  = 2,NVERT
         PBAR(J)  = P(J-1)/1.0E6
         PRESS(J) = (P(J-1) + P(J)) * 0.5
         DPG(J-1) = (PRESS(J)-PRESS(J-1)) / G
 45   CONTINUE
      PBAR(NLAYER)  = P(NVERT)/1.0E6
      PRESS(NLAYER) = PBOT
      DPG(NVERT)    = (PRESS(NLAYER)-PRESS(NVERT)) / G
C
C     SKIN TEMPERATURE
C
      TT(NLAYER) = TGRND
C
C     AMOUNT OF WATER VAPOR ABOVE MODEL DOMAIN (GM / CM**2)
c     From 1976 U.S. Standard Atmosphere, mid-latitude sounding:
c
c     For z_top = 5 km
c     RDH2O(1) = .13
c
c     For z_top = 2.6 km
c     RDH2O(1) = .64
C
      RDH2O(1)   = h2ocol_aerad
C
C     INTERPOLATE TEMPERATURE FROM LAYER CENTER (T) TO LAYER EDGE (TT)
C
      TT(1) = T(1)
      DO 12 J = 2, NVERT
         TT(J) = T(J-1) * (PRESS(J)/P(J-1)) **
     1           (log(T(J)/T(J-1))/log(P(J)/P(J-1)))
   12 CONTINUE
C
C     DEFINE MASS MIXING RATIOS. O3MIX TAKEN FROM U.S. STANDARD ATMOS-
C     PHERE, MID-LATITUDE SOUNDING
C
      O2MIX          =   0.22*O2MOL/AM
      CO2MIX         =   3.5E-4*CO2MOL/AM
C
C     OZONE COLUMN ABUNDANCE O3C (#/CM**2) ABOVE PTOP WAS CALCULATED
C     FROM THE U.S. STANDARD ATMOSPHERE MID-LATITUDE PROFILE.
C
c    This is for z_top = 5 km
c     O3C            =   9.02E18
c    This is for z_top = 3 km
c     O3C            =   9.2E18
c    This is for z_top = 1 km
      O3C            =   9.3E18
C
C     CONVERT O3C TO MASS MIXING RATIO O3MIX2.
C
      O3MIX2 = O3C*O3MOL*G/(PTOP*AVG)
C
C     CALCULATE CLOUD AND AEROSOL OPACITY.
C
      DO 100  J             = 1,NLAYER
          DO 50 L           = 1,NTOTAL
               TAUAER(L,J)  = 0.
               TAUCLD(L,J)  = 0.
               WOL(L,J)     = 0.
               GOL(L,J)     = 0.
               GCLD(L,J)    = 0.
               TAURAY(L)    = 0.
 50        CONTINUE
100   CONTINUE
C
      DO 120 L           =   NSOLP+1,NTOTAL
         LTEMP(L-NSOLP)  =   NPROB(L) - NSOL
 120  CONTINUE
C
      X                  =   ALOS/AVG
C
C     CONVERT SOLAR ABSORPTION COEFFICIENTS TO CM**2/GM.
C
      DO 140 L           =   1,NSOLP
         ACO2(L)         =   ACO2(L)/(X*CO2MOL)
         AO2(L)          =   AO2(L)/(X*O2MOL)
         AO3(L)          =   AO3(L)/(X*O3MOL)
 140  CONTINUE
C
C     CALCULATE ABSORPTION COEFFICIENTS
C
      DO 160 J           =   1,NLAYER
          DO 150 L       =   1,NSOLP
            PAH2O(L,J)   =   AH2O(L)*PBAR(J)**PSH2O(L)
            PACO2(L,J)   =   ACO2(L)*PBAR(J)**PSCO2(L)
            PAO2(L,J)    =   AO2(L)*PBAR(J)**PSO2(L)
            PAO3(L,J)    =   AO3(L)*PBAR(J)**PSO3(L)
 150      CONTINUE
 160  CONTINUE
C
      DO 180 J           =   1,NLAYER
         DO 170 L        =   NSOLP+1,NTOTAL
           PAO2(L,J)     =   0.0
           PACO2(L,J)    =   0.0
           PAO3(L,J)     =   0.0
           PAH2O(L,J)    =   0.0
 170     CONTINUE
 180  CONTINUE
C
      PS                 =   0.0
      DO 300 J           =   1,NLAYER
          DO 200 I       =   1,6
            II           =   I
            IF(PBAR(J) .GT. PJ(I)) GO TO 202
 200      CONTINUE
 202      CONTINUE
C
            if( ii .eq. 1 ) ii = 2
            DP           =   log(PJ(II-1)/PJ(II))
            if( pbar(j) .gt. pj(6) )then
              ij = ii - 1
            else
              ij = ii
            endif
            PS           =   PBAR(J)/PJ(IJ)
C
            IF (J.NE.1)
     1      O3MIX(J) =O3MIXP(IJ)*PS**(log(O3MIXP(II-1)/O3MIXP(II))/DP)
C
            DO 210 L           =   1,31
              PAH2O(NSOLP+L,J) = AKH2O(L,IJ)*PS**(log
     1                           (AKH2O(L,II-1)/AKH2O(L,II))/DP)
 210        CONTINUE
C
            DO 220 L           =   32,35
              PAH2O(NSOLP+L,J) = AKH2O(32,IJ)*PS**(log
     1                           (AKH2O(32,II-1)/AKH2O(32,II))/DP)
              PAO3(NSOLP+L,J)  = AKO3(L-31,IJ)*PS**(log
     1                           (AKO3(L-31,II-1)/AKO3(L-31,II))/DP)
 220        CONTINUE
C
            PAH2O(NSOLP+36,J)  = AKH2O(33,IJ)*PS**(log
     1                           (AKH2O(33,II-1)/AKH2O(33,II))/DP)
C
            DO 230 L           =   37,40
              PACO2(NSOLP+L,J) = AKCO2(1,IJ)*PS**(log
     1                           (AKCO2(1,II-1)/AKCO2(1,II))/DP)
              PAH2O(NSOLP+L,J) = AKH2O(L-3,IJ)*PS**(log
     1                           (AKH2O(L-3,II-1)/AKH2O(L-3,II))/DP)
 230        CONTINUE
C
            DO 240 L           =   41,44
             PACO2(NSOLP+L,J)  =   AKCO2(2,IJ)*PS**(log
     1                             (AKCO2(2,II-1)/AKCO2(2,II))/DP)
             PAH2O(NSOLP+L,J)  =   AKH2O(L-7,IJ)*PS**(log
     1                             (AKH2O(L-7,II-1)/AKH2O(L-7,II))/DP)
 240        CONTINUE
C
            DO 250 L           =   45,48
              PACO2(NSOLP+L,J) = AKCO2(3,IJ)*PS**(log
     1                           (AKCO2(3,II-1)/AKCO2(3,II))/DP)
              PAH2O(NSOLP+L,J) = AKH2O(L-11,IJ)*PS**(log
     1                           (AKH2O(L-11,II-1)/AKH2O(L-11,II))/DP)
 250        CONTINUE
C
            DO 260 L           =   49,51
              PACO2(NSOLP+L,J) = AKCO2(4,IJ)*PS**(log
     1                           (AKCO2(4,II-1)/AKCO2(4,II))/DP)
              PAH2O(NSOLP+L,J) = AKH2O(L-11,IJ)*PS**(log
     1                           (AKH2O(L-11,II-1)/AKH2O(L-11,II))/DP)
 260        CONTINUE
C
            DO 270 L           =   52,54
              PACO2(NSOLP+L,J) = AKCO2(5,IJ)*PS**(log
     1                           (AKCO2(5,II-1)/AKCO2(5,II))/DP)
              PAH2O(NSOLP+L,J) = AKH2O(L-14,IJ)*PS**(log
     1                           (AKH2O(L-14,II-1)/AKH2O(L-14,II))/DP)
 270        CONTINUE
C
            DO 280 L           =   55,57
              PACO2(NSOLP+L,J) = AKCO2(6,IJ)*PS**(log
     1                           (AKCO2(6,II-1)/AKCO2(6,II))/DP)
              PAH2O(NSOLP+L,J) = AKH2O(L-17,IJ)*PS**(log
     1                           (AKH2O(L-17,II-1)/AKH2O(L-17,II))/DP)
 280        CONTINUE
C
            DO 290 L           =   58,NIRP
              PAH2O(NSOLP+L,J) = AKH2O(L-17,IJ)*PS**(log
     1                           (AKH2O(L-17,II-1)/AKH2O(L-17,II))/DP)
 290  CONTINUE
 300  CONTINUE
C
C     STORE O3MIX2 IN O3MIX(1)
C
      O3MIX(1) = O3MIX2
C
C     HERE WE FIND TAUGAS. IT IS TAUCO2+TAUO2+TAUO3.
C
      PM             =   PTOP/G
      DO 305 L       =   1,NTOTAL
         TAUGAS(L,1) =   PM*(O2MIX*PAO2(L,1)+CO2MIX*PACO2(L,1)+
     $                         O3MIX(1)*PAO3(L,1))

 305  CONTINUE
C
      DO 308   J     =   2,NLAYER
         PM          =   DPG(J-1)
         DO 306 L    =   1,NTOTAL
            TAUGAS(L,J) = PM*( O2MIX*PAO2(L,J)+CO2MIX*PACO2(L,J)+
     $                         O3MIX(J)*PAO3(L,J))

 306  CONTINUE
 308  CONTINUE
C
C     WAVE MUST BE IN MICRONS
C     CALCULATE RAYLEIGH OPTICAL DEPTH PARAMETERS.
C
      DO 310 L      =   1,NTOTAL
C
          WVO       =    WAVE(NPROB(L))
          TAURAY(L) =   (8.46E-9/WVO**4) *
     1                  ( 1.+0.0113/WVO**2+0.00013/WVO**4 )
 310  CONTINUE
C
C     WE DO NOT INCLUDE RAYLEIGHT SCATTERING IN INFRARED
C
      DO 330 J          =   1,NVERT
         DO 320 L       =   1,NTOTAL
           if( L .LE. NSOLP ) then
             PARAY(L,J+1) = TAURAY(L)*DPG(J)*G
           else
             PARAY(L,J+1) = 0.
           endif
 320     CONTINUE
C
         DO 325 L       =   1,NTOTAL
           if( L .LE. NSOLP ) then
             PARAY(L,1) = TAURAY(L)*PTOP
           else
             PARAY(L,1) = 0.
           endif
 325     CONTINUE
C
 330  CONTINUE
C
C **********************************************************************
C
C            CALCULATE THE AEROSOL EXTINCTION CROSS SECTIONS
C
C **********************************************************************
C
c     Get <is_grp_ice> and radius grid from interface common block
c     and calculate cross-sectional area for each bin.
c
      do ig = 1, NGROUP
        is_grp_ice(ig) = is_grp_ice_aerad(ig)
        DO I = 1, NRAD
           R(I,ig)      = r_aerad(i,ig)
           XSECTA(I,ig) = PI * R(I,ig)**2.
           RUP(I,ig)    = rup_aerad(i,ig)
        ENDDO
      enddo
c
c     Set <i_mie> = I_READ to write the mie coefficients to a data file,
c                 = I_WRITE to read them
c
      i_mie = I_READ
      i_mie = I_WRITE

      if ( i_mie .eq. I_READ ) then
c
c     Read extinction and scattering coefficients
c
        open(LUNMIE,file='../data/mie.data',form='formatted')
c
c     Check that input file is consistent with radius and
c     wavelength grids
c
        read(LUNMIE,*) mwave,mrad,mgroup
        all_ok = (mwave .eq. NWAVE) .and. (mrad .eq. NRAD) .and.
     1           (mgroup .eq. NGROUP) 
c
        read(LUNMIE,*) (rmin(i),i=1,mgroup)
c
        do ig = 1, NGROUP
c...dbg:
c          all_ok = all_ok .and. ( r(1,ig) .eq. rmin(ig) )
        enddo
c
        if ( .not. all_ok )then
          write(LUNOPRT,*) ' SETUPRAD: mie.data grid(s) bad: '
          write(LUNOPRT,*) ' In mie.data, mwave, mrad, mgroup = ',
     1           mwave, mrad, mgroup
          write(LUNOPRT,*) ' and rmin = ', rmin
          stop 1
        endif
c
        do ig = 1,NGROUP
         do i = 1,NRAD
          do L = 1,NWAVE
            read(LUNMIE,*) rdqext(i,ig,L), qscat(i,ig,L), qbrqs(i,ig,L)
          enddo
         enddo
        enddo
c
        close(LUNMIE)
c
      else
c
c     Calculate extinction and scattering coefficients
c
      do ig = 1,NGROUP
c
c     Select ice/liquid index of refractive index array
c
      if( is_grp_ice(ig) )then
        irefr = 2
      else
        irefr = 1
      endif
c
c     <thetd> is angle between incident and scattered radiation
c     <j_thetd> is number of <thetd> values to consider
c
      thetd = 0.0
      n_thetd = 1

      DO 110 L           =   1,NWAVE
C
       REAL = TREAL(irefr,L)
       TMAG = TTMAG(irefr,L)
C
C      CALCULATE THE CENTER OF THE WAVELENGTH INTERVAL OF AN IR INTERVAL
C
        if( L .le. NSOL ) then
          awave = wave(L)
        else
          awave = 0.5*(WAVE(L)+WAVE(L+1))
        endif
       WVNO              =   2.*PI/(AWAVE*1.0E-4)
C
       DO 107 I          =   1,NRAD
        IF(I .EQ. 1) THEN
          DDR            =   (RUP(1,ig)-R(1,ig))/5.
          RR             =   R(1,ig)
        ELSE
          DDR            =   (RUP(I,ig)-RUP(I-1,ig))/5.
          RR             =   RUP(I-1,ig)
        ENDIF
C
        RDQEXT(I,ig,L)   =   0.0
        QSCAT(I,ig,L)    =   0.0
        QBRQS(I,ig,L)    =   0.0
C
        DO 104 J         =   1,6
C 
c
c  Limit x=2*pi/wave to no larger 100 to avoid anguish in the Mie
c  code.
c
         if( wvno*rr .gt. 100. ) rr = 100./wvno

         CALL DMIESS(RR,REAL,TMAG,thetd,n_thetd,QEXTD,QSCATD,CTBRQS,
     1               CORERAD,COREREAL,COREIMAG,WVNO)
C 
         RDQEXT(I,ig,L)     =   RDQEXT(I,ig,L)+QEXTD/6.
         QSCAT(I,ig,L)      =   QSCAT(I,ig,L)+QSCATD/6.
         QBRQS(I,ig,L)      =   QBRQS(I,ig,L)+CTBRQS/6.
         RR                 =   RR+DDR
C 
 104    CONTINUE

 107   CONTINUE
 110  CONTINUE
C
      enddo      ! ig=1,NGROUP
c
      endif   
c
      if ( i_mie .eq. I_WRITE ) then

c     Write extinction and scattering coefficients to data file
c
        open(LUNMIE,file='../data/mie.data',form='formatted')
c
        write(LUNMIE,*) NWAVE,NRAD,NGROUP
        write(LUNMIE,*) (r(1,ig),ig=1,NGROUP)
c
        do ig = 1,NGROUP
         do i = 1,NRAD
          do L = 1,NWAVE
            write(LUNMIE,*) rdqext(i,ig,L), qscat(i,ig,L), qbrqs(i,ig,L)
          enddo
         enddo
        enddo
c
        close(LUNMIE)
c
      endif   
c
c     Write some values to print file
c
      do ig = 1, NGROUP
        WRITE (LUNOPRT,500) ig
        DO I = 1, NRAD
          sizparm6 = 2.*pi*r(i,ig)/(wave(6)*1.e-4)
          sizparm24 = 2.*pi*r(i,ig)/(wave(24)*1.e-4)
          WRITE(LUNOPRT,505) I,R(I,ig),rdqext(i,ig,6),sizparm6,
     1       rdqext(i,ig,24),sizparm24
        ENDDO
      enddo
C
 500  FORMAT(/," SETUPRAD: igroup = ",i4,//,
     1       "   i     r(cm)   rdqext(6)      x6",
     2       "   rdqext(24)      x24",/)
 505  FORMAT(I4,5(1PE11.2))
C
C *********************************************************************
C
C                              CHECK SUM OF WEIGHTS
C
C **********************************************************************
C
      SUM           =   0.0
      SUM1          =   0.0
      SUM2          =   0.0
      DO 340 L      =   1,NSOLP
         SUM        =   SUM+WEIGHT(L)
 340  CONTINUE
      DO 350 L      =   NSOLP+1,NTOTAL
         SUM1       =   SUM1+WEIGHT(L)
 350  CONTINUE
      SUM2          =   SUM+SUM1
C
      IF ( ABS(NWAVE-SUM2) .GT. 1.E-3 ) WRITE(LUNOPRT,355) SUM,SUM1,SUM2
C
 355  FORMAT(//,"SETUPRAD: Error in weights ",/,
     1        " Sum of weights for solar = ",1PE15.5,/,
     2        " sum of weights for ir = ",1PE15.5,/,
     3        " total sum = ",1PE15.5)
C
      DO 360 L   =   1,NSOLP
         SOL(L)  =   SOLFX(NPROB(L)) * WEIGHT(L)
 360  CONTINUE
C
C *********************************************************************
C
C     COMPUTE PLANCK FUNCTION TABLE. WAVE IS IN UNITS OF MICRONS.
C
C **********************************************************************
C
      NLOW          =   12500
      NHIGH         =   32500
      NCOUNT        =   NHIGH - NLOW
c
c     Set <blackbody_above> = 1 to include a source of radiation
c     at the top of the radiative transfer model domain
c    
      blackbody_above = 0
      t_above = 150.
c
c     Set <beyond_spectrum> = 1 to include blackbody radiation at
c     wavelengths longer than WAVE(NWAVE+1) in PLANK(NWAVE+1-NSOL) and
c     at wavelengths shorter than WAVE(NSOL+1) in PLANK(1)
c
      beyond_spectrum = 1
c
      DO 380 J  =   1,NCOUNT
         JJ     =   NLOW+J
         T1     =   0.01 * FLOAT(JJ)

         if( beyond_spectrum .eq. 1 )then

           plank(nwave+1-nsol,j) = t1**4
           DO I =   NSOL+2,NWAVE
              K =   I-NSOL
              V =   1.438E4  /  WAVE(I)
              CALL PLNK(V,T1,PLANK(K,J))
           ENDDO

         else

           DO I =   NSOL+1,NWAVE+1
              K =   I-NSOL
              V =   1.438E4  /  WAVE(I)
              CALL PLNK(V,T1,PLANK(K,J))
           ENDDO

         endif

 380  CONTINUE
C
      DO 410 J   =   1,NCOUNT

         if( beyond_spectrum .eq. 1 )then

           plank(1,j) = plank(2,j)*sbk/pi
           DO L  =   NSOL+2,NWAVE
              K  =   L-NSOL
              PLANK(K,J) = (PLANK(K+1,J)-PLANK(K,J))*SBK/PI
           ENDDO

         else

           DO L  =   NSOL+1,NWAVE
              K  =   L-NSOL
              PLANK(K,J) = (PLANK(K+1,J)-PLANK(K,J))*SBK/PI
           ENDDO

         endif

 410  CONTINUE
C
      RETURN
      END
