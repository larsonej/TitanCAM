c
c  Declare and define symbolic constants
c   (Implicit typing in effect unless explicitly specified)
c
c
c  Include implicit declarations
c
      include 'precision.h'
c
c
c  Include symbolic constants shared between aerosol and radiation models
c
      include 'aerad.h'
c
c
C DEFINE THE DIMENSIONS USED BY THE RADIATION MODEL
C                                                                       
C NVERT  = MAXIMUM NUMBER OF LAYERS;                                    
C NLAYER = MAXIMUM NUMBER OF LAYER BOUNDARIES                           
C NDBL   = TWICE THE MAXIMUM NUMBER OF LAYER BOUNDARIES                 
C NRAD   = MAXIMUM NUMBER OF AEROSOL RADIUS BINS;                       
C                                                                       
      PARAMETER ( NVERT = NZ_RAD )                                              
      PARAMETER ( NRAD = NBIN )                                              
      PARAMETER ( NLAYER = NVERT+1 )
      PARAMETER ( NDBL = 2*NLAYER )                         
      PARAMETER ( NRADVER = NRAD*NVERT )                                    
      PARAMETER ( NRADLAY = NRAD*NLAYER )                                   
C                                                                       
C NTOTAL = TOTAL NUMBER OF PROBABILITY INTERVALS;                       
C NSOLP  = NUMBER OF SOLAR PROBABILITY INTERVALS;                       
C NIRP   = NUMBER OF INFRARED PROBABILITY INTERVALS;                    
C                                                                       
      PARAMETER ( NSOLP = 77 )
      PARAMETER ( NIRP = 71 )
      PARAMETER ( NTOTAL = NSOLP + NIRP )          
C                                                                       
C NGAUSS = TOTAL NUMBER OF GAUSS QUADRATURE POINTS;                     
C                                                                       
      PARAMETER ( NGAUSS = 3)                                           
C                                                                       
C NCOUNT = USED TO CALCULATE PLANK FUNCTION.                            
C                                                                       
      PARAMETER ( NLOW = 12500 )
      PARAMETER ( NHIGH = 32500 )
      PARAMETER ( NCOUNT = NHIGH - NLOW )
c
c Define values of flag used for reading/writing Mie coefficients
c
       parameter ( I_READ = 0 )
       parameter ( I_WRITE = 1 )
C                                                                       
C CONSTANT PARAMETERS THAT MIGHT BE SPECIFIED BY AN EXTERNAL MODEL      
C                                                                       
       logical is_grp_ice

      COMMON /irad1/                                                  
     3 O3MIX(NLAYER), O3MIXP(6), O3C, VRAT,                         
     2 PTOP, PBOT, RMIN(NGROUP), R(NRAD,NGROUP),
     3 is_grp_ice(NGROUP)
C                                                                       
C TIME-DEPENDENT VARIABLES THAT MIGHT BE SPECIFIED BY AN EXTERNAL MODEL 
C                                                                       
      COMMON /irad2/                                                  
     1 U0EXT,                                                       
     2 ALBEDO_SFC, EMISIR, 
     3 P(NVERT), T(NVERT), Q(NVERT)
C                                                                       
C OUTPUT VARIABLES, CALCULATED BY THE RADIATION MODEL, THAT MIGHT       
C BE USED BY THE EXTERNAL MODEL                                         
C                                                                       
      COMMON /irad3/                                                  
     1 HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER),                      
     2 SOLNET, XIRDOWN, XIRUP                                           
C                                                                       
C INITIATED IN SETUPRAD FOR RADIATION CALCULATION                       
C
      COMMON /irad4/                                                  
     1  LLA, LLS, JDBLE, JN, EPSILON, EXPMAX,                           
     2  TPI, SQ3, SBK,                                                  
     2  AM, AVG, ALOS, G, PI, SCDAY, RGAS,                              
     3  GANGLE(NGAUSS), GWEIGHT(NGAUSS), SFLX(NSOL), WVLN(NSOL),        
     4  GRATIO(NGAUSS),
     2  EMIS(NTOTAL), RSFX(NTOTAL),LTEMP(NTOTAL),NPROB(NTOTAL),         
     3  SOL(NTOTAL),TAURAY(NTOTAL),WEIGHT(NTOTAL),                      
     4  GCLD(  NTOTAL,NLAYER),   GOL(NTOTAL,NLAYER),                    
     5  PARAY( NTOTAL,NLAYER),TAUAER(NTOTAL,NLAYER),                    
     6  WCLD(NTOTAL,NLAYER),                                            
     6  TAUCLD(NTOTAL,NLAYER),WOL(   NTOTAL,NLAYER),                    
     7  TREAL(2,NWAVE), TTMAG(2,NWAVE),
     8  contnm(NIRP), nprobi(NWAVE,2)
C                                                                       
      COMMON/irad5/                                                  
     1       ACO2(NTOTAL),           AH2O(NTOTAL),                      
     2       AO2(NTOTAL),            AO3(NTOTAL),                       
     6       PACO2(NTOTAL,NLAYER),   PAH2O(NTOTAL,NLAYER),              
     7       PAO2(NTOTAL,NLAYER),    PAO3(NTOTAL,NLAYER),               
     8       PLANK(NIR+1,NCOUNT),                                       
     9       PSCO2(NTOTAL),          PSH2O(NTOTAL),                     
     1       PSO2(NTOTAL),           PSO3(NTOTAL),                      
     2       SOLFX(NSOL),            WAVE(NWAVE+1),                     
     3       TAUGAS(NTOTAL,NLAYER),  
     4       XSECTA(NRAD,NGROUP),    RUP(NRAD,NGROUP),                         
     5       QSCAT(NRAD,NGROUP,NWAVE), blackbody_above,
     6       QBRQS(NRAD,NGROUP,NWAVE), t_above,     
     7       RDQEXT(NRAD,NGROUP,NWAVE)                                         
C                                                                       
      COMMON/irad6/   CO2(NLAYER), RDH2O(NLAYER),   O2(NLAYER),         
     1                O3(NLAYER), CAER(NRAD,NLAYER,NGROUP), 
     2                PRESS(NLAYER), PBAR(NLAYER),                                     
     3                DPG(NLAYER), TT(NLAYER), Y3(NTOTAL,NGAUSS,NLAYER),
     4                TGRND,  U0,  FDEGDAY, ISL, IR, IRS
C                                                                       
C DEFINED IN 'OPPROP'                                                   
C
      COMMON /irad7/                                                   
     1  WOT, GOT,                                                       
     2  PTEMPG(NTOTAL),        PTEMPT(NTOTAL),
     3  G0(   NTOTAL,NLAYER),  OPD( NTOTAL,NLAYER),                     
     3  PTEMP(NTOTAL,NLAYER),  TAUL(NTOTAL,NLAYER),
     5  TAUH2O(NTOTAL,NLAYER), TAUS(NWAVE,NLAYER),
     6  TAUA(NWAVE,NLAYER),    G01(NWAVE,NLAYER),
     7  uG0(   NTOTAL,NLAYER), uTAUL(NTOTAL,NLAYER),
     8  W0(   NTOTAL,NLAYER),                                           
     9  uW0(   NTOTAL,NLAYER) , uopd(NTOTAL,NLAYER)
C
C DEFINED IN 'TWOSTR'                                                   
C
      COMMON /irad8/                                                  
     1  U1S( NTOTAL),           U1I( NTOTAL),                           
     2  ACON(NTOTAL,NLAYER),   AK(  NTOTAL,NLAYER),                     
     3  BCON(NTOTAL,NLAYER),   B1(  NTOTAL,NLAYER),                     
     4  B2(  NTOTAL,NLAYER),   EE1( NTOTAL,NLAYER),                     
     5  EM1(NTOTAL,NLAYER),                                             
     6  EM2(NTOTAL,NLAYER),    EL1( NTOTAL,NLAYER),                     
     7  EL2(NTOTAL,NLAYER),    GAMI(NTOTAL,NLAYER),                     
     8  AF(NTOTAL,NDBL), BF(NTOTAL,NDBL), EF(NTOTAL,NDBL)               
C
C DEFINED IN 'ADD'                                                      
C
      COMMON /irad9/                                                   
     1  SFCS(NTOTAL),                                                   
     2  B3(  NTOTAL,NLAYER),   CK1(   NTOTAL,NLAYER),                   
     3  CK2( NTOTAL,NLAYER),   CP(    NTOTAL,NLAYER),                   
     4  CPB( NTOTAL,NLAYER),   CM(    NTOTAL,NLAYER),                   
     5  CMB( NTOTAL,NLAYER),   DIRECT(NTOTAL,NLAYER),                   
     6  EE3( NTOTAL,NLAYER),   EL3(   NTOTAL,NLAYER),                   
     7  FNET(NTOTAL,NLAYER),   TMI(   NTOTAL,NLAYER),                   
     8  AS(  NTOTAL,NDBL),     DF(    NTOTAL,NDBL),                     
     9  DS(  NTOTAL,NDBL),     XK(    NTOTAL,NDBL)                      
C
C DEFINED IN 'NEWFLUX1'                                                 
C
      COMMON /irad10/                                                   
     1  WEIT(   NTOTAL),                                                
     2  DIREC(  NTOTAL,NLAYER), DIRECTU(NTOTAL,NLAYER),                 
     3  SLOPE(  NTOTAL,NLAYER),                                         
     4  DINTENT(NTOTAL,NGAUSS,NLAYER),                                  
     5  UINTENT(NTOTAL,NGAUSS,NLAYER),                                  
     6  TMID(NTOTAL,NLAYER),TMIU(NTOTAL,NLAYER)
c
c printed in 'radout' (defined in 'radtran')
c
      common /irad11/
     1  tslu,tsld,alb_tot,tiru,firu(NIR),firn(NIR),fsLu(NSOL),
     2  fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL),
     3  fupbs(NLAYER),fdownbs(NLAYER),fnetbs(NLAYER),
     4  fupbi(NLAYER),fdownbi(NLAYER),fnetbi(NLAYER),
     5  qrad(NLAYER),alb_tomi,alb_toai
c
c ensure all rad local variables are stored statically
c
      save
