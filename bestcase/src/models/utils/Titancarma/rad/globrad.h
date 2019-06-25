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
C IVERT  = MAXIMUM NUMBER OF LAYERS;                                    
C ILAYER = MAXIMUM NUMBER OF LAYER BOUNDARIES                           
C IDBL   = TWICE THE MAXIMUM NUMBER OF LAYER BOUNDARIES                 
C IRAD   = MAXIMUM NUMBER OF AEROSOL RADIUS BINS;                       
C NGROUP = number of aerosol types
c IACAP  = size of the complex work array used (by the Mie scattering code)
c          ACAP (which is not described)
C                                                                       
      PARAMETER ( IVERT = NZ_RAD )                                              
      PARAMETER ( IRAD = NBIN_AERAD )                                              
      PARAMETER ( NGROUP = NGROUP_AERAD )                                          
      parameter ( IACAP = 7000 )
      PARAMETER ( ILAYER = IVERT+1, IDBL = 2*ILAYER )                         
      PARAMETER ( IRADVER = IRAD*IVERT )                                    
      PARAMETER ( IRADLAY = IRAD*ILAYER )                                   
C                                                                       
C DEFINE THE DIMENSIONS USED ONLY BY THE RADIATION MODEL                
C                                                                       
C ISOL   = TOTAL NUMBER OF SOLAR WAVELENGTHS                            
C IIR    = TOTAL NUMBEROF INFRARED WAVELENGTHS                          
C IWAVE  = ISOL + IIR                                                   
C                                                                       
      PARAMETER ( ISOL = NSOL_AERAD )
      PARAMETER ( IIR = NIR_AERAD )
      PARAMETER ( IWAVE = NWAVE_AERAD )
C                                                                       
C ITOTAL = TOTAL NUMBER OF PROBABILITY INTERVALS;                       
C ISOLP  = NUMBER OF SOLAR PROBABILITY INTERVALS;                       
C IIRP   = NUMBER OF INFRARED PROBABILITY INTERVALS;                    
C                                                                       
      PARAMETER ( ISOLP = 77, IIRP = 71, ITOTAL = ISOLP + IIRP)         
C                                                                       
C IGAUSS = TOTAL NUMBER OF GAUSS QUADRATURE POINTS;                     
C                                                                       
      PARAMETER ( IGAUSS = 3)                                           
C                                                                       
C ICOUNT = USED TO CALCULATE PLANK FUNCTION.                            
C                                                                       
      PARAMETER ( ICOUNT = 20000 )
C                                                                       
C LUNOPRT = UNIT NUMBER FOR PRINT OUTPUT                                   
C                                                                       
       PARAMETER ( LUNOPRT = LUNOPRT_AERAD )                                              
C                                                                       
C LUNMIE = unit number for input and output of Mie coefficients
C                                                                       
       parameter ( LUNMIE = LUNMIE_AERAD )                                              
c
c
c  Define logical unit number for radiation submodel print output
c
      parameter( LUNORAD = LUNORAD_AERAD )
c
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
     3 O3MIX(ILAYER), O3MIXP(6), O3C, VRAT,                         
     2 PTOP, PBOT, RMIN(NGROUP), R(IRAD,NGROUP),
     3 is_grp_ice(NGROUP), NVERT, NLAYER, NRAD
C                                                                       
C TIME-DEPENDENT VARIABLES THAT MIGHT BE SPECIFIED BY AN EXTERNAL MODEL 
C                                                                       
      COMMON /irad2/                                                  
     1 U0EXT,                                                       
     2 ALBEDO_SFC, EMISIR, 
     3 P(IVERT), T(IVERT), Q(IVERT)
C                                                                       
C OUTPUT VARIABLES, CALCULATED BY THE RADIATION MODEL, THAT MIGHT       
C BE USED BY THE EXTERNAL MODEL                                         
C                                                                       
      COMMON /irad3/                                                  
     1 HEATI(ILAYER), HEATS(ILAYER), HEAT(ILAYER),                      
     2 SOLNET, XIRDOWN, XIRUP                                           
C                                                                       
C INITIATED IN SETUPRAD FOR RADIATION CALCULATION                       
C
      COMMON /irad4/                                                  
     1  LLA, LLS, JDBLE, JN, EPSILON, EXPMAX,                           
     2  TPI, SQ3, SBK,                                                  
     2  AM, AVG, ALOS, G, PI, SCDAY, RGAS,                              
     3  NSOL, NIR, NSOLP, NIRP, NTOTAL, NWAVE,                          
     3  GANGLE(IGAUSS), GWEIGHT(IGAUSS), SFLX(ISOL), WVLN(ISOL),        
     4  GRATIO(IGAUSS),
     2  EMIS(ITOTAL), RSFX(ITOTAL),LTEMP(ITOTAL),NPROB(ITOTAL),         
     3  SOL(ITOTAL),TAURAY(ITOTAL),WEIGHT(ITOTAL),                      
     4  GCLD(  ITOTAL,ILAYER),   GOL(ITOTAL,ILAYER),                    
     5  PARAY( ITOTAL,ILAYER),TAUAER(ITOTAL,ILAYER),                    
     6  WCLD(ITOTAL,ILAYER),                                            
     6  TAUCLD(ITOTAL,ILAYER),WOL(   ITOTAL,ILAYER),                    
     7  TREAL(2,IWAVE), TTMAG(2,IWAVE),
     8  contnm(iirp), nprobi(iwave,2)
C                                                                       
      COMMON/irad5/                                                  
     1       ACO2(ITOTAL),           AH2O(ITOTAL),                      
     2       AO2(ITOTAL),            AO3(ITOTAL),                       
     6       PACO2(ITOTAL,ILAYER),   PAH2O(ITOTAL,ILAYER),              
     7       PAO2(ITOTAL,ILAYER),    PAO3(ITOTAL,ILAYER),               
     8       PLANK(IIR+1,ICOUNT),                                       
     9       PSCO2(ITOTAL),          PSH2O(ITOTAL),                     
     1       PSO2(ITOTAL),           PSO3(ITOTAL),                      
     2       SOLFX(ISOL),            WAVE(IWAVE+1),                     
     3       TAUGAS(ITOTAL,ILAYER),  NLOW, NHIGH,                       
     4       XSECTA(IRAD,NGROUP),    RUP(IRAD,NGROUP),                         
     5       QSCAT(IRAD,NGROUP,IWAVE), blackbody_above,
     6       QBRQS(IRAD,NGROUP,IWAVE), t_above,     
     7       RDQEXT(IRAD,NGROUP,IWAVE)                                         
C                                                                       
      COMMON/irad6/   CO2(ILAYER), RDH2O(ILAYER),   O2(ILAYER),         
     1                O3(ILAYER), CAER(IRAD,ILAYER,NGROUP), 
     2                PRESS(ILAYER), PBAR(ILAYER),                                     
     3                DPG(ILAYER), TT(ILAYER), Y3(ITOTAL,IGAUSS,ILAYER),
     4                TGRND,  U0,  ISL, IR, IRS, NGAUSS, FDEGDAY       
C                                                                       
C DEFINED IN 'OPPROP'                                                   
C
      COMMON /irad7/                                                   
     1  WOT, GOT,                                                       
     2  PTEMPG(ITOTAL),        PTEMPT(ITOTAL),
     3  G0(   ITOTAL,ILAYER),  OPD( ITOTAL,ILAYER),                     
     3  PTEMP(ITOTAL,ILAYER),  TAUL(ITOTAL,ILAYER),
     5  TAUH2O(ITOTAL,ILAYER), TAUS(IWAVE,ILAYER),
     6  TAUA(IWAVE,ILAYER),    G01(IWAVE,ILAYER),
     7  uG0(   ITOTAL,ilayer), uTAUL(ITOTAL,ilayer),
     8  W0(   ITOTAL,ILAYER),                                           
     9  uW0(   ITOTAL,ilayer) , uopd(itotal,ilayer)
C
C DEFINED IN 'TWOSTR'                                                   
C
      COMMON /irad8/                                                  
     1  U1S( ITOTAL),           U1I( ITOTAL),                           
     2  ACON(ITOTAL,ILAYER),   AK(  ITOTAL,ILAYER),                     
     3  BCON(ITOTAL,ILAYER),   B1(  ITOTAL,ILAYER),                     
     4  B2(  ITOTAL,ILAYER),   EE1( ITOTAL,ILAYER),                     
     5  EM1(ITOTAL,ILAYER),                                             
     6  EM2(ITOTAL,ILAYER),    EL1( ITOTAL,ILAYER),                     
     7  EL2(ITOTAL,ILAYER),    GAMI(ITOTAL,ILAYER),                     
     8  AF(ITOTAL,IDBL), BF(ITOTAL,IDBL), EF(ITOTAL,IDBL)               
C
C DEFINED IN 'ADD'                                                      
C
      COMMON /irad9/                                                   
     1  SFCS(ITOTAL),                                                   
     2  B3(  ITOTAL,ILAYER),   CK1(   ITOTAL,ILAYER),                   
     3  CK2( ITOTAL,ILAYER),   CP(    ITOTAL,ILAYER),                   
     4  CPB( ITOTAL,ILAYER),   CM(    ITOTAL,ILAYER),                   
     5  CMB( ITOTAL,ILAYER),   DIRECT(ITOTAL,ILAYER),                   
     6  EE3( ITOTAL,ILAYER),   EL3(   ITOTAL,ILAYER),                   
     7  FNET(ITOTAL,ILAYER),   TMI(   ITOTAL,ILAYER),                   
     8  AS(  ITOTAL,IDBL),     DF(    ITOTAL,IDBL),                     
     9  DS(  ITOTAL,IDBL),     XK(    ITOTAL,IDBL)                      
C
C DEFINED IN 'NEWFLUX1'                                                 
C
      COMMON /irad10/                                                   
     1  WEIT(   ITOTAL),                                                
     2  DIREC(  ITOTAL,ILAYER), DIRECTU(ITOTAL,ILAYER),                 
     3  SLOPE(  ITOTAL,ILAYER),                                         
     4  DINTENT(ITOTAL,IGAUSS,ILAYER),                                  
     5  UINTENT(ITOTAL,IGAUSS,ILAYER),                                  
     6  TMID(ITOTAL,ILAYER),TMIU(ITOTAL,ILAYER)
c
c printed in 'radout' (defined in 'radtran')
c
      common /irad11/
     1  tslu,tsld,alb_tot,tiru,firu(iir),firn(iir),fsLu(isoL),
     2  fsLd(isoL),fsLn(isoL),alb_toa(isol),
     3  fupbs(ilayer),fdownbs(ilayer),fnetbs(ilayer),
     4  fupbi(ilayer),fdownbi(ilayer),fnetbi(ilayer),
     5  qrad(ilayer),alb_tomi,alb_toai
c
c ensure all rad local variables are stored statically
c
      save
