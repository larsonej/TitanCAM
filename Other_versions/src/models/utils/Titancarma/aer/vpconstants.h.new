!  @(#) vpconstants.h  Barth  Jul-2002
!  This is the include file for the constants used in calculating
!  vapor pressures and latent heats
!
!  This file is needed for the <vaporp> and <setupgrow> subroutines  
!
!
!  Define multiplication factors for conversion to cgs 
!
!  from mmHg:
      parameter( RmmHg2cgs = 1.333d3 )
 
!  from atm:
      parameter( Ratm2cgs = 1.013d6 )
!
!  from Pa:
      parameter( RPa2cgs = 10.)
!
!  from Joules:
      parameter( RJ2cgs = 1.d7 )
!
!
!  Define coefficients in saturation vapor pressures
!  (see <vaporp> for references)
!
!   Ethane liquid
      parameter( vplA_C2H6 = 20.6973 )
      parameter( vplB_C2H6 = -1.1341d3 )
      parameter( vplC_C2H6 = -5.2514 )
      parameter( vplD_C2H6 = -9.8774d-11 )
      parameter( vplE_C2H6 = 6.7329d-6 )
!
!   Ethane solid
      parameter( vpiA_C2H6 = 10.01 )
      parameter( vpiB_C2H6 = 1085.0 )
      parameter( vpiC_C2H6 = 0.561 )
!
!   Methane liquid
      parameter( vplA_CH4 = 3.901408 )
      parameter( vplB_CH4 = 437.54809 )
      parameter( vplC_CH4 = 1598.8512 )
      parameter( vplD_CH4 = 154567.02 )
!
!   Methane solid
      parameter( vpiA_CH4 = 4.425070 )
      parameter( vpiB_CH4 = 453.92414 )
      parameter( vpiC_CH4 = 4055.6016 )
      parameter( vpiD_CH4 = 115352.19 )
      parameter( vpiE_CH4 = 1165560.7 )
!
!   Nitrogen liquid
      parameter( vplA_N2 = 1.d5 )
      parameter( vplB_N2 = 3.95 )
      parameter( vplC_N2 = 306. )
!
!
!   Constants from temperature dependent eqns 
!    in Thompson, Zollweg, & Gabis
      parameter( TZGa0 = 0.8096 )
      parameter( TZGa1 = -52.07 )
      parameter( TZGa2 = 5443.0 )
      parameter( TZGb0 = -0.0829 )
      parameter( TZGb1 = 9.34 )
      parameter( TZGc0 = 0.0720 )
      parameter( TZGc1 = -6.27 )
!
!
      dimension phi_N2(50), phi_CH4(50), X_N2(50), X_CH4(50),
     $          Hc_N2(50), Hc_CH4(50)

      dimension Pcrit(2), Ptrip(2), Tcrit(2), Ttrip(2),
     $          a4IS(2), b0IS(2), b1IS(2), psatl(2)
!
!
!    Fugacity coefficients
      data phi_N2/   0.895, 0.907, 0.915, 0.924, 0.932, 0.938, 0.944,
     $ 0.949, 0.954, 0.958, 0.962, 0.965, 0.967, 0.969, 0.971, 0.973,
     $ 0.974, 0.975, 0.976, 0.976, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /
 
      data phi_CH4/  0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.996,
     $ 0.997, 0.997, 0.998, 0.998, 0.998, 0.999, 0.999, 0.999, 0.999,
     $ 0.999, 0.999, 0.999, 0.999, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /
 
!    Saturation mole fractions 
      data X_N2/     0.179, 0.204, 0.214, 0.232, 0.250, 0.264, 0.272,
     $ 0.282, 0.292, 0.300, 0.295, 0.289, 0.286, 0.274, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /
 
      data X_CH4/    0.947, 0.918, 0.905, 0.879, 0.860, 0.840, 0.824, 
     $ 0.807, 0.790, 0.774, 0.775, 0.784, 0.783, 0.795, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /

!    Enthalpy of condensation [J/mol]
      data Hc_N2/ -4952.5, -5068.8, -5146.8, -5224.3, -5296.7, -5359.1, 
     $   -5411.1, -5462.9, -5510.8, -5557.0, -5592.0, -5626.6, -5660.0, 
     $   -5687.2, -5712.3, -5732.0, -5750.9, -5761.3, -5770.3, -5778.4,     
     $       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
     $       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
     $       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
     $       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     $       0.0,     0.0 /

      data Hc_CH4/ -8685.7, -8739.5, -8783.9, -8822.0, -8860.0, -8898.0, 
     $    -8930.6, -8962.5, -8992.4, -9021.5, -9044.9, -9067.0, -9087.8, 
     $    -9105.4, -9121.5, -9133.6, -9145.6, -9152.9, -9158.8, -9164.1,
     $       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
     $       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
     $       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
     $       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     $       0.0,     0.0 /
 
!    Critical and triple point values for N2 and CH4
      data Pcrit/ 34.002, 45.955 /   ! bar
      data Ptrip/ 0.1252, 0.1170 /   ! bar
      data Tcrit/ 126.20, 190.53 /   ! Kelvin
      data Ttrip/ 63.15,   90.68 /   ! Kelvin
 
      data a4IS/ 3.065972, 3.159023 /
      data b0IS/ -23.52451, -19.36816 /
      data b1IS/ 6103.604, 8799.140 /

