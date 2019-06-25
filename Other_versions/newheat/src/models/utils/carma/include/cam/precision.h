!  @(#) precision.h  McKie  Oct-1995
!
!  This file should be included in the global include file
!  as well as any source code routines that do not include
!  the global include file.
!
!  It declares implicit types for integer & floating pt variables.
!  It also defines precision-dependent global parameters.
!
!   Note:
!    Logical & character variables should be explicitly
!    declared, but floating points and integers should not be
!    explicitly declared so that the following implicit types
!    will automatically be assigned to floats and ints.
!
!
!  This file is organized into 3 sections -- A, B, & C.
!  It is intended that the statements in one of these sections
!  be uncommented to select a particular precision, and the
!  other two sections are then commented.
!
!  The two possible precisions supported here are 32 bit & 64 bit.
!
!
!  Typical recommended precision for various systems:
!
!   For 32 bit precision:
!
!     sunos:  uncomment Section A, comment Sections B & C
!     sgi:    uncomment Section A, comment Sections B & C
!     linux:  uncomment Section A, comment Sections B & C
!     hp:     uncomment Section A, comment Sections B & C
!     unicos: (32 bit precision not available on unicos)
!
!   For 64 bit precision:
!
!     sunos:  uncomment Section B, comment Sections A & C
!     sgi:    uncomment Section B, comment Sections A & C
!     linux:  uncomment Section B, comment Sections A & C
!     hp:     uncomment Section B, comment Sections A & C
!     unicos: uncomment Section C, comment Sections A & B
!
!
!=============================================================================
!
!  Section A:  32 bit word length machine, using 32 bit single precision
!
!
!  Define implicit type for floating point variables
!
!     implicit real ( a-h, o-z )
!
!
!  Define implicit type for integer variables
!
!     implicit integer ( i-n )
!
!
!  Define double precision 1.0 to be multiplied by literal constants
!  that are embedded within parentheses.
!
!     parameter( ONE = 1.d0 )
!
!
!  Define smallest possible number such that ONE + ALMOST_ZERO > ONE
!
!     parameter( ALMOST_ZERO = 1.d-7 )
!
!
!  Define value slightly < 1.0 to test for total evaporation [ dimensionless ]
!
!     parameter( ALMOST_ONE = ONE - ALMOST_ZERO )
!
!
!  Define value for maximum exponential argument [ dimensionless ]
!
!     parameter( POWMAX = 85.d0 )
!
!
!  Define small particle number concentration [ # / x_units / y_units / z_units ]
!
!     parameter( SMALL_PC = 1.d-25 )
!
!
!  Define (arbitrary) symbols to specify precision of netcdf integers
!
!     parameter( NCDF_LONG = 1 )
!     parameter( NCDF_SHORT = 2 )
!
!
!  Define (arbitrary) symbols to specify precision of netcdf floating points
!
!     parameter( NCDF_DOUBLE = 1 )
!     parameter( NCDF_FLOAT = 2 )
!
!
!  Define symbols to control the netcdf types for integers & floating points
!  that are selected at run time in the outhis_ncdf.f routine when <do_ncdf>
!  is .true.
!
!     parameter( NCDF_INT = NCDF_LONG )
!     parameter( NCDF_FP  = NCDF_FLOAT )
!
!
!=============================================================================
!
!  Section B:  32 bit word length machine, using 64 bit double precision
!
!
!  Define implicit type for floating point variables
!
      implicit real ( a-h, o-z )
!
!
!  Define implicit type for integer variables
!
      implicit integer ( i-n )
!
!
!  Define double precision 1.0 to be multiplied by literal constants
!  that are embedded within parentheses.
!
      parameter( ONE = 1.d0 )
!
!
!  Define smallest possible number such that ONE + ALMOST_ZERO > ONE
!
      parameter( ALMOST_ZERO = 1.d-15 )
!
!
!  Define value slightly < 1.0 to test for total evaporation [ dimensionless ]
!
      parameter( ALMOST_ONE = ONE - ALMOST_ZERO )
!
!
!  Define value for maximum exponential argument [ dimensionless ]
!
      parameter( POWMAX = 700.d0 )
!
!
!  Define small particle number concentration [ # / x_units / y_units / z_units ]
!
      parameter( SMALL_PC = 1.d-50 )
!
!
!  Define (arbitrary) symbols to specify precision of netcdf integers
!
      parameter( NCDF_LONG = 1 )
      parameter( NCDF_SHORT = 2 )
!
!
!  Define (arbitrary) symbols to specify precision of netcdf floating points
!
      parameter( NCDF_DOUBLE = 1 )
      parameter( NCDF_FLOAT = 2 )
!
!
!  Define symbols to control the netcdf types for integers & floating points
!  that are selected at run time in the outhis_ncdf.f routine when <do_ncdf>
!  is .true.
!
      parameter( NCDF_INT = NCDF_LONG )
      parameter( NCDF_FP  = NCDF_DOUBLE )
!
!
!=============================================================================
!
!  Section C:  64 bit word length machine, using 64 bit single precision
!
!
!--inactive--c  Define implicit type for floating point variables
!--inactive--c
!--inactive--      implicit real ( a-h, o-z )
!--inactive--c
!--inactive--c
!--inactive--c  Define implicit type for integer variables
!--inactive--c
!--inactive--      implicit integer ( i-n )
!--inactive--c
!--inactive--c
!--inactive--c  Define double precision 1.0 to be multiplied by literal constants
!--inactive--c  that are embedded within parentheses.
!--inactive--c
!--inactive--      parameter( ONE = 1.d0 )
!--inactive--c
!--inactive--c
!--inactive--c  Define smallest possible number such that ONE + ALMOST_ZERO > ONE
!--inactive--c
!--inactive--      parameter( ALMOST_ZERO = 1.d-15 )
!--inactive--c
!--inactive--c
!--inactive--c  Define value slightly < 1.0 to test for total evaporation [ dimensionless ]
!--inactive--c
!--inactive--      parameter( ALMOST_ONE = ONE - ALMOST_ZERO )
!--inactive--c
!--inactive--c
!--inactive--c  Define value for maximum exponential argument [ dimensionless ]
!--inactive--c
!--inactive--      parameter( POWMAX = 5500.d0 )
!--inactive--c
!--inactive--c
!--inactive--c  Define small particle number concentration [ # / x_units / y_units / z_units ]
!--inactive--c
!--inactive--      parameter( SMALL_PC = 1.d-50 )
!--inactive--c
!--inactive--c
!--inactive--c  Define (arbitrary) symbols to specify precision of netcdf integers
!--inactive--c
!--inactive--      parameter( NCDF_LONG = 1 )
!--inactive--      parameter( NCDF_SHORT = 2 )
!--inactive--c
!--inactive--c
!--inactive--c  Define (arbitrary) symbols to specify precision of netcdf floating points
!--inactive--c
!--inactive--      parameter( NCDF_DOUBLE = 1 )
!--inactive--      parameter( NCDF_FLOAT = 2 )
!--inactive--c
!--inactive--c
!--inactive--c  Define symbols to control the netcdf types for integers & floating points
!--inactive--c  that are selected at run time in the outhis_ncdf.f routine when <do_ncdf>
!--inactive--c  is .true.
!--inactive--c
!--inactive--      parameter( NCDF_INT = NCDF_LONG )
!--inactive--      parameter( NCDF_FP  = NCDF_DOUBLE )