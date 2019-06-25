c  @(#) precision.h  McKie  Oct-1995
c
c  This file should be included in the global include file
c  as well as any source code routines that do not include
c  the global include file.
c
c  It declares implicit types for integer & floating pt variables.
c  It also defines precision-dependent global parameters.
c
c   Note:
c    Logical & character variables should be explicitly
c    declared, but floating points and integers should not be
c    explicitly declared so that the following implicit types
c    will automatically be assigned to floats and ints.
c
c
c  This file is organized into 3 sections -- A, B, & C.
c  It is intended that the statements in one of these sections
c  be uncommented to select a particular precision, and the
c  other two sections are then commented.
c
c  The two possible precisions supported here are 32 bit & 64 bit.
c
c
c  Typical recommended precision for various systems:
c
c   For 32 bit precision:
c
c     sunos:  uncomment Section A, comment Sections B & C
c     sgi:    uncomment Section A, comment Sections B & C
c     linux:  uncomment Section A, comment Sections B & C
c     hp:     uncomment Section A, comment Sections B & C
c     unicos: (32 bit precision not available on unicos)
c
c   For 64 bit precision:
c
c     sunos:  uncomment Section B, comment Sections A & C
c     sgi:    uncomment Section B, comment Sections A & C
c     linux:  uncomment Section B, comment Sections A & C
c     hp:     uncomment Section B, comment Sections A & C
c     unicos: uncomment Section C, comment Sections A & B
c
c
c=============================================================================
c
c  Section A:  32 bit word length machine, using 32 bit single precision
c
c
c  Define implicit type for floating point variables
c
      implicit real ( a-h, o-z )
c
c
c  Define implicit type for integer variables
c
      implicit integer ( i-n )
c
c
c  Define double precision 1.0 to be multiplied by literal constants
c  that are embedded within parentheses.
c
      parameter( ONE = 1.d0 )
c
c
c  Define smallest possible number such that ONE + ALMOST_ZERO > ONE
c
      parameter( ALMOST_ZERO = 1.d-7 )
c
c
c  Define value slightly < 1.0 to test for total evaporation [ dimensionless ]
c
      parameter( ALMOST_ONE = ONE - ALMOST_ZERO )
c
c
c  Define value for maximum exponential argument [ dimensionless ]
c
      parameter( POWMAX = 85.d0 )
c
c
c  Define tolerance for series convergence in Mie calculations:
c  should be 1.d-7 for 32 bit machines and 1.d-14 for 64 bit machines,
c  regardless of implicit precision defined above
c
      parameter( EPSILON_MIE = 1.d-7 )
c
c
c  Define small particle number concentration [ # / x_units / y_units / z_units ]
c
      parameter( SMALL_PC = 1.d-25 )
c
c
c  Define (arbitrary) symbols to specify precision of netcdf integers
c
      parameter( NCDF_LONG = 1 )
      parameter( NCDF_SHORT = 2 )
c
c
c  Define (arbitrary) symbols to specify precision of netcdf floating points
c
      parameter( NCDF_DOUBLE = 1 )
      parameter( NCDF_FLOAT = 2 )
c
c
c  Define symbols to control the netcdf types for integers & floating points
c  that are selected at run time in the outhis_ncdf.f routine when <do_ncdf>
c  is .true.
c
      parameter( NCDF_INT = NCDF_LONG )
      parameter( NCDF_FP  = NCDF_FLOAT )
c
c
c=============================================================================
c
c  Section B:  32 bit word length machine, using 64 bit double precision
c
c
c--inactive--c  Define implicit type for floating point variables
c--inactive--c
c--inactive--      implicit double precision ( a-h, o-z )
c--inactive--c
c--inactive--c
c--inactive--c  Define implicit type for integer variables
c--inactive--c
c--inactive--      implicit integer ( i-n )
c--inactive--c
c--inactive--c
c--inactive--c  Define double precision 1.0 to be multiplied by literal constants
c--inactive--c  that are embedded within parentheses.
c--inactive--c
c--inactive--      parameter( ONE = 1.d0 )
c--inactive--c
c--inactive--c
c--inactive--c  Define smallest possible number such that ONE + ALMOST_ZERO > ONE
c--inactive--c
c--inactive--      parameter( ALMOST_ZERO = 1.d-15 )
c--inactive--c
c--inactive--c
c--inactive--c  Define value slightly < 1.0 to test for total evaporation [ dimensionless ]
c--inactive--c
c--inactive--      parameter( ALMOST_ONE = ONE - ALMOST_ZERO )
c--inactive--c
c--inactive--c
c--inactive--c  Define value for maximum exponential argument [ dimensionless ]
c--inactive--c
c--inactive--      parameter( POWMAX = 700.d0 )
c--inactive--c
c--inactive--c
c--inactive--c  Define tolerance for series convergence in Mie calculations:
c--inactive--c  should be 1.d-7 for 32 bit and 1.d-14 for 64 bit machines,
c--inactive--c  regardless of implicit precision defined above
c--inactive--c
c--inactive--      parameter( EPSILON_MIE = 1.d-7 )
c--inactive--c
c--inactive--c
c--inactive--c  Define small particle number concentration [ # / x_units / y_units / z_units ]
c--inactive--c
c--inactive--      parameter( SMALL_PC = 1.d-50 )
c--inactive--c
c--inactive--c
c--inactive--c  Define (arbitrary) symbols to specify precision of netcdf integers
c--inactive--c
c--inactive--      parameter( NCDF_LONG = 1 )
c--inactive--      parameter( NCDF_SHORT = 2 )
c--inactive--c
c--inactive--c
c--inactive--c  Define (arbitrary) symbols to specify precision of netcdf floating points
c--inactive--c
c--inactive--      parameter( NCDF_DOUBLE = 1 )
c--inactive--      parameter( NCDF_FLOAT = 2 )
c--inactive--c
c--inactive--c
c--inactive--c  Define symbols to control the netcdf types for integers & floating points
c--inactive--c  that are selected at run time in the outhis_ncdf.f routine when <do_ncdf>
c--inactive--c  is .true.
c--inactive--c
c--inactive--      parameter( NCDF_INT = NCDF_LONG )
c--inactive--      parameter( NCDF_FP  = NCDF_DOUBLE )
c
c
c=============================================================================
c
c  Section C:  64 bit word length machine, using 64 bit single precision
c
c
c--inactive--c  Define implicit type for floating point variables
c--inactive--c
c--inactive--      implicit real ( a-h, o-z )
c--inactive--c
c--inactive--c
c--inactive--c  Define implicit type for integer variables
c--inactive--c
c--inactive--      implicit integer ( i-n )
c--inactive--c
c--inactive--c
c--inactive--c  Define double precision 1.0 to be multiplied by literal constants
c--inactive--c  that are embedded within parentheses.
c--inactive--c
c--inactive--      parameter( ONE = 1.d0 )
c--inactive--c
c--inactive--c
c--inactive--c  Define smallest possible number such that ONE + ALMOST_ZERO > ONE
c--inactive--c
c--inactive--      parameter( ALMOST_ZERO = 1.d-15 )
c--inactive--c
c--inactive--c
c--inactive--c  Define value slightly < 1.0 to test for total evaporation [ dimensionless ]
c--inactive--c
c--inactive--      parameter( ALMOST_ONE = ONE - ALMOST_ZERO )
c--inactive--c
c--inactive--c
c--inactive--c  Define value for maximum exponential argument [ dimensionless ]
c--inactive--c
c--inactive--      parameter( POWMAX = 5500.d0 )
c--inactive--c
c--inactive--c
c--inactive--c  Define tolerance for series convergence in Mie calculations:
c--inactive--c  should be 1.d-7 for 32 bit machines and 1.d-14 for 64 bit machines,
c--inactive--c  regardless of implicit precision defined above
c--inactive--c
c--inactive--      parameter( EPSILON_MIE = 1.d-14 )
c--inactive--c
c--inactive--c
c--inactive--c  Define small particle number concentration [ # / x_units / y_units / z_units ]
c--inactive--c
c--inactive--      parameter( SMALL_PC = 1.d-50 )
c--inactive--c
c--inactive--c
c--inactive--c  Define (arbitrary) symbols to specify precision of netcdf integers
c--inactive--c
c--inactive--      parameter( NCDF_LONG = 1 )
c--inactive--      parameter( NCDF_SHORT = 2 )
c--inactive--c
c--inactive--c
c--inactive--c  Define (arbitrary) symbols to specify precision of netcdf floating points
c--inactive--c
c--inactive--      parameter( NCDF_DOUBLE = 1 )
c--inactive--      parameter( NCDF_FLOAT = 2 )
c--inactive--c
c--inactive--c
c--inactive--c  Define symbols to control the netcdf types for integers & floating points
c--inactive--c  that are selected at run time in the outhis_ncdf.f routine when <do_ncdf>
c--inactive--c  is .true.
c--inactive--c
c--inactive--      parameter( NCDF_INT = NCDF_LONG )
c--inactive--      parameter( NCDF_FP  = NCDF_DOUBLE )
