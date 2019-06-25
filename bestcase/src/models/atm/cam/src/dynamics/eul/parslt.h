!
! $Id: parslt.h 17 2006-12-11 21:50:24Z hpc $
! $Author: hpc $
!
!
! Parameters common to many SLT routines
!
      integer ppdy              ! length of interpolation grid stencil
      logical plimdr            ! flag to limit derivatives
!
      parameter(ppdy = 4,  plimdr = .true.)
!
 
