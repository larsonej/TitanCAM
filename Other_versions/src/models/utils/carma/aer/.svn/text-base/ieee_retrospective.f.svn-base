      subroutine ieee_retrospective
c
c
c  @(#) ieee_retrospective.f  McKie Feb-1991
c  This routine takes the place of the system ieee floating point
c  uncleared exception reporting routine.  It reports only the
c  serious (so-called "common") ieee floating point exceptions
c  which have not been cleared.  This routine can be called explicitly
c  at any time by the user program, or is called implicitly at the
c  end of a normal fortran program by the fortran run-time system. 
c
c  This version is for Sun SunOS 4.1, SunFORTRAN v1.4
c
c  Define sun4 ieee fl pt exception bits  (See /usr/include/sys/ieeefp.h)
c
      parameter(IINEXACT=1)
      parameter(IDIVISION=2)
      parameter(IUNDERFLOW=4)
      parameter(IOVERFLOW=8)
      parameter(IINVALID=16)
c
c
c  Declare local storage
c
      character*(16) out
c
c
c  Define formats
c
    1 format(a)
    2 format(a,$)
c
c
c  Get bit pattern showing which ieee fl pt exceptions are uncleared
c
      ibits = ieee_flags('get','exception','all',out)
c
c
c  Define bit mask for serious ieee fl pt exceptions
c
      mask = IDIVISION + IOVERFLOW + IINVALID
c
c
c  List serious exceptions which are uncleared
c
      if( and(ibits,mask) .ne. 0 )then
       write(*,1) ' '
       write(*,1) 'Uncleared serious ieee fl pt exceptions:'
       if( and(ibits,IDIVISION) .ne. 0 ) write(*,2) ' Divide by zero; '
       if( and(ibits,IOVERFLOW) .ne. 0 ) write(*,2) ' Overflow; '
       if( and(ibits,IINEXACT) .ne. 0 ) write(*,2) ' Inexact; '
       write(*,1) ' '
      endif
c
c
c  Clear all exception flags
c
      ibits = ieee_flags('clear','exception','all',out)
c
c
c  Return to caller with retrospective complete 
c
      return
      end
