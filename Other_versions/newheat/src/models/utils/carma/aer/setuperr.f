      subroutine setuperr
c
c
c  @(#) setuperr.f  McKie  Nov-1999
c
c
c  This routine is called to handle any run-time set up of error exception
c  handling.  Normally, it is a do-nothing routine.  But system dependent
c  error handling code could be inserted here, e.g. temporarily for debugging.
c
c
      return
      end
