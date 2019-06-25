      subroutine prtsep
c
c
c  @(#) prtsep.f  McKie  Sep-1997
c  This routine outputs separator line to the print file.
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Define formats
c
    1 format(/,79('='),/)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter prtsep'
c
c
c  Output the separator line
c
      write(LUNOPRT,1)
c
c
c  Return to caller with separator line output to print file
c
      return
      end
