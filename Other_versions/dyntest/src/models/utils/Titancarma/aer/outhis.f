      subroutine outhis
c
c
c  @(#) outhis.f  McKie  Jan-1997
c  This routine controls history file output.
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
    1 format('History output # ',i6,
     $       ' at itime: ',i6,3x,'time: ',f12.2)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter outhis'
c
c
c  Output current timepoint's history state to netcdf or binary file as requested
c
      if( do_netcdf )then
       call outhis_ncdf
      else
       call outhis_bin
      endif
c
c
c  Increment the counter of history timepoints output
c
      khist = khist + 1
c
c
c  Report history output
c
cc    call prtsep
cc    write(LUNOPRT,1) khist, itime, time
c
c
c  Return to caller with history output complete
c
      return
      end
