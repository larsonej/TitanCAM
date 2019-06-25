      subroutine quit
c
c
c  @(#) quit.f  McKie  Oct-1995
c  This routine terminates a normal model run.
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
    1 format(a,' terminating normally')
    2 format('Last timestep index: ',i6)
    3 format('Last time: ',f12.2)
    4 format(a,' output was to file ',a)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter quit'
c
c
c  Report last timestep & time
c
      call prtsep
      write(LUNOPRT,1)
      write(LUNOPRT,2) itime
      write(LUNOPRT,3) time
c
c
c  Close files for time step diagnostic output
c
c     close(unit=LUNOSTEP)
c
c
c  Report name of output history file
c
c      call dblank(hisofil, ns)
c      write(LUNOPRT,4) 'History', hisofil(1:ns)
c
c
c  Report name of output restart file
c
c      call dblank(resofil, ns)
c      write(LUNOPRT,4) 'Restart', resofil(1:ns)
c
c
c  Close output restart file
c
c      close(unit=LUNORES)
c
c
c  Close output history file
c
c      if( do_netcdf )then
c       call ncsnc(ncdf_file, ierr)
c       call ncclos(ncdf_file, ierr)
c      else
c       close(unit=LUNOHIS)
c      endif
c
c
c  Close output print file
c
c      call prtsep
c      close(unit=LUNOPRT)
c
c
c  Terminate the program's execution
c
      stop
      end
