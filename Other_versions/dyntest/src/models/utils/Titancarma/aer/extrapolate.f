      subroutine extrapolate
c
c
c  @(#) extrapolate.f  Barth  Jan-2001
c  This routine calculates new gas concentrations by extrapolating from 
c  previous runs over some timestep
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
c-----------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter extrapolate'
c
c-----------------------------------------------------------------------------
c
c
      extime = tlex + trecov 

      dt_small = .001 * dtime

      if( (time+dt_small) .ge. extime ) then

        do igas=1,NGAS
         do k=1,NZ
          gc3(k,igas) = 2.*gc3(k,igas) - gcb(k,igas)
          gcb(k,igas) = gc3(k,igas)
         enddo
        enddo

        tlex = time
        time = 2.*time - tbegin
        tbegin = time

      else

        do igas=1,NGAS
         do k=1,NZ
          gcb(k,igas) = gc3(k,igas)
         enddo
        enddo
        tbegin = time
        prt_year = .true.

      endif
c
c
c  Return to caller with gc updated 
c
      return
      end
