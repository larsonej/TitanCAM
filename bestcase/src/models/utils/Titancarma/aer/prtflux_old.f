      subroutine prtflux(open_file,close_file,prt_cloud)
c
c
c  @(#) prtflux.f  Barth  May-2001
c  This routine outputs diffusion and advection fluxes at the current 
c  timestep to an output print file.
c
c  Argument list input:
c    open_file
c    close_file
c    prt_cloud
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
    1 format(2(1pe13.6))
    2 format(1pe13.6)
c
c
       character*(50) fluxfile
       logical open_file,close_file,prt_cloud
c
c
c---------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter prtflux'
      return
c
c
c   ..Write out fluxes to files:
c 
       if(prt_cloud) then

         if(open_file) then

c   ..Initial file names:
c
c   ..First file name tag to be appended to initial file name
           ieins=ieins+1
           if (ieins.gt.9) then
             ieins=0
c   ..Second file name tag
             izwei=izwei+1
             if (izwei.gt.9) then
               izwei=0
               idrei=idrei+1
               if (idrei.gt.9) then
                 idrei=0
                 ivier=ivier+1
               endif
             endif
           endif
c
c  ..Construct file names..
c
           fluxfile='Files/Output/mflux'//ext//char(ivier+48)//
     $       char(idrei+48)//char(izwei+48)//char(ieins+48)//'.p'
c
c   ..Report file name writen to...
c
           print*,"Writing to:"
           print*,fluxfile

           open(unit=LUNOMFLX,file=fluxfile,STATUS='UNKNOWN')

           open_file = .false.

c        ...Time for matching to particle/gas concentrations
           write(LUNOMFLX,*),time+dtime

c        ...Diffusion parameters are independent of particle size
           do k=1,NZ
             write(LUNOMFLX,1) vertdifu(k),vertdifd(k)
           enddo

         endif
 
         do k = 1,NZ
           write(LUNOMFLX,2) vertadvd(k)
c          write(LUNOMFLX,2) fluxcl3(k)
         enddo

       else

c          write(LUNOMFLX,1) 0.00
c        do k = 2,NZ
         do k = 1,NZ
c          cvert_bnd = ( cvert(k-1) + cvert(k) ) / 2.d0
c          write(LUNOMFLX,1) vertdifu(k)*cvert_bnd,vertdifd(k)*cvert_bnd
           write(LUNOMFLX,1) vertdifu(k),vertdifd(k)
         enddo

         if(close_file) CLOSE(unit=LUNOMFLX)
       endif
c
c
c  Return to caller with flux info output to print file
c
      return
      end
