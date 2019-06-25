      subroutine step
c
c
c  @(#) step.f  McKie  Oct-1995
c  This routine performs all calculations necessary to
c  take one timestep.
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
    5 format(a,i4,2x,'gc/vol:',f15.10,2x,'gc/cld:',f15.10,2x,
     $       'cmf:',f15.10,2x,'cloud(pc):',1pe13.6,i4)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter step'
c
c
c  Do pre-timestep processing
c
      call prestep
c
c
c  Update model state at new time
c
      call newstate
c
c
c  Modify time-step if necessary
c
      if( do_varstep ) call varstep

c
c----------------------------------------------------------
c  Write out fractions for debugging
c
c     if(do_write(15)) then
c     do k=1,NZ
c      do i=1,NBIN
c       volpart = pc3(k,i,2)*rmass(i,2) -  pc3(k,i,3)
c       if( pc3(k,i,5) .ge. volpart .or.
c    $       pc3(k,i,5) .ge. rmass(i,2)*pc3(k,i,2) ) then
c        write(*,*) 'gcmf > 1, after varstep for alt',k
c        iwa = k
c       endif
c      enddo
c     enddo
c            do i=1,NBIN
c              volpart = pc3(iwa,i,2)*rmass(i,2)
c    $               -  pc3(iwa,i,3)
c              write(*,5) 'fractions(varstp):',i,
c    $           pc3(iwa,i,5)/volpart,
c    $           pc3(iwa,i,5)/(rmass(i,2)*pc3(iwa,i,2)),
c    $           pc3(iwa,i,3)/(rmass(i,2)*pc3(iwa,i,2)),
c    $           pc3(iwa,i,2),iwa
c            enddo
c            write(*,*) '-----------------------------------'
c     endif
c----------------------------------------------------------
c
c
c
c  Do post-timestep processing
c
      if( do_step ) call postep 
c
c
c  Return to caller with one timestep taken
c
      return
      end
