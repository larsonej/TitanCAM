      subroutine restart
c
c
c  @(#) init.f  Barth  Apr-2001
c  This routine restarts the model from saved data when restart file
c  cannot be used.
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
c  Declare local variables
c
        character*(25) partfile
        character*(25) gasfile
c
c  Define formats
c
    1 format(3(1pe13.6,3x))
    2 format(3(1pe8.2,3x))
    9 format(/,'Doing a restart initialization from a prev run')
   21 format(5(1pe12.2,3x))
   22 format(4(1pe12.2,3x))
   24 format(1pe14.8)
   25 format(2(1pe14.8,3x))

c
c---------------------------------------------------------------------------
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter initres'
c
c
c  Report that a restart is being done
c
      write(LUNOPRT,9)
c
c
c  Initialize all non-time-dependent quantities
c
      call initnew
c
c  Initialize time-step values from end of last run
c
      dtime = 3.6d3
      dtmax = 3.6d3
      itime = ibtime + 1
c
c  Read in time-dependent values from files
c
c  ..particle file
c
        partfile='Files/Output/ptcljc914.p'
c
      open(unit=4,file=partfile,status='old')
          read(4,24) time
      do k = 1,NZ,2
        do j = 1,NBIN

          read(4,21) rjunk,pc3(k,j,4),pc3(k,j,1),pc3(k,j,2),
     $               pc3(k,j,3)
        enddo

      enddo
      close(unit=4)

c   Average to get pc values in between 2km layers
      do k = 3,NZ,2
        do i = 1,NELEM
          do j = 1,NBIN

            pc3(k-1,j,i) = ( pc3(k-2,j,i) + pc3(k,j,i) ) / 2.d0

          enddo
        enddo
      enddo
c
c  ..gas file
c
        gasfile='Files/Output/gasjb914.p'
c
      open(unit=4,file=gasfile,status='old')
          do k = 1,NZ,2
            read(4,22) rjunk,gc(1,1,k,1),gcmix,rhc2h6i
            supsati3(k,1) = rhc2h6i/100. - ONE
      enddo
      close(unit=4)

c   Average to get gc values in between 2km layers
      do k = 3,NZ,2

        gc3(k-1,1) = ( gc3(k-2,1) + gc3(k,1) ) / 2.d0

      enddo
c
c
c  Return to caller with restart initializations complete
c
      return
      end
