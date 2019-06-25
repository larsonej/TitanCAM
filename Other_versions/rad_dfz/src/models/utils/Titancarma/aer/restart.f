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
        character*(25) fluxfile

        dimension ftopsmall(NBIN)
c
c  Define formats
c
    1 format(3(1pe13.6,3x))
    2 format(3(1pe8.2,3x))
    3 format(2(i4),5(3x,e13.7))
    4 format(i4,6(3x,e13.7))
    5 format(2(i4),8(3x,e13.7))
    6 format(2(i4),11(3x,e13.7))
    7 format(2(i4),12(3x,e13.7))
    8 format(e13.7)
    9 format(/,'Doing a restart initialization from a prev run')
   10 format(2(i4),4(3x,e13.7))
   11 format(2(i4),7(3x,e13.7))
   12 format(i4,10(3x,e13.7))
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
c  Initialize all non-time-dependent quantities
c
      call initnew
c
c  Initialize time-step values from end of last run
c
      dtime = dtmin
      itime = ibtime + 1
c
c  Read in time-dependent values from files
c
       !original output was netcdf format
       !converted to ascii using idl
c
c  min nucleated particle
c
        inucmin(1) = 20
        inucmin(2) = 10
        inucmin(3) = 13
c
c  ..particle file
c
        partfile='Files/ptclae187.p'
        time=108770.d0*60.d0**2*24.d0
cccccccccccccccc

c  For temperature perturbation cases, allow model to run for 
c  several days before inititating sinusoidal T wave

        tpstart = time + 8.64d4

c.....................Possibly change endtime.........................
c  Reset endtime to be (start) time + x year (+1 day)
ccc     write(LUNOPRT,*) 'RESTART: Resetting endtime to 2 years later!'
ccc     endtime = time + 2.d0*YEAR + 8.64d4

cc      write(LUNOPRT,*) 
cc   $        'RESTART: Resetting endtime to 100 days (+1 year) later!'
cc      endtime = time + 1.d0*YEAR + 100.d0*8.64d4

c        write(LUNOPRT,*) 'RESTART: Resetting endtime to 20 days later!'
c        endtime = time + 20.d0 * 8.64d4

cc      write(LUNOPRT,*) 'RESTART: Resetting endtime to 125 days later!'
cc      endtime = time + 125.d0 * 8.64d4

cc      write(*,*) 'RESTART: Resetting ietime to 500 later!!'
cc      ietime = ibtime + 500
c.......................................................................

       open(unit=4,file=partfile,status='old')
       do k = 1,NZ
        do j = 1,NBIN

cccccc 11 element format
         if( NELEM .eq. 11 ) then
          read(4,6) junk,junk,pc3(k,j,1),pc3(k,j,2),
     $               pc3(k,j,3),pc3(k,j,4),pc3(k,j,5),
     $               pc3(k,j,6),pc3(k,j,7),pc3(k,j,8),
     $               pc3(k,j,9),pc3(k,j,10),pc3(k,j,11)
cccccc 8 element format
         else if( NELEM .eq. 8 ) then
          read(4,5) junk,junk,pc3(k,j,1),pc3(k,j,2),
     $               pc3(k,j,3),pc3(k,j,4),pc3(k,j,5),
     $               pc3(k,j,6),pc3(k,j,7),pc3(k,j,8)
cccccc 12 element format
         else if( NELEM .eq. 12 ) then
          read(4,7) junk,junk,pc3(k,j,1),pc3(k,j,2),
     $               pc3(k,j,3),pc3(k,j,4),pc3(k,j,5),
     $               pc3(k,j,6),pc3(k,j,7),pc3(k,j,8),
     $               pc3(k,j,9),pc3(k,j,10),pc3(k,j,11),
     $               pc3(k,j,12)
cccccc 7 element format
         else if( NELEM .eq. 7 ) then
          read(4,11) junk,junk,pc3(k,j,1),pc3(k,j,2),
     $               pc3(k,j,3),pc3(k,j,4),pc3(k,j,5),
     $               pc3(k,j,6),pc3(k,j,7)
cccccc 5 element format
         else if( NELEM .eq. 5 ) then
          read(4,3) junk,junk,pc3(k,j,1),pc3(k,j,2),
     $               pc3(k,j,3),pc3(k,j,4),pc3(k,j,5)
cccccc 4 element format
         else if( NELEM .eq. 4 ) then
          read(4,10) junk,junk,pc3(k,j,1),pc3(k,j,2),
     $               pc3(k,j,3),pc3(k,j,4)
         else
          write(*,*) 'Wrong number of elements',NELEM
          stop
         endif

        enddo !bins
       enddo  !alt
       close(unit=4)
c
c  Don't allow particle concentrations to get too small.
c
      do ixyz = 1,NXYZ
        do ibin = 1,NBIN
          do ielem = 1,NELEM
            call smallconc(ibin,ielem)
          enddo
        enddo
      enddo

c
c  ..gas file
c
        gasfile='Files/gasae187.p'

       open(unit=4,file=gasfile,status='old')
       do k = 1,NZ
         if( NGAS.eq.2 ) then
            read(4,12) junk,gc3(k,1),gc3(k,2),rhc2h6i,
     $                 rhch4i,t3(k),rhch4l,gfluxc2h6,
     $                 gfluxch4,puddle(1),puddle(2)
            supsati3(k,1) = rhc2h6i/100.d0 - ONE
            supsati3(k,2) = rhch4i/100.d0 - ONE
            supsatl3(k,2) = rhch4l/100.d0 - ONE
            if( k.eq.1 ) fbotevap(1) = gfluxc2h6
            if( k.eq.1 ) fbotevap(2) = gfluxch4
         elseif( NGAS.eq.1 ) then
            read(4,4) junk,gc3(k,1),
     $                 rhch4i,t3(k),rhch4l,
     $                 gfluxch4,puddle(1)
            supsati3(k,1) = rhch4i/100.d0 - ONE
            supsatl3(k,1) = rhch4l/100.d0 - ONE
            if( k.eq.1 ) fbotevap(1) = gfluxch4
         else
           write(*,*) '<restart> Wrong NGAS',NGAS
           stop
         endif
       enddo
       close(unit=4)

c      write(*,*) 'No gas read in !!!!!!!!!!!!'
c
c  ..possibly read in particle flux
c 
c      Save fluxes for 'small' tholin particles (loaded in <initaer>),
c      but offset to match new radius bins
cc        do j=11,NBIN
c           ftopsmall(j-10) = ftoppart(1,j,1)
c         enddo
c         do j=1+NBIN-10,NBIN
c           ftopsmall(j) = ZERO
c         enddo
c      fluxfile = 'Files/fluxt_062.p'
c      open(unit=4,file=fluxfile,status='old')
c         do j=1,NBIN
c          if(j .le. 48) then
c           read(4,8) ftoppart(1,j,1)
c           ftoppart(1,j,1) = ftoppart(1,j,1) 
c    $                        + ftopsmall(j) 
c          else
c           ftoppart(1,j,1) = ZERO
c          endif
c         enddo
c      close(unit=4)
c      write(*,*) 'Reading in new tholin fluxes !!!!!!!!!!!!'

c
c
c  Return to caller with restart initializations complete
c
      return
      end
