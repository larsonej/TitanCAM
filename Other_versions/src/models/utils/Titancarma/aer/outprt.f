      subroutine outprt
c
c
c  @(#) outprt.f  McKie  Oct-1995
c  This routine outputs information about the current timestep
c  to the output print file.
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
      dimension totn(NX,NY,NZ,NGROUP), rn(NX,NY,NZ,NGROUP), 
     $          rsig(NX,NY,NZ,NGROUP)

      character*(6) cwave
      logical rep_gas, rep_mass, rep_logfit, rep_part
c
c
c  Define formats
c
    1 format('Timestep itime: ',i10,3x,'time: ',2(1pe14.6))
    2 format(/,'Total particle mass, total solute mass [g]: ',1p,2e14.6)
    3 format(/,'Particle size distributions at (ix,iy,iz) = ',
     $       '(', i4, ',', i4, ',', i4, ')' )
    4 format(i6,3x,i6,3x,1p,e12.2,3x,1p,e13.6)
    5 format(/,'Gas concentrations for ',a,' at (ix,iy) = ',
     $  '(', i4, ',', i4, ')', //,
     $  a3, 1x, 4(a11,4x), /)
    6 format(i3,1x,1p,3(e11.3,4x),0p,f11.3)
    7 format(/,'Water vapor, condensed, total [g]: ',1p,3e14.6)
    8 format(/,'Total particle number: ',1p,1e14.6)
    9 format(/,'Total CN mass, total {in,}volatile core mass [g]: ',
     $  1p,3e14.6)
   10 format(/, 2x, a, /)
   11 format(18x, a17, 18x, a19, /, 18x, 33('-'),3x,32('-'))
   12 format(a3, a12, 2(a12,a12,a12), /)
   13 format(i3, 1p, e12.3, 6e12.3, /, (3x, 12x, 1p, 6e12.3) )
   14 format(3x,'Ielem',5x,'ibin',4x,'r (microns)',3x,'N or fraction')
   15 format(/,a, ' at (ix,iy) = (', i4, ',', i4, ')' ,/)
   16 format(/,'Albedo (',a,', ',a,') [%]: x across, y down',//,20i10)
   17 format(i4,20f10.1)
   18 format(/,'Optical Depth (',a,', ',a,'): x across, y down',//,
     $  20i10)
   19 format(i4,20(1pe10.2))
   20 format('Detailed radiative transfer output for ix, iy = ',2i4)
   21 format(5(1pe12.2,3x))
   22 format(4(1pe12.2,3x))
   23 format(f12.2,i4,f6.2,i4,8(1pe10.2,3x))
   24 format(1pe14.8)
   25 format(2(1pe14.8,3x))
c
c
        character*(50) partfile
        character*(50) gasfile
c       integer thirdbit
c
c
c---------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter outprt'
c
c
c  Report current timestep index and current simulation time
c
c     call prtsep
      write(LUNOPRT,1) itime, time
c     write(LUNOTME,24) time
c
c
      write(*,1) itime,time,time/60.**2/24.     !days
c     write(LUNOTEMP,1) itime,time,time/60.**2/24.      !days

c
c   ..Write out particle & gas data to files:
      if( .true.) then
c  
c   ..Initial file names:
c
c   ..First file name tag to be appended to initial file name
         iuno=iuno+1
         if (iuno.gt.9) then
                iuno=0
c   ..Second file name tag
c         If(thirdbit.lt.1) itre = 0
                iduo=iduo+1
                if (iduo.gt.9) then
                    itre=itre+1
                    iduo=0
c                   thirdbit = 2
                endif
         endif
c
c  ..Construct file names..
c
        partfile='Files/Output/ptcl'//ext//char(itre+48)//
     $   char(iduo+48)//char(iuno+48)//'.p'
        gasfile='Files/Output/gas'//ext//char(itre+48)//
     $   char(iduo+48)//char(iuno+48)//'.p'
c
c   ..Report file name writen to...
c
         print*,"Writing to:"
         print*,partfile
         print*,gasfile
c
c  ..Write out particle distributions:  N(r,z)
c
        OPEN(unit=2,file=partfile,STATUS='UNKNOWN')

        write(2,25) time,dtime
      do k = 1,NZ
        do j = 1,NBIN

c         rcore=(3.*vcore(k,j,4)/(4.*PI))**(1./3.)
c         if(rcore = 'NaN') rcore=1.e-30
c       coreshell(k,j,2)=rmshell(k,j,2)/
c    $        (4.*PI*rhoelem(2)*r(j,2)**2) + rcore

          write(2,21) zl(1,1,k),pc3(k,j,4),pc3(k,j,1),
     $                  pc3(k,j,2),pc3(k,j,3)
          enddo

      enddo
        CLOSE(unit=2)
c
c  ..Write out gas and RH vs. altitude
c
      OPEN(unit=3,file=gasfile,STATUS='UNKNOWN')
          do k = 1,NZ
           gcmix=gc(1,1,k,1)/rhoa(1,1,k)*1.e6
           relhi = 100.*(1. + supsati3(k,1))
c        write(3,22) zl3(k),gc(1,1,k,1),gcmix,relhi
         write(3,22) zl3(k),gc(1,1,k,1),t3(k),relhi
      enddo
      CLOSE(unit=3)
c
      endif !write output to ptcl and gas files
c
c
c
c  Return to caller with timestep info output to print file
c
      return
      end
