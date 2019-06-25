      program main
c
c
c  @(#) testhis.f  McKie  Sep-1996
c  This program reads through an aerosol model history file.
c
c  Usage:  testhis  history_file_name
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Define local symbolic constants
c
      parameter(LUNI=1)
c
c
c  Declare local variables
c
      character*(80) progvers
      character*(50) hisfile
c
c
c  Define formats
c
    1 format(/,'Enter name of history file ',$)
    2 format(a)
    3 format(/,'Begin header record input ...')
    4 format(' Input ',a)
    5 format(/,'New timepoint: itime, time = ', i8, 2x,f14.4)
    6 format(/,'Begin timepoint input')
    7 format(/,'NXY, NXYZ, NXYZP1 = ', 3(1x,i8) )
    8 format(/,'Examining history file: ',a)
c
c
c  Get name of input history file to be examined from single cmd line arg
c
c <--  write(*,1)
c <--  read(*,2) hisfile
      if( iargc() .ne. 1 )then
       write(*,*) 'Error--Need input history file name on command line'
       stop
      endif
      call getarg(1, hisfile)
c
c
c  Open input history file
c
      write(*,8) hisfile
      open(LUNI,file=hisfile,status='old',form='unformatted')
c
c
c  Input history header info
c
      write(*,3)
    
      write(*,4) 'progvers'
      read(LUNI) progvers
      write(*,*) progvers

      write(*,4) 'simtitle'
      read(LUNI) simtitle
      write(*,*) simtitle

      write(*,4) 'ibtime, ietime, nhist'
      read(LUNI) ibtime, ietime, nhist
      write(*,*) ibtime, ietime, nhist

      write(*,4) 'NX, NAXNY, NZ, NBIN, NELEM, ' // 
     $  'NGROUP, NGAS, NSOLUTE'
      read(LUNI)
     $  iNX,    iNY,    iNZ,
     $  iNBIN, iNELEM,
     $  iNGROUP,  iNGAS,  iNSOLUTE
      write(*,*) iNX, iNY, iNZ, iNBIN, iNELEM,
     $  iNGROUP, iNGAS, iNSOLUTE

      write(*,4) 'NX, NY, NZ, NBIN, NELEM, NGROUP, NGAS, NSOLUTE'
      read(LUNI)
     $  NX,       NY,       NZ,
     $  NBIN,    NELEM,
     $  NGROUP,  NGAS,     NSOLUTE
      write(*,*) NX, NY, NZ, NBIN, NELEM, NGROUP, NGAS, NSOLUTE

      write(*,4) 'itype'
      read(LUNI) ( itype(ie), ie=1,NELEM )
      write(*,4) 'igelem'
      read(LUNI) ( igelem(ie), ie=1,NELEM )
      write(*,4) 'nelemg'
      read(LUNI) ( nelemg(ig), ig=1,NGROUP )
      write(*,4) 'ncore'
      read(LUNI) ( ncore(ig), ig=1,NGROUP )
      write(*,4) 'ienconc'
      read(LUNI) ( ienconc(ig), ig=1,NGROUP )

      write(*,4) 'groupname'
      read(LUNI) ( groupname(ig), ig=1,NGROUP )
      write(*,4) 'elemname'
      read(LUNI) ( elemname(ie), ie=1,NELEM )
      write(*,4) 'gasname'
      read(LUNI) ( gasname(igas), igas=1,NGAS )
      write(*,4) 'solname'
      read(LUNI) ( solname(isol), isol=1,NSOLUTE )

      write(*,4) 'r'
      read(LUNI) (( r(ib,ig),ib=1,NBIN),ig=1,NGROUP )
      write(*,4) 'dr'
      read(LUNI) (( dr(ib,ig),ib=1,NBIN),ig=1,NGROUP )
      write(*,4) 'rmass'
      read(LUNI) (( rmass(ib,ig),ib=1,NBIN),ig=1,NGROUP )
      write(*,4) 'dm'
      read(LUNI) (( dm(ib,ig),ib=1,NBIN),ig=1,NGROUP )
      write(*,4) 'zc3'
      read(LUNI) ( zc3(i),i=1,NXYZ )
c
c
c  Compute constants that are functions of input values
c
      write(*,7) NXY, NXYZ, NXYZP1
c
c
c  Input next timepoint of history data until eof
c
      write(*,6)
 4100 continue
       read(LUNI,end=7100) itime, time
       write(*,5) itime, time
       write(*,4) 'pc3'
       read(LUNI) ((( pc3(i,ib,ie),i=1,NXYZ),ib=1,NBIN),ie=1,NELEM )
       write(*,4) 'gc3'
       read(LUNI) (( gc3(i,igas),i=1,NXYZ),igas=1,NGAS )
       write(*,4) 'p3'
       read(LUNI) ( p3(i),i=1,NXYZ )
       write(*,4) 't3'
       read(LUNI) ( t3(i),i=1,NXYZ )
       write(*,4) 'ptc3'
       read(LUNI) ( ptc3(i),i=1,NXYZ )
       write(*,4) 'supsati3'
       read(LUNI) (( supsati3(i,igas),i=1,NXYZ),igas=1,NGAS )
       write(*,4) 'w3'
       read(LUNI) ( w3(i),i=1,NXYZP1 )
       write(*,4) 'dkz3'
       read(LUNI) ( dkz3(i),i=1,NXYZP1 )
       goto 4100
 7100 continue
      write(*,*) 'EOF found on input history file'
c
c
c  Close input file
c
      close(unit=LUNI)
c
c
c  Terminate normally
c
      end
