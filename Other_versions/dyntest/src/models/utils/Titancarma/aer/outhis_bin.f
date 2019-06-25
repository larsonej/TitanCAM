      subroutine outhis_bin
c
c
c  @(#) outhis.f  McKie  Oct-1995
c  This routine outputs the current model state to the history file.
c  This version outputs into a traditional Fortran unformatted binary file.
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
      character*(80) text
c
c
c  Define formats
c
    1 format(/,'History output # ',i6,
     $       ' at itime: ',i6,3x,'time: ',f12.2)
c
c
c  Compute constants
c
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter outhis_bin'
c
c
c  Output history header info if this is 1st call to this routine in this run
c
      if( khist .eq. 0 )then

       text = PROGNAM // ' (' // PROGTAG // ')'

       write(LUNOHIS) text
       write(LUNOHIS) simtitle

       write(LUNOHIS) ibtime, ietime, nhist

       write(LUNOHIS)
     $  NX,       NY,       NZ,
     $  NBIN,    NELEM,
     $  NGROUP,  NGAS,     NSOLUTE

       write(LUNOHIS) ( itype(ie), ie=1,NELEM )
       write(LUNOHIS) ( igelem(ie), ie=1,NELEM )
       write(LUNOHIS) ( nelemg(ig), ig=1,NGROUP )
       write(LUNOHIS) ( ncore(ig), ig=1,NGROUP )
       write(LUNOHIS) ( ienconc(ig), ig=1,NGROUP )

       write(LUNOHIS) ( groupname(ig), ig=1,NGROUP )
       write(LUNOHIS) ( elemname(ie), ie=1,NELEM )
       write(LUNOHIS) ( gasname(igas), igas=1,NGAS )
       write(LUNOHIS) ( solname(isol), isol=1,NSOLUTE )

       write(LUNOHIS) (( r(ib,ig),ib=1,NBIN),ig=1,NGROUP )
       write(LUNOHIS) (( dr(ib,ig),ib=1,NBIN),ig=1,NGROUP )
       write(LUNOHIS) (( rmass(ib,ig),ib=1,NBIN),ig=1,NGROUP )
       write(LUNOHIS) (( dm(ib,ig),ib=1,NBIN),ig=1,NGROUP )

       write(LUNOHIS) dom_llx, dom_lly, dom_urx, dom_ury
       write(LUNOHIS) rlon0, rlat0, rlat1, rlat2, hemisph
       write(LUNOHIS) igridv, igridh
       write(LUNOHIS) gridname
       write(LUNOHIS) ((( zl(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZP1 )
       write(LUNOHIS) ((( zc(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( xc(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( yc(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( dx(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( dy(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( dz(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( xmet(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( ymet(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( zmet(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )

       write(LUNOHIS) ((( u(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
       write(LUNOHIS) ((( v(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
      endif
c
c
c  Output principal history state info to history file for current time
c
      write(LUNOHIS) itime, time
      write(LUNOHIS) ((((( pc(ix,iy,iz,ib,ie),ix=1,NX),iy=1,NY),
     $                   iz=1,NZ),ib=1,NBIN),ie=1,NELEM )
      write(LUNOHIS) (((( gc(ix,iy,iz,igas),ix=1,NX),iy=1,NY),iz=1,NZ)
     $                   ,igas=1,NGAS )
      write(LUNOHIS) ((( p(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
      write(LUNOHIS) ((( t(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
      write(LUNOHIS) ((( ptc(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZ )
      write(LUNOHIS) (((( supsatl(ix,iy,iz,igas),ix=1,NX),iy=1,NY),
     $                   iz=1,NZ),igas=1,NGAS )
      write(LUNOHIS) (((( supsati(ix,iy,iz,igas),ix=1,NX),iy=1,NY),
     $                   iz=1,NZ),igas=1,NGAS )
      write(LUNOHIS) ((( w(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZP1 )
      write(LUNOHIS) ((( dkz(ix,iy,iz),ix=1,NX),iy=1,NY),iz=1,NZP1 )
c
c
c  Return to caller with history output complete
c
      return
      end
