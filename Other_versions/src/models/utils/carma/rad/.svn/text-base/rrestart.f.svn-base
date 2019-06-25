      SUBROUTINE rrestart(iswitch)

      include 'globrad.h'
 
C     *************************************************************
c     IF (iswitch == 0) THEN read common blocks
c     ELSE write them
c    
c     The success of this routine depends upon coordinated managment
c     of the common blocks in file RADCOMM.H.  The algorithm is
c     based upon the contiguous storage of common block variables,
c     and requires that the first and last variables in each block
c     be identified in this subroutine, in order to determine their
c     addresses in memory.  The difference between the two addresses
c     determines the number of words in the common block, which are
c     read (or written) through an equivalenced array.
c
c     This subroutine is called by RESTART, where the file is
c     opened and closed.
c
c     It is assumed that a system library routine LOC exists, and
c     returns the location of a variable in memory.
C     *************************************************************

c  IDUMVAL is a dummy value used for error checking.
      parameter ( idumval = 69 )

      integer idum1(2), idum2(2),
     1        idum3(2), idum4(2), idum5(2), idum6(2), idum7(2)

      equivalence (idum1(1),irbeg1),
     2            (idum2(1),irbeg2),
     3            (idum3(1),irbeg3),
     4            (idum4(1),irbeg4),
     5            (idum5(1),irbeg5),
     6            (idum6(1),irbeg6),
     7            (idum7(1),irbeg7)

      data irbeg1 / idumval /,
     2     irbeg2 / idumval /,
     3     irbeg3 / idumval /,
     4     irbeg4 / idumval /,
     5     irbeg5 / idumval /,
     6     irbeg6 / idumval /,
     7     irbeg7 / idumval /

      data irend1 / idumval /,
     2     irend2 / idumval /,
     3     irend3 / idumval /,
     4     irend4 / idumval /,
     5     irend5 / idumval /,
     6     irend6 / idumval /,
     7     irend7 / idumval /

c IUNITS is the unit of storage.

      iunits = loc(idum1(2)) - loc(idum1(1))

c These are number of words to write for each block.

      nw1 = (loc(irend1) - loc(irbeg1))/iunits + 1
      nw2 = (loc(irend2) - loc(irbeg2))/iunits + 1
      nw3 = (loc(irend3) - loc(irbeg3))/iunits + 1
      nw4 = (loc(irend4) - loc(irbeg4))/iunits + 1
      nw5 = (loc(irend5) - loc(irbeg5))/iunits + 1
      nw6 = (loc(irend6) - loc(irbeg6))/iunits + 1
      nw7 = (loc(irend7) - loc(irbeg7))/iunits + 1

c Error checking.

      if (idum1(1).ne.idumval .or. idum1(nw1).ne.idumval)
     1  stop 'block 1'
      if (idum2(1).ne.idumval .or. idum2(nw2).ne.idumval)
     1  stop 'block 2'
      if (idum3(1).ne.idumval .or. idum3(nw3).ne.idumval)
     1  stop 'block 3'
      if (idum4(1).ne.idumval .or. idum4(nw4).ne.idumval)
     1  stop 'block 4'
      if (idum5(1).ne.idumval .or. idum5(nw5).ne.idumval)
     1  stop 'block 5'
      if (idum6(1).ne.idumval .or. idum6(nw6).ne.idumval)
     1  stop 'block 6'
      if (idum7(1).ne.idumval .or. idum7(nw7).ne.idumval)
     1  stop 'block 7'

      if (iswitch.eq.0) then
        read(3) (idum1(j),j=1,nw1)
        read(3) (idum2(j),j=1,nw2)
        read(3) (idum3(j),j=1,nw3)
        read(3) (idum4(j),j=1,nw4)
        read(3) (idum5(j),j=1,nw5)
        read(3) (idum6(j),j=1,nw6)
        read(3) (idum7(j),j=1,nw7)
      else
        write(3) (idum1(j),j=1,nw1)
        write(3) (idum2(j),j=1,nw2)
        write(3) (idum3(j),j=1,nw3)
        write(3) (idum4(j),j=1,nw4)
        write(3) (idum5(j),j=1,nw5)
        write(3) (idum6(j),j=1,nw6)
        write(3) (idum7(j),j=1,nw7)
      endif

      RETURN
      END
