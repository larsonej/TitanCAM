c
c  @(#) dumnetcdf.f  McKie  Sep-1999
c
c
c  This is a set of dummy netcdf routines called in carma source.
c  This f77 source is intended to be compiled into an object
c  library, to be linked in the case when netcdf history
c  output is not turned on (i.e. regular Fortran binary history
c  output is turned on).  These routines allow the program linking
c  to complete without unsatisfied externals in that case.  If any of the
c  netcdf routines used by carma are accidentally called at run time after
c  this library is used at link time, then an error message is written
c  to standard output, and the program is aborted.
c
c  Note that only those netcdf routines known to be used in
c  carma source are defined here.
c
c
      function nccre(c1,i2,ier)
       integer nccre
       character*(*) c1
       integer i2
       integer ier
       nccre = 1
       call nc_abort
      end
c
      function ncddef(i1,c2,i3,ier)
       integer ncddef
       integer i1
       character*(*) c2
       integer i3
       integer ier
       ncddef = 1
       call nc_abort
      end
c
      function ncdid(i1,c2,i3,i4,i5,ier)
       integer ncdid
       integer i1
       character*(*) c2
       integer i3
       integer i4
       integer i5
       integer ier
       ncdid = 1
       call nc_abort
      end
c
      function ncvdef(i1,c2,i3,i4,i5,ier)
       integer ncvdef
       integer i1
       character*(*) c2
       integer i3
       integer i4
       integer i5
       integer ier
       ncvdef = 1
       call nc_abort
      end
c
      function ncvid(i1,c2,ier)
       integer ncvid
       integer i1
       character*(*) c2
       integer ier
       ncvid = 1
       call nc_abort
      end
c
      subroutine ncvpt(i1,i2,i3,i4,i5,ier)
       integer i1
       integer i2
       integer i3
       integer i4
       integer i5
       integer ier
       call nc_abort
      end
c
      subroutine ncvptc(i1,i2,i3,i4,c5,i6,ier)
       integer i1
       integer i2
       integer i3
       integer i4
       character*(*) c5
       integer i6
       integer ier
       call nc_abort
      end
c
      subroutine ncsnc(i1,ier)
       integer i1
       integer ier
       call nc_abort
      end
c
      subroutine ncclos(i1,ier)
       integer i1
       integer ier
       call nc_abort
      end
c
      subroutine ncendf(i1,ier)
       integer i1
       integer ier
       call nc_abort
      end
c
      subroutine ncaptc(i1,i2,c3,i4,i5,c6,ier)
       integer i1
       integer i2
       character*(*) c3
       integer i4
       integer i5
       character*(*) c6
       integer ier
       call nc_abort
      end
c
      subroutine nc_abort
       write(*,*) 'Dummy netcdf library routine called unexpectedly.'
       write(*,*) 'Either turn off netcdf history output,'
       write(*,*) 'Or specify a real netcdf library in top Makefile.'
       write(*,*) 'Aborting.'
       stop
      end
