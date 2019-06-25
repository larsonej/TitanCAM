      subroutine print2d
     $  (
     $   f, idim1,idim2,
     $   ixb,ixe, iyb,iye, izb,ize,
     $   lun, title
     $  )
c
c
c  @(#) print2d.f  McKie  Jul-1997
c  This routine prints a 2-D subset of a 3-D array.
c
c  Note:
c    1st index of f is x direction
c    2nd index of f is y direction
c    3rd index of f is z direction
c
c  Input:
c    f(idim1,idim2,*) = field to be printed
c    idim1,idim2      = 1st,2nd dimensions of f
c    ixb,ixe          = begin,end x index range
c    iyb,iye          = begin,end y index range
c    izb,ize          = begin,end z index range
c    lun              = output lun assumed open
c    title            = title text
c
c  Note:
c    One of the ranges nxb,nxe, nyb,nye, nzb,nze must be 1 index.
c    If i?b > i?e, then that index is reverse printed.
c    The printout "plane" is one of:  (x vs y), (x vs z), (y vs z).
c    Field is printed as NDIGITS digit integers, scaled by a multiplier.
c    Must change formats 1 & 2 if NDIGITS change (E.g. i7 for NDIGITS=7)
c
c
c  Output:  none
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
      integer idim1,idim2
      dimension f(idim1,idim2,*)
      integer ixb,ixe
      integer iyb,iye
      integer izb,ize
      integer lun
      character*(*) title
c
c
c  Define local symbolic constants
c
      parameter(NXPAGE=10)
      parameter(NDIGITS=7)
      parameter(BIG=1.d+35)
      parameter(SMALL=1.d-30)
      parameter(I_XVSY=1)
      parameter(I_XVSZ=2)
      parameter(I_YVSZ=3)
c
c
c  Declare local variables
c
      character*(6) orient
c
c
c  Define formats for this routine
c
    1 format(1x,i2,':',18(1x,i7))
    2 format(1x,3x,18(1x,i7))
    4 format(/,1x,50('-'))
    5 format(' ')
    7 format(1x,a,3x,'(',a,')',2x,a,'=',i3)
   10 format(1x,'[min,max]: [',1p,e10.3,',',e10.3,']',/,
     $  1x,'Field was multiplied by ',e9.1,' for printing',0p,/,
     $  1x,'Page',i2,' of ',i2)
c
c
c  Copy requested index ranges to local variable
c
      ib = ixb
      ie = ixe
      jb = iyb
      je = iye
      kb = izb
      ke = ize
c
c
c  Define index increments in 3 dimensions
c
      if( ib .le. ie )then
       idi = 1
      else
       idi = -1
      endif
      if( jb .le. je )then
       idj = 1
      else
       idj = -1
      endif
      if( kb .le. ke )then
       idk = 1
      else
       idk = -1
      endif
c
c
c  Decide if this is a (x vs y), (x vs z), or (y vs z) print orientation
c   Do a (x vs y) orientation with lowest z index if none of x,y,z is 1 index range.
c
      if( kb .eq. ke )then
       ior = I_XVSY
       orient = 'x vs y'
      else if( jb .eq. je )then
       ior = I_XVSZ
       orient = 'x vs z'
      else if( ib .eq. ie )then
       ior = I_YVSZ
       orient = 'y vs z'
      else
       ior = I_XVSY
       ke = kb  
       orient = 'x vs y'
      endif
c
c
c  Compute min, max of field
c
      fmin = BIG
      fmax = -fmin
      do k=kb,ke,idk
       do j=jb,je,idj
        do i=ib,ie,idi
         fmin = min( fmin, f(i,j,k) )
         fmax = max( fmax, f(i,j,k) )
        enddo
       enddo
      enddo
c
c
c  Compute a power of 10 multiplier (to get 6 significant digit integers)
c
      call getmult(fmin,fmax,NDIGITS,SMALL, rmult)
c
c
c--------------------------------------------------------------------------
c
c  Do (x vs y) orientation array printout if requested
c
      if( ior .eq. I_XVSY )then
       write(lun,4)
       nx = abs( ie - ib ) + 1
       npage = nx / NXPAGE
       if( ( npage * NXPAGE ) .lt. nx ) npage = npage + 1
       k = kb
       do m=1,npage
        i1 = ( m - 1 ) * NXPAGE + ib
        i2 = min( i1 + NXPAGE - 1 , ie )
        write(lun,5)
        write(lun,7) title, orient, 'k',k
        write(lun,10) fmin,fmax,rmult,m,npage
        write(lun,5)
        write(lun,2) (i,i=i1,i2,idi)
        write(lun,5)
        do j=jb,je,idj
         write(lun,1) j,(nint(rmult*f(i,j,k)),i=i1,i2,idi)
        enddo
        write(lun,5)
        write(lun,2) (i,i=i1,i2,idi)
       enddo
c
c
c--------------------------------------------------------------------------
c
c  Do (x vs z) orientation array printout if requested
c
      else if( ior .eq. I_XVSZ )then
       write(lun,4)
       nx = abs( ie - ib ) + 1
       npage = nx / NXPAGE
       if( ( npage * NXPAGE ) .lt. nx ) npage = npage + 1
       j = jb
       do m=1,npage
        i1 = ( m - 1 ) * NXPAGE + ib
        i2 = min( i1 + NXPAGE - 1 , ie )
        write(lun,5)
        write(lun,7) title, orient, 'j',j
        write(lun,10) fmin,fmax,rmult,m,npage
        write(lun,5)
        write(lun,2) (i,i=i1,i2,idi)
        write(lun,5)
        do k=kb,ke,idk
         write(lun,1) k,(nint(rmult*f(i,j,k)),i=i1,i2,idi)
        enddo
        write(lun,5)
        write(lun,2) (i,i=i1,i2,idi)
       enddo
c
c
c--------------------------------------------------------------------------
c
c  Do (y vs z) orientation array printout if requested
c
      else if( ior .eq. I_YVSZ )then
       write(lun,4)
       ny = abs( je - jb ) + 1
       npage = ny / NXPAGE
       if( ( npage * NXPAGE ) .lt. ny ) npage = npage + 1
       i = ib
       do m=1,npage
        j1 = ( m - 1 ) * NXPAGE + jb
        j2 = min( j1 + NXPAGE - 1 , je )
        write(lun,5)
        write(lun,7) title, orient, 'i',i
        write(lun,10) fmin,fmax,rmult,m,npage
        write(lun,5)
        write(lun,2) (j,j=j1,j2,idj)
        write(lun,5)
        do k=kb,ke,idk
         write(lun,1) k,(nint(rmult*f(i,j,k)),j=j1,j2,idj)
        enddo
        write(lun,5)
        write(lun,2) (j,j=j1,j2,idj)
       enddo
      endif
c
c
c--------------------------------------------------------------------------
c
c  Return to caller with array subset printed
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine getmult(fmin,fmax,ndigits,tiny, rmult)
c
c
c  @(#) getmult.f  McKie  Jul-1988
c  This routine computes a power of 10 multiplier which scales a range
c  into integers with a specified number of significant digits.
c
c  Input:
c          fmin,fmax = range min,max
c            ndigits = max # digits in integers
c               tiny = small number
c
c  Output:
c              rmult = the multiplier
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
      absmax = max(abs(fmin),abs(fmax))
      if( absmax .lt. tiny )then
       rmult = 0.
      else
       rmult = log10(absmax)
       iexp = int(rmult)
       rmult = 10.**(ndigits - 2 - iexp)
      endif
      return
      end
