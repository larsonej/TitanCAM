      program main
c
c
c  @(#) tprint2d.f  McKie  Jul-1997
c  This program tests the print2d routine.
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Define local symbolic constants
c
      parameter(NX=25)
      parameter(NY=14)
      parameter(NZ=12)
c
c  Declare local variables
c
      dimension f(NX,NY,NZ)
      character*(60) title
c
c
c  Open output file
c
      lun = 1
      open(unit=lun,file='tprint2d.p',status='unknown')
c
c
c  Define elements of <f> that indicate its indices
c
      do k=1,NZ
       do j=1,NY
        do i=1,NX
         f(i,j,k) = float( i*10000 + j*100 + k )
        enddo
       enddo
      enddo
c
c
c  Print a x vs y subset of f
c
      nxb = 3
      nxe = 21
      nyb = 5
      nye = 12
      nzb = 2
      nze = 2
      title = 'TITLE IS X VS Y,  x:[3,21] y:[5,12] z:[2,2]'
      call print2d
     $ (
     $  f,NX,NY,
     $  nxb,nxe, nyb,nye, nzb,nze,
     $  lun, title
     $ ) 
c
c
c  Print a x vs z subset of f
c
      nxb = 3
      nxe = 19
      nyb = 5
      nye = 5
      nzb = 10
      nze = 2
      title = 'TITLE IS X VS Y,  x:[3,19] y:[5,5] z:[10,2]'
      call print2d
     $ (
     $  f,NX,NY,
     $  nxb,nxe, nyb,nye, nzb,nze,
     $  lun, title
     $ ) 
c
c
c  Print a y vs z subset of f
c
      nxb = 11
      nxe = 11
      nyb = 8
      nye = 1
      nzb = 7
      nze = 2
      title = 'TITLE IS X VS Y,  x:[11,11] y:[8,1] z:[7,2]'
      call print2d
     $ (
     $  f,NX,NY,
     $  nxb,nxe, nyb,nye, nzb,nze,
     $  lun, title
     $ ) 
c
c
c  Close output file
c
      close(unit=lun)
c
c
c  Terminate normally
c
      stop
      end
