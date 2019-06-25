       subroutine initaer
c
c
c  @(#) initaer.f  Jensen  Oct-1995
c  This routine initializes the particle concentrations
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
      include 'globals.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter initaer'
c
c
c  Initialize particle number densities 
c  Core mass is assumed to be 100% of particle mass
c
      do ielem = 1,nelem
       ig = igelem(ielem)
       ip = ienconc(ig)
       do j = 1,nbins
        do k=1,nz
         do iy=1,ny
          do ix=1,nx

           if( ielem .eq. ip )then
c
c  Particle number concentration [#/cm^3]
c
             pc(ix,iy,k,j,ielem) = SMALL_PC

           elseif( itype(ielem) .eq. 2 )then
c
c  Core mass concentration [g/cm^3]
c
              pc(ix,iy,k,j,ielem) = pc(ix,iy,k,j,ip) 
     $                            * rmass(j,ig) / 10.

           elseif( itype(ielem) .eq. 3 )then
c
c  Second moment of core mass distribution [ (g/cm^3)^2 ]
c
              pc(ix,iy,k,j,ielem) = pc(ix,iy,k,j,ip) 
     $                            * (rmass(j,ig)/10.)**2

           endif
 
          enddo
         enddo
        enddo
       enddo
      enddo
c
c
c  Initial particle size distribution: log-normal
c
c     r0 = r( 5, 1 )
c     sig = 1.3
c     sc1 = 100.
c
c     do ie = 1, nelem
c      ig = igelem(ie)
c      ip = ienconc(ig)
c      do j = 1,nbins
c       do k=1,nz
c        do m1=1,ny
c         do m2=1,nx
c
c          if( ie .eq. 1 ) then
c            pc(m2,m1,k,j,ie) = sc1*exp(
c    $            -( (log(r(j,ig)/r0/ig))**2 ) /
c    $            ( 2.*( (log(sig))**2) ) ) /
c    $            ( sqrt(2.*PI) * r(j,ig) * log(sig) ) * dr(j,ig)
c          endif
c
c         enddo
c        enddo
c       enddo
c      enddo
c     enddo
c
c
c  Initialize first CN bins with some particles.
c
c     pc(1,1,1,1,1) = 400.
c     pc(1,1,1,2,2) = 200.
c     pc(1,1,1,3,3) = 100.
      pc(1,1,1,1,1) = 100.
      pc(1,1,1,1,2) = 100.
      pc(1,1,1,1,3) = 100.
c
c
c  Return to caller with particle concentrations initialized.
c
      return
      end
