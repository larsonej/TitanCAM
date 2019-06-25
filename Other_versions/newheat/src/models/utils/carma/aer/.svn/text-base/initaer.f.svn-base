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
      include 'globaer.h'
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
      do ielem = 1,NELEM
       ig = igelem(ielem)
       ip = ienconc(ig)
       do j = 1,NBIN
        do iz = 1,NZ
         do iy = 1,NY
          do ix = 1,NX

           if( ielem .eq. ip )then
c
c  Particle number concentration [#/cm^3]
c
             pc(ix,iy,iz,j,ielem) = 0.
c            pc(ix,iy,iz,j,ielem) = SMALL_PC

           elseif( itype(ielem) .eq. I_COREMASS )then
c
c  Core mass concentration [g/cm^3]
c
             pc(ix,iy,iz,j,ielem) = 0.
c            pc(ix,iy,iz,j,ielem) = pc(ix,iy,iz,j,ip) *
c    $                              rmass(j,ig) * FIX_COREF

           elseif( itype(ielem) .eq. I_CORE2MOM )then
c
c  Second moment of core mass distribution [ (g/cm^3)^2 ]
c
             pc(ix,iy,iz,j,ielem) = 0.
c            pc(ix,iy,iz,j,ielem) = pc(ix,iy,iz,j,ip) *
c    $                              (rmass(j,ig)*FIX_COREF)**2

           endif

c          pc(ix,iy,iz,j,ielem) = small_val(j,ielem)
 
          enddo
         enddo
        enddo
       enddo
      enddo
c
c
c  Initial particle distribution: log-normal size distribution 
c  for first particle group (which has only one particle element)
c  in a single column
c
      ig = 1
      ie = ienconc(ig)

      do ix = 1,NX
      do iy = 1,NY
      do iz = 1,NZ
c
c
c  Calculate index <jz> that increases with altitude
c
        if( igridv .eq. I_CART )then
           jz = iz
        else
           jz = NZ + 1 - iz
        endif
c
c  Log-normal parameters:
c  
c    r0   = number mode radius
c    rsig = geometric standard deviation
c    totn = total number concentration
c
        r0   = 1.e-5
        rsig = 1.2
        totn = 100.
c
c
c  Smooth gradient of <totn> 
c
c       if( NZ .gt. 1 )then
c         totn = 1.e3 / 10.**( (jz-1) / float(NZ-1) )
c       else
c         totn = 1.e2
c       endif
c
c
c  Particles in the top layer only
c
c       if( jz .eq. NZ )then
c         totn = 1.e3
c       else
c         totn = 0.
c       endif
c
c
c  Adjust prefactor to yield particle number concentration <ntot>
c
        sum = 0.
        do j = 1,NBIN
          arg1 = dr(j,ig) / ( sqrt(2.*PI) * r(j,ig) * log(rsig) ) 
          arg2 = -log( r(j,ig) / r0 )**2 / ( 2.*log(rsig)**2 )
          sum  = sum + arg1 * exp( arg2 )
        enddo
        totn = totn / sum

        do j = 1,NBIN

          arg1 = totn * dr(j,ig) / ( sqrt(2.*PI) * r(j,ig) * log(rsig) ) 
          arg2 = -log( r(j,ig) / r0 )**2 / ( 2.*log(rsig)**2 )

c         pc(ix,iy,iz,j,ie) = max( arg1 * exp( arg2 ), SMALL_PC )
          pc(ix,iy,iz,j,ie) = arg1 * exp( arg2 )

        enddo

      enddo
      enddo
      enddo
c
c
c  Log-normal parameters for hydrometeors
c
c
c      do ix = 1,NX
c      do iy = 1,NY
c      do iz = 1,NZ
c
c        ig = 2
c        ie = ienconc(ig)
c        r0   = 50.e-4
c        rsig = 1.4
c        totn = 1000.
c
c        sum = 0.
c        do j = 1,NBIN
c          arg1 = dr(j,ig) / ( sqrt(2.*PI) * r(j,ig) * log(rsig) ) 
c          arg2 = -log( r(j,ig) / r0 )**2 / ( 2.*log(rsig)**2 )
c          sum  = sum + arg1 * exp( arg2 )
c        enddo
c        totn = totn / sum
c
c        do j = 1,NBIN
c
c          arg1 = totn*dr(j,ig) / ( sqrt(2.*PI)*r(j,ig)*log(rsig) ) 
c          arg2 = -log( r(j,ig) / r0 )**2 / ( 2.*log(rsig)**2 )
c
c          pc(ix,iy,iz,j,ie) = max( arg1 * exp( arg2 ), SMALL_PC )
c          pc(ix,iy,iz,j,ie+1) = pc(ix,iy,iz,j,ie) *
c    $                           rmass(j,ig) * FIX_COREF
c
c        enddo
c
c      enddo
c      enddo
c      enddo
c
c
c  Initialize a single bin with some particles.
c
c      pc(1,1,NZ,1,1) = 100.
c      pc(1,1,1,21,4) = 20.
c      pc(1,1,1,17,3) = 10.*rmass(20,2)*FIX_COREF
c      pc(1,1,1,21,5) = 20.*rmass(21,3)*FIX_COREF
c
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
c
c  Specify fluxes at top and bottom of model
c  [#/cm^2/s for particle number types, g/cm^2/s for core mass, etc.]
c
      do ixy = 1,NXY
        do ie = 1, NELEM
          do j = 1,NBIN
            ftoppart(ixy,j,ie) = 0.
            fbotpart(ixy,j,ie) = 0.
          enddo
        enddo
      enddo
c
c
c  Scale particle concentrations and boundary fluxes from
c  cartesian coordinates to coordinate system specified by <igrid>
c
c
c  Pick indices for top and bottom layers
c
      if( igridv .eq. I_CART )then
        iztop = NZ
        izbot = 1
      else
        iztop = 1
        izbot = NZ
      endif

      do ie = 1, NELEM
        do j = 1,NBIN

          do ixyz = 1,NXYZ
            pc3(ixyz,j,ie) = pc3(ixyz,j,ie) * 
     $                       ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )
          enddo

          do ixy = 1,NXY
            ftoppart(ixy,j,ie) = ftoppart(ixy,j,ie) *
     $                           ( xmet2(ixy,iztop)*ymet2(ixy,iztop) )
            fbotpart(ixy,j,ie) = fbotpart(ixy,j,ie) *
     $                           ( xmet2(ixy,izbot)*ymet2(ixy,izbot) )
          enddo

        enddo
      enddo
c
c
c  Specify the values of <pc> assumed just above(below) the top(bottom)
c  of the model domain.
c  
      do ie = 1, NELEM
        do j = 1,NBIN
          do ixy=1,NXY
            pc_topbnd(ixy,j,ie) = pc2(ixy,NZ,j,ie)
            pc_botbnd(ixy,j,ie) = pc2(ixy,1,j,ie)
          enddo
        enddo
      enddo
c
c
c  Return to caller with particle concentrations initialized.
c
      return
      end
