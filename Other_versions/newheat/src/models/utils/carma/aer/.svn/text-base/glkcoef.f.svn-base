       subroutine glkcoef(k)
c
c
c  @(#) glkcoef.f  Jensen  Apr-1997
c  This routine calculates coefficients needed for horizontal
c  advection scheme (see htranglk.f).
c  Currently only valid for uniform grid and periodic B.C.s
c
c  Argument list input:
c    k
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables.
c
      include 'globaer.h'
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter glkcoef'
c
c
c  Calculate the coefficients used for horizontal transport.
c  x-direction first.
c
      idir = IDIRX

      do iy = 1,NY       ! Loop over latitudes

        do ix = 1,NX
          htrans(ix) = u(ix,iy,k)
        enddo

        do ix = 2,NX-1

          x1 = 2.*htrans(ix-1) + htrans(ix)
          x2 = -htrans(ix-1) + htrans(ix+1)
          x3 = 2.*htrans(ix+1) + htrans(ix)

          ca(1,ix,iy) = dx(ix,iy,k)/dtime - 0.5*x1
          cb(1,ix,iy) = 2.*( dx(ix,iy,k) + dx(ix,iy,k) )/dtime + 0.5*x2
          cd(1,ix,iy) = dx(ix,iy,k)/dtime + 0.5*x3
          ce(1,ix,iy) = dx(ix,iy,k)/dtime + 0.5*x1
          cf(1,ix,iy) = 2.*( dx(ix,iy,k) + dx(ix,iy,k) )/dtime - 0.5*x2
          cg(1,ix,iy) = dx(ix,iy,k)/dtime - 0.5*x3

        enddo
c
c  Periodic BCs
c
        ix = 1
        x1 = 2.*htrans(NX) + htrans(ix)
        x2 = -htrans(NX) + htrans(ix+1)
        x3 = 2.*htrans(ix+1) + htrans(ix)
        ca(1,ix,iy) = dx(ix,iy,k)/dtime - 0.5*x1
        cb(1,ix,iy) = 2.*( dx(ix,iy,k) + dx(ix,iy,k) ) / dtime + 0.5*x2
        cd(1,ix,iy) = dx(ix,iy,k)/dtime + 0.5*x3
        ce(1,ix,iy) = dx(ix,iy,k)/dtime + 0.5*x1
        cf(1,ix,iy) = 2.*( dx(ix,iy,k) + dx(ix,iy,k) ) / dtime - 0.5*x2
        cg(1,ix,iy) = dx(ix,iy,k)/dtime - 0.5*x3

        ix = NX
        x1 = 2.*htrans(ix-1) + htrans(ix)
        x2 = -htrans(ix-1) + htrans(1)
        x3 = 2.*htrans(1) + htrans(ix)
        ca(1,ix,iy) = dx(ix,iy,k)/dtime - 0.5*x1
        cb(1,ix,iy) = 2.*( dx(ix,iy,k) + dx(ix,iy,k) ) / dtime + 0.5*x2
        cd(1,ix,iy) = dx(ix,iy,k)/dtime + 0.5*x3
        ce(1,ix,iy) = dx(ix,iy,k)/dtime + 0.5*x1
        cf(1,ix,iy) = 2.*( dx(ix,iy,k) + dx(ix,iy,k) ) / dtime - 0.5*x2
        cg(1,ix,iy) = dx(ix,iy,k)/dtime - 0.5*x3

      enddo
c
c
c  y-direction next
c
      idir = IDIRY

      do ix = 1,NX       ! Loop over latitudes

        do iy = 1,NY
          htrans(iy) = v(ix,iy,k)
        enddo

        do iy = 2,NY-1

          x1 = 2.*htrans(iy-1) + htrans(iy)
          x2 = -htrans(iy-1) + htrans(iy+1)
          x3 = 2.*htrans(iy+1) + htrans(iy)

          ca(2,ix,iy) = dy(ix,iy,k)/dtime - 0.5*x1
          cb(2,ix,iy) = 2.*( dy(ix,iy,k) + dy(ix,iy,k) )/dtime + 0.5*x2
          cd(2,ix,iy) = dy(ix,iy,k)/dtime + 0.5*x3
          ce(2,ix,iy) = dy(ix,iy,k)/dtime + 0.5*x1
          cf(2,ix,iy) = 2.*( dy(ix,iy,k) + dy(ix,iy,k) )/dtime - 0.5*x2
          cg(2,ix,iy) = dy(ix,iy,k)/dtime - 0.5*x3

        enddo
c
c  Periodic BCs
c
        iy = 1
        x1 = 2.*htrans(NY) + htrans(iy)
        x2 = -htrans(NY) + htrans(iy+1)
        x3 = 2.*htrans(iy+1) + htrans(iy)
        ca(2,ix,iy) = dy(ix,iy,k)/dtime - 0.5*x1
        cb(2,ix,iy) = 2.*( dy(ix,iy,k) + dy(ix,iy,k) ) / dtime + 0.5*x2
        cd(2,ix,iy) = dy(ix,iy,k)/dtime + 0.5*x3
        ce(2,ix,iy) = dy(ix,iy,k)/dtime + 0.5*x1
        cf(2,ix,iy) = 2.*( dy(ix,iy,k) + dy(ix,iy,k) ) / dtime - 0.5*x2
        cg(2,ix,iy) = dy(ix,iy,k)/dtime - 0.5*x3

        iy = NY
        x1 = 2.*htrans(iy-1) + htrans(iy)
        x2 = -htrans(iy-1) + htrans(1)
        x3 = 2.*htrans(1) + htrans(iy)
        ca(2,ix,iy) = dy(ix,iy,k)/dtime - 0.5*x1
        cb(2,ix,iy) = 2.*( dy(ix,iy,k) + dy(ix,iy,k) ) / dtime + 0.5*x2
        cd(2,ix,iy) = dy(ix,iy,k)/dtime + 0.5*x3
        ce(2,ix,iy) = dy(ix,iy,k)/dtime + 0.5*x1
        cf(2,ix,iy) = 2.*( dy(ix,iy,k) + dy(ix,iy,k) ) / dtime - 0.5*x2
        cg(2,ix,iy) = dy(ix,iy,k)/dtime - 0.5*x3

      enddo
c
c
c  Return to caller with new coefficients for horizontal transport.
c
      return
      end
