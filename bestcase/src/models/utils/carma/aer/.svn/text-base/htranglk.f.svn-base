       subroutine htranglk(idir,iloop)
c
c
c  @(#) htranglk.f  Jensen  Apr-1997
c  This routine calculates horizontal advection rates using
c  Galerkin method with Chapeau functions.
c  Currently only valid for uniform grid and periodic B.C.s
c
c  Argument list input:
c    idir, iloop
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
c  Local declarations
c
      dimension cal(NXORNY), cbl(NXORNY), cdl(NXORNY),
     $  cel(NXORNY), cfl(NXORNY), cgl(NXORNY),
     $  den(NXORNY), tpl(NXORNY), spl(NXORNY),
     $  vpl(NXORNY), dl(NXORNY), el(NXORNY),
     $  fl(NXORNY)

c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter htranglk'
c
c
c  <idir> = IDIRX: x-direction; <idir>=IDIRY: y-direction
c
      if( idir .eq. IDIRX ) then
        nadv = NX
        do ix = 1,NX
          cal(ix) = ca(1,ix,iloop)
          cbl(ix) = cb(1,ix,iloop)
          cdl(ix) = cd(1,ix,iloop)
          cel(ix) = ce(1,ix,iloop)
          cfl(ix) = cf(1,ix,iloop)
          cgl(ix) = cg(1,ix,iloop)
        enddo
      end if
      if( idir .eq. IDIRY ) then
        nadv = NY
        do iy = 1,NY
          cal(iy) = ca(2,iloop,iy)
          cbl(iy) = cb(2,iloop,iy)
          cdl(iy) = cd(2,iloop,iy)
          cel(iy) = ce(2,iloop,iy)
          cfl(iy) = cf(2,iloop,iy)
          cgl(iy) = cg(2,iloop,iy)
        enddo
      end if
c
c
c  Source terms
c
      do i = 2,nadv-1
        dl(i) = cel(i)*chor(i-1) + cfl(i)*chor(i) +
     $           cgl(i)*chor(i+1)
      enddo

      i = 1
      dl(i) = cel(i)*chor(nadv) + cfl(i)*chor(i) +
     $         cgl(i)*chor(i+1)

      i = nadv
      dl(i) = cel(i)*chor(i-1) + cfl(i)*chor(i) +
     $         cgl(i)*chor(1)
c
c
c  Now use tridiagonal solver
c
      el(1) = dl(1) / cbl(1)
      vpl(1) = 0.
      vpl(nadv) = 0.
      tpl(1) = 0.
      spl(1) = -cal(1) / cbl(1)
      fl(1) = -cdl(1) / cbl(1)

      do i = 2,nadv

        den(i) = 1./(cbl(i) + cal(i)*fl(i-1))
        fl(i) = -cdl(i)*den(i)
        spl(i) = -cal(i)*spl(i-1)*den(i)
        eltemp = dl(i) - cal(i)*el(i-1)
          if( abs(eltemp) .le. 1.e-11*abs(dl(i)) ) eltemp = 0.
        el(i) = eltemp*den(i)

      enddo

      tpl(nadv) = 1.

      do i = nadv-1, 1, -1
        tpl(i) = tpl(i+1)*fl(i) + spl(i)
        vpl(i) = vpl(i+1)*fl(i) + el(i)
      enddo

      divn = den(nadv) * ( cdl(nadv)*tpl(1) + cal(nadv)*spl(nadv-1) )
     $       + 1.
      chor(nadv) = ( el(nadv) - vpl(1)*cdl(nadv)*den(nadv) ) / divn

      do i = nadv, 2, -1
        xnum1 = fl(i-1)*chor(i) + spl(i-1)*chor(nadv)
        xnum2 = el(i-1) + xnum1
        chor(i-1) = xnum2
          if( abs(xnum2) .le. 1.e-11*abs(xnum1) ) chor(i-1) = 0.
      enddo
c
c
c  Remove negative concentrations from the grid
c  using a downstream borrowing technique  
c
      qprob1  = 0.
      qprob2  = 0.
      qsum    = 0.
      qcor    = 0.
      do i = 1, nadv-1
        qsum = qsum + chor(i) * ( dhor(i) + dhor(i+1) )/2.
      enddo
      qsum = qsum + chor(nadv) * ( dhor(nadv) + dhor(1) )/2.

      do i = 1, nadv-1

        if( i .eq. nadv-1 ) then
          dxx = (dhor(i)+dhor(i+1))/(dhor(1)+dhor(i+1))
        else if( i .ne. nadv-1 ) then
          dxx = (dhor(i)+dhor(i+1))/(dhor(i+1)+dhor(i+2))
        endif

        if( chor(i) .lt. 0. )
     $      chor(i+1) = chor(i+1) + chor(i)*dxx
        if( chor(i) .lt. 0. ) chor(i) = 0.

      enddo

      i = nadv
      ibot = 1
      dxx = (dhor(i)+dhor(1))/(dhor(1)+dhor(ibot+1))
      if( chor(i) .lt. 0. )
     $    chor(1) = chor(1) + chor(i)*dxx
      if( chor(i) .lt. 0. ) chor(i) = 0.

      do i = 1, nadv-1

        if( i .eq. nadv-1 ) then
          dxx = (dhor(i)+dhor(i+1))/(dhor(1)+dhor(i+1))
        else if( i .ne. nadv-1 ) then
          dxx = (dhor(i)+dhor(i+1))/(dhor(i+1)+dhor(i+2))
        endif

        if( chor(i) .lt. 0. )
     $      chor(i+1) = chor(i+1) + chor(i)*dxx
        if( chor(i) .lt. 0. ) chor(i) = 0.

      enddo

      dxx = 0.5*( dhor(1) + dhor(nadv) )
      qcor = 0.
      if( chor(nadv) .lt. 0. ) qcor = chor(nadv)*dxx
      if( chor(nadv) .lt. 0. ) chor(nadv) = 0.

c...Rescale to conserve mass

      if( ( (qsum .eq. 0.) .and. (qcor .eq. 0.) ) .or.
     $         (qsum-qcor) .eq. 0. ) then
        qcor1 = 0.
      else if( ( (qsum .ne. 0.) .or. (qcor .ne. 0.) ) .and.
     $         (qsum-qcor) .ne. 0. ) then
        qcor1 = qsum / ( qsum - qcor )
      endif

      do i = 1,nadv
        chor(i) = chor(i) * qcor1
      enddo
c
c
c  Return to caller with new concentrations after horizontal transport.
c
      return
      end
