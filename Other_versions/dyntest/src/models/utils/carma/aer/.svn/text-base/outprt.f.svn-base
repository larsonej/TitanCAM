      subroutine outprt
c
c
c  @(#) outprt.f  McKie  Oct-1995
c  This routine outputs information about the current timestep
c  to the output print file.
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
      dimension totn(NX,NY,NZ,NGROUP), rn(NX,NY,NZ,NGROUP), 
     $          rsig(NX,NY,NZ,NGROUP)

      character*(6) cwave
c
c
c  Define formats
c
    1 format('Timestep itime: ',i6,3x,'time: ',1pe14.6)
    2 format(/,'Total particle mass, total solute mass [g]: ',
     $   1p,2e14.6)
    3 format(/,'Particle size distributions at (lon,lat,iz) = ',
     $       '(', f8.3, ',', f8.3, ',', i4, ')' )
    4 format(i6,3x,i6,3x,1p,e12.2,3x,1p,e13.6)
    5 format(/,'Gas concentrations for ',a,' at (lon,lat) = ',
     $  '(', f8.3, ',', f8.3, ')', //,
     $  a3, 1x, 4(a11,4x), /)
    6 format(i3,1x,1p,3(e11.3,4x),0p,f11.3)
    7 format(/,'Water vapor, condensed, total [g]: ',1p,3e14.6)
    8 format(/,'Total particle number: ',1p,1e14.6)
    9 format(/,'Total CN mass, total {in,}volatile core mass [g]: ',
     $  1p,3e14.6)
   10 format(/, 2x, a, /)
   11 format(18x, a17, 18x, a19, /, 18x, 33('-'),3x,32('-'))
   12 format(a3, a12, 2(a12,a12,a12), /)
c  13 format(i3, 1p, e12.3, 6e12.3, /, (3x, 12x, 1p, 6e12.3) )
   13 format(i3, 1p, e12.3, 6e12.3 )
   14 format(3x,'Ielem',5x,'ibin',4x,'r (microns)',3x,'N or fraction')
   15 format(/,a, ' at (lat,lon) = (', f8.3, ',', f8.3, ')' ,/)
   16 format(/,'Albedo (',a,', ',a,') [%]: x across, y down',//,20i10)
   17 format(i4,20f10.1)
   18 format(/,'Optical Depth (',a,', ',a,'): x across, y down',//,
     $  20i10)
   19 format(i4,20(1pe10.2))
   20 format('Detailed radiative transfer output for lon, lat = ',
     $  2f8.3)
c
c---------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter outprt'
c
c
c  Report current timestep index and current simulation time
c
      call prtsep
      write(LUNOPRT,1) itime, time
c
c
c  Set vertical loop index to increment downwards
c
      if( igridv .eq. I_CART )then
        kb  = NZ
        ke  = 1
        idk = -1
      else
        kb  = 1
        ke  = NZ
        idk = 1
      endif
c
c
c  Radiative transfer output 
c
      if( do_rad )then
c       
c       
c  Detailed radiative transfer output requires the radiative transfer calculations
c  to be computed for each column 
c
        do ix = 1, NX
          do iy = 1, NY
 
            ixy = NX * ( iy - 1 ) + ix 
 
            write(LUNOPRT,20) xc(ix,iy,1), yc(ix,iy,1)
 
            call prerad
            call radtran
            call postrad
            call radout
 
          enddo
        enddo
c
c
c  2-D albedo maps
c
        write(LUNOPRT,'(a,f10.3)') 'u0 = ',u0

        write(LUNOPRT,16) 'integrated','toa',(ix,ix=1,NX)
        do iy = NY,1,-1
          write(LUNOPRT,17) iy,(100.*alb_toai(ix,iy),ix=1,NX)
        enddo

        write(LUNOPRT,16) 'integrated','tom',(ix,ix=1,NX)
        do iy = NY,1,-1
          write(LUNOPRT,17) iy,(100.*alb_tomi(ix,iy),ix=1,NX)
        enddo

        iwave = 9
        write(cwave,'(f3.1,'' um'')') wave(iwave)

        write(LUNOPRT,16) cwave,'toa',(ix,ix=1,NX)
        do iy = NY,1,-1
          write(LUNOPRT,17) iy,(100.*alb_toa(ix,iy,iwave),ix=1,NX)
        enddo

        write(LUNOPRT,18) cwave,'tom',(ix,ix=1,NX)
        do iy = NY,1,-1
          write(LUNOPRT,19) iy,(opd(ix,iy,iwave),ix=1,NX)
        enddo

      endif
c
c  
c  Report some gas concentrations and supersaturations at horiz space pt (ix,iy)
c
      ix = 1
      iy = 1
      ixy = NX * ( iy - 1 ) + ix

      do igas = 1,NGAS
 
        write(LUNOPRT,5) gasname(igas), xc(ix,iy,1), yc(ix,iy,1),
     $                   'iz','zc','gc [g/cm^3]','supsat','T [K]'
 
        do iz = kb,ke,idk
          xyzmet = xmet2(ixy,iz)*ymet2(ixy,iz)*zmet2(ixy,iz)
          write(LUNOPRT,6) iz,zc2(ixy,iz),gc2(ixy,iz,igas)/xyzmet,
     $                     supsatl2(ixy,iz,igas),t2(ixy,iz)
        enddo

      enddo
c
c
c  Report total particle mass and total solute mass
c
      rmcn = 0.
      rmdrop = 0.
      rmcore_inv = 0.
      rmcore_vol = 0.
      rntot = 0.

      do ixyz = 1,NXYZ
       do ie = 1,NELEM

        ig = igelem(ie)
        ip = ienconc(ig)

        do j = 1,NBIN

          if( itype(ie) .eq. I_INVOLATILE )then

            rmcn = rmcn + pc3(ixyz,j,ie)*rmass(j,ig) *
     $             dx3(ixyz)*dy3(ixyz)*dz3(ixyz)

          elseif( itype(ie) .eq. I_VOLATILE )then

            rmdrop = rmdrop + pc3(ixyz,j,ie)*rmass(j,ig) *
     $               dx3(ixyz)*dy3(ixyz)*dz3(ixyz)

          elseif( itype(ie) .eq. I_COREMASS )then

            if( itype(ip) .ne. I_INVOLATILE )then

              rmcore_vol = rmcore_vol + pc3(ixyz,j,ie) *
     $                     dx3(ixyz)*dy3(ixyz)*dz3(ixyz)
            else
              rmcore_inv = rmcore_inv + pc3(ixyz,j,ie) *
     $                     dx3(ixyz)*dy3(ixyz)*dz3(ixyz)
            endif

          endif

          if( itype(ie) .eq. I_INVOLATILE .or. 
     $        itype(ie) .eq. I_VOLATILE )then

            rntot = rntot + pc3(ixyz,j,ie) *
     $              dx3(ixyz)*dy3(ixyz)*dz3(ixyz)
          endif

        enddo
       enddo
      enddo

      rmcond = rmdrop - rmcore_vol

      igas = 1
      rmvap = 0.
      do ixyz = 1,NXYZ
        rmvap = rmvap + gc3(ixyz,igas) *
     $          dx3(ixyz)*dy3(ixyz)*dz3(ixyz)
      enddo
      write(LUNOPRT,7) rmvap, rmcond, rmvap+rmcond

      write(LUNOPRT,2) rmcn+rmdrop, rmcn+rmcore_vol
      write(LUNOPRT,9) rmcn, rmcore_inv, rmcore_vol
      write(LUNOPRT,8) rntot
c
c
c  Report parameters from log-normal fit to particle size distributions for
c  column at horizontal space point (ix,iy).
c
      call lognormal( totn, rn, rsig )

      ix = 1
      iy = 1

      write(LUNOPRT,15)
     $   'Log-normal fits to particle size distributions ', 
     $    xc(ix,iy,1), yc(ix,iy,1)

      write(LUNOPRT,11) 'CN', 'Drops'
      write(LUNOPRT,12) 'iz', 'zc', 'totn', 'rn', 'rsig',
     $                              'totn', 'rn', 'rsig'

      do iz = kb,ke,idk
        write(LUNOPRT,13) iz, zc(ix,iy,iz),
     $       ( totn(ix,iy,iz,ig), rn(ix,iy,iz,ig), rsig(ix,iy,iz,ig), 
     $         ig = 1,NGROUP )
      enddo
c
c
c  Report particle size distribution for a particular space point (ix,iy,iz)
c
c    <itype> = I_INVOLATILE or I_VOLATILE: number and mass concentrations
c    <itype> = I_COREMASS: core mass fraction and mass concentrations 
c    <itype> = I_CORE2MOM: ratio of core second moment to the bin second moment
c
      ix = 1
      iy = 1
      iz = 1

      ixyz = NXY * ( iz - 1 ) + NX * ( iy - 1 ) + ix

      write(LUNOPRT,3) xc(ix, iy, iz), yc(ix, iy, iz), iz

      xyzmet = xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)

      do ie = 1,NELEM

        write(LUNOPRT,10) elemname(ie)
        write(LUNOPRT,14)

        ig = igelem(ie)

        do j = 1,NBIN
          ip = ienconc(ig)

          if( ie .eq. ip )then
           outpc = pc3(ixyz,j,ie) / xyzmet
          else
           if( itype(ie) .eq. I_COREMASS .or.
     $         itype(ie) .eq. I_VOLCORE )then
            if( pc3(ixyz,j,ip) .eq. 0. )then
              outpc = 0.
            else
              outpc = pc3(ixyz,j,ie) / ( pc3(ixyz,j,ip)*rmass(j,ig))
            endif
           else
            if( pc3(ixyz,j,ip) .eq. 0. )then
              outpc = 0.
            else
              outpc = pc3(ixyz,j,ie) / (pc3(ixyz,j,ip)*rmass(j,ig)**2)
            endif
           endif
          endif

          write(LUNOPRT,4) ie, j, r(j,ig)*1.e4, outpc
 
        enddo
      enddo
c
c
c  Return to caller with timestep info output to print file
c
      return
      end
