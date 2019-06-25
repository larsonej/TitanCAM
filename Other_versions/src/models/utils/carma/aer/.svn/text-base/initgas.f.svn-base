       subroutine initgas
c
c
c  @(#) initgas.f  Ackerman  Dec-1995
c  This routine initializes the atmospheric profiles of all gases.
c
c    gc       Gas concentration at layer mid-point [g/cm^3]
c
c  Presently the only vertical coordinate is altitude, zl [cm].
c
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
c   Define formats
c
    1 format(/,'Gas concentrations for ',a,'(initgas)',//,
     $  a3, 1x, 4(a11,4x), /) 
    2 format(i3,1x,1p,3(e11.3,4x),0p,f11.3)

c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter initgas'
c
c-------------------------------------------------------------------------------
c
c
c   Calculate vapor pressures
c
      do ixyz = 1,NXYZ
        call vaporp
      enddo
c
c
c   Parameters for water vapor profile
c
      igas = 1
c
c
c   <rh_init> is initial relative humidity [%]
c
      rh_init = 95.0d0
c
c
c   Gas constant for water vapor
c
      rvap = RGAS/gwtmol(igas)
c
c
c   Loop over all spatial dimensions and gases
c
      do igas = 1,NGAS
       do ix = 1,NX
        do iy = 1,NY
         do iz = 1,NZ

           if( igas .eq. 1 )then
c
c
c   Water vapor concentration from relative humidity and vapor pressure
c
              gc(ix,iy,iz,igas) = rh_init/100.*pvapl(ix,iy,iz,igas)
     $                        / ( rvap*t(ix,iy,iz) ) 

           else
             write(LUNOPRT,'(/,a)') 'invalid <igas> in initgas.f'
             call endcarma
           endif

         enddo
        enddo
       enddo
      enddo
c
c
c  Specify fluxes at top and bottom of model [g/cm^2/s]
c
      do igas = 1,NGAS
        do ixy = 1,NXY
          ftopgas(ixy,igas) = 0.
          fbotgas(ixy,igas) = 0.
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

      do igas = 1, NGAS

        do ixyz = 1,NXYZ
          gc3(ixyz,igas) = gc3(ixyz,igas) * 
     $                     ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )
        enddo

        do ixy = 1,NXY
          ftopgas(ixy,igas) = ftopgas(ixy,igas) *
     $                        ( xmet2(ixy,iztop)*ymet2(ixy,iztop) )
          fbotgas(ixy,igas) = fbotgas(ixy,igas) *
     $                        ( xmet2(ixy,izbot)*ymet2(ixy,izbot) )
        enddo
      enddo
c
c
c  Initialize <supsati> and <supsatl>
c
      do ixyz = 1,NXYZ
        call supersat
      enddo
c
c
c  Print gas concentrations at 1 horizontal grid point
c
      ix = 1
      iy = 1
      ixy = NX * ( iy - 1 ) + ix 
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

      do igas = 1,NGAS

        write(LUNOPRT,1) gasname(igas),
     $                   'iz','zc','gc [g/cm^3]','supsat','T [K]'

        do iz = kb,ke,idk
          xyzmet = xmet2(ixy,iz)*ymet2(ixy,iz)*zmet2(ixy,iz)
          write(LUNOPRT,2) iz, zc2(ixy,iz), gc2(ixy,iz,igas)/xyzmet,
     $       supsatl2(ixy,iz,igas), t2(ixy,iz)
        enddo

      enddo
c
c
c  Specify the values of <gc> assumed just above(below) the top(bottom)
c  of the model domain.  
c
      do igas = 1,NGAS
        do ixy=1,NXY
          gc_topbnd(ixy,igas) = gc2(ixy,NZ,igas)
          gc_botbnd(ixy,igas) = gc2(ixy,1,igas)
        enddo
      enddo
c
c  Return to caller with gas concentrations initialized.
c
      return
      end
