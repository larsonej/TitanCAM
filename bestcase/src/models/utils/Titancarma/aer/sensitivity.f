      subroutine sensitivity
c
c
c  @(#) sensitivity.f  Barth  Jun-2002
c  This routine changes some parameters (defined in the <init...> and 
c  <setup...> routines) for performing sensitivity tests on the model
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter sensitivity'
c
c
c  Changes to eddy diffusion profile (default is Toon et al. 1992)
c  --originally defined in <initatm>
c
c   Eddy diffusion from Lara et al. 
 
       do k=1,NZP1
         do ix=1,NX
           do iy=1,NY

              dkz_toublanc = 2000.
              zl_km = zl3(k)/1.d5
              dkz_strobel = 4.d8 * exp( 1.76d5/1.6 *(zl_km-985.)/
     $                                (RTITAN+985.)/(RTITAN+zl_km) )
 
              dkz(ix,iy,k) = (dkz_toublanc**2 * dkz_strobel)**(1./3.)
 
           enddo
         enddo
       enddo
c
c    Change eddy diffusion by a constant factor

       do k=1,NZP1
         do ix=1,NX
           do iy=1,NY

             dkz(ix,iy,k) = 5.*dkz(ix,iy,k)
 
           enddo
         enddo
       enddo
c
c
c  Changes to contact parameter (default is ct=0.986 for scrit=1.15)
c  --originally defined in <setupnuc>
c
       ct(1)=0.9736     !scrit = 1.3
       ct(1)=0.966      !scrit = 1.4
       ct(1)=0.9589     !scrit = 1.5
c
c
c  Changes to gas flux into top of model (default is from Yung et al. 1984)
c  --originally defined in <initgas>
c
      do ixy = 1,NXY
        ftopgas(ixy,igas) = 2.d0 * ftopgas(ixy,igas)
      enddo
c
c
c  Changes to tholin profile
c
c     pc(ix,iy,iz,...) = 
c
c
c  Write any changes to output file...
c
       write(LUNOPRT,*) 'Eddy Diffusion Profile'
       do k=1,NZ
          write(LUNOPRT,*) zl3(k)/1.d5, dkz(k)
       enddo

       write(LUNOPRT,*) 'Ethane Flux'
       write(LUNOPRT,*) ftopgas(1,1)

       write(LUNOPRT,*) 'Contact Parameter'
       write(LUNOPRT,*) ct(1)
c
c  Return to caller with changes to ...
c
      return
      end
