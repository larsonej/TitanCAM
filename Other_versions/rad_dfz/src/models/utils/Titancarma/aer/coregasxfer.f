      subroutine coregasxfer
c
c
c  @(#) coregasxfer.f 
c  This routine calculates... 
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Define Formats
c
    1 format(a,i3,3x,f5.2,' km  ',f8.6,' days  ratio: ',
     $ 1pe13.6,'  cloud: ',1pe13.6)
    2 format('growfrac: ',1pe11.4,'mgcr: ',1pe11.4,' mcld: ',1pe11.4,
     $ ' mcor: ',1pe11.4,' cldfrac: ',1pe11.4)
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter coregasxfer'
c
c
c-------------------------------------------------------------------------------
c
c
      do igroup = 1,NGROUP
       if( pconmax(ixyz,igroup) .gt. FEW_PC ) then

        iepart = ienconc(igroup)
        if( itype(iepart) .eq. I_VOLATILE ) then

         ipgas = igrowgas(iepart)  !cloud pc gas

c  Calculate production/loss of growcores from vapor-droplet exchange

         do ibin = 1,NBIN
          do ielemg = 2, nelemg(igroup)

           ielem = ielemg - 1 + iepart
           igas = igrowgas(ielem)  !growcore gas

           cloudmass = rmass(ibin,igroup)*pc3(ixyz,ibin,iepart)

           if( itype(ielem) .eq. I_GROWCORE ) then

            if( dfdt(ibin,ielem) .lt. 0.d0 ) then
           ! Loss [s-1]
              xle = -dfdt(ibin,ielem) / 
     $              ( pc3(ixyz,ibin,ielem) / cloudmass )

              xpe = 0.d0

            else
           ! Production [ g/cm3/s ]
              xle = 0.d0

              xpe = dfdt(ibin,ielem) * cloudmass
 
          ! check that prod doesn't exceed bin size
              if( dfdt(ibin,ielem) * dtime .gt. 1.d0 ) then
                write(70,*) 'growcore larger than cloud particle',ibin,
     $                  time/60.d0**2/24.d0
c               stop 1
              endif

            endif
            
            pc3(ixyz,ibin,ielem) = ( pc3(ixyz,ibin,ielem) + dtime*xpe )/
     $                             ( ONE + xle*dtime )            

              if( pc3(ixyz,ibin,ielem) .gt. 
     $             rmass(ibin,igroup)*pc3(ixyz,ibin,iepart) ) then
                 if( time .gt. 51092) then
                   write(*,1) 'growcore > cloud',ibin,zl3(ixyz)/1.d5,
     $                  time/60.d0**2/24.d0,pc3(ixyz,ibin,ielem)/
     $                  rmass(ibin,igroup)/pc3(ixyz,ibin,iepart),
     $                  pc3(ixyz,ibin,iepart)
                 endif
c                print *,pc3(ixyz,ibin,ielem),
c    $                rmass(ibin,igroup)*pc3(ixyz,ibin,iepart)
c                stop 1

                 pc3(ixyz,ibin,ielem) = 0.5 * ( pcl(ixyz,ibin,ielem)
     $              + pc3(ixyz,ibin,iepart)*rmass(ibin,igroup) )

                 xpe = (pc3(ixyz,ibin,ielem) - pcl(ixyz,ibin,ielem))/
     $                   dtime
              endif

          ! smallpc check ???

c  Calculate production/loss of gases from exchange

            exchange = xle * pc3(ixyz,ibin,ielem) - xpe

           !growcore gas
            gprod_ex(igas) = gprod_ex(igas)
     $                               + exchange

           !cloud particle gas
            gprod_ex(ipgas) = gprod_ex(ipgas)
     $                               - exchange

           endif !GROWCORE
          enddo !elements in group
         enddo !bins

   ! update Latent heat ???


c  Calculate relative fractions of growcore gases present in cloud

         do ibin = 1,NBIN 

          coretot = 0.d0
          cloudmass = rmass(ibin,igroup)*pc3(ixyz,ibin,iepart)

          do ielemg = 2,nelemg(igroup)
           ielem = ielemg - 1 + iepart   !use icorelem ?????

c        Total the mass of the non-growing cores
           if( itype(ielem) .eq. I_COREMASS .or.
     $         itype(ielem) .eq. I_VOLCORE ) then
          
             coretot = coretot + pc3(ixyz,ibin,ielem)

           endif
          enddo !elements in group

          xferfrac(ibin,iepart) = 1.d0

          do ielemg = 2,nelemg(igroup)
           ielem = ielemg - 1 + iepart
           
           if( itype(ielem) .eq. I_GROWCORE ) then

             xferfrac(ibin,ielem) = pc3(ixyz,ibin,ielem) /
     $                               (cloudmass - coretot)

             xferfrac(ibin,iepart) = xferfrac(ibin,iepart)
     $                                - xferfrac(ibin,ielem)


             if(time .gt. 51092) then
              write(*,*) 'FRAC ',ibin
              write(*,2) xferfrac(ibin,ielem),pc3(ixyz,ibin,ielem),
     $                   cloudmass,coretot,xferfrac(ibin,iepart)
               
             endif
           endif
          enddo !elements in group
         enddo !bins
         
        endif !VOLATILE

       endif !.gt. FEW_PC
      enddo !groups
c
c
c  Return to caller with <> evaluated.
c
      return
      end
