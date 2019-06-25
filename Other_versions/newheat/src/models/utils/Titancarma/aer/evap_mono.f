      subroutine evap_mono(ibin,ig,iavg,ieto,igto)
c
c
c  @(#) evap_mono.f  Ackerman  Aug-2001
c  This routine calculates particle source terms <evappe> due to total
c  evaporation from bin <ibin> group <ig> into a monodisperse
c  distribution.
c
c  Distinct evaporation of cores has not been treated.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c
      logical conserve_mass
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter evap_mono'
c
c
c-------------------------------------------------------------------------------
c
c
c  Define option to conserve mass or number when a choice must be made
c  during monodisperse total evaporation beyond CN grid -- should be done in setupaer()
c
      conserve_mass = .true.
c
c  Set automatic flag for total evaporation used in gasexchange()
c
      totevap(ibin,ig) = .true.
c
c
c  Possibly put all of core mass into largest, smallest, or
c  smallest nucelated CN bin 
c
      if( too_big .or. too_small .or. nuc_small )then

        if( too_big )then
          if( igto.eq.1 ) then
           jbin = NCCNBIN
          else
           jbin = NBIN
          endif
        elseif( too_small )then
          jbin = 1
        else
          jbin = inucmin(igto)
        endif

        if( conserve_mass )then
          factor = coreavg/rmass(jbin,igto)
        else
          factor = ONE
        endif
c
c  First the CN number concentration element
c 
        evappe(jbin,ieto) = evappe(jbin,ieto) + factor*evdrop

c       if( ieto.eq.1 ) write(*,*) itime,' evap_mono:',jbin,ixyz,
c    $       evappe(jbin,ieto),factor,evdrop

        if(evappe(jbin,ieto) .lt. 0.) then
          print *,'ev1',evappe(jbin,ieto),ieto
          stop 1
        endif

c
c  Now the CN cores
c
        do ic = 1, ncore(ig)
          iecore = icorelem(ic,ig)
          ie2cn  = ievp2elem(iecore)

          if( ie2cn .ne. 0 ) then
            if( ie2cn .ne. ieto ) then   !! Don't double count core that
                                         !! evaps to #conc element
              evappe(jbin,ie2cn) = evappe(jbin,ie2cn) +
     $           factor*evcore(ic)*rmass(jbin,igto)

              if(evappe(jbin,ie2cn) .lt. 0.) then
                print *,'ev2',evappe(jbin,ie2cn),ie2cn
                stop 1
              endif
 
              if( if_sec_mom(ig) .and. itype(iecore).eq.I_COREMASS) then
                iemom = imomelem(ig)
                ie2c2m = ievp2elem(iemom) 
                evappe(jbin,ie2c2m) = evappe(jbin,ie2c2m) +
     $            pc3(ixyz,jbin,iemom)*evaplg(jbin,ig)

c               write(*,*) 'jbin',jbin,igto,coreavg,rmass(NBIN,igto),
c    $               too_big,too_small,nuc_small,inucmin(igto),
c    $               evappe(jbin,ie2c2m)
              endif

            endif !(ie2cn .ne. ieto)

          else

c    Return growcore component to the gas phase for a total
c    evaporating particle

            if( itype(iecore) .eq. I_GROWCORE ) then
              igas = igrowgas(iecore)
              gprod_grow(ig,igas) = gprod_grow(ig,igas) 
     $               + pc3(ixyz,ibin,iecore) * evaplg(ibin,ig)
     $               + pc3(ixyz,ibin,iecore+1) * evaplg(ibin,ig)

c            if(ixyz.eq.iwa)
c    $        write(*,*) 'gprod_grow(mono): ',gprod_grow(ig,igas)

              gprod_mono(ig,igas) = gprod_mono(ig,igas)
     $               + pc3(ixyz,ibin,iecore) * evaplg(ibin,ig)
     $               + pc3(ixyz,ibin,iecore+1) * evaplg(ibin,ig)

            endif     
          endif
 
        enddo

      else
c
c
c  Partition core mass between two CN bins, conserving total core mass
c  and number.  The number will be subdivided into bins <iavg> and <iavg>-1.
c
       if( iavg .le. 1 .or. iavg .gt. NBIN )then
         print*,' stop in evap_mono: bad iavg = ', iavg
         call endcarma
       endif

       fracmass = ( rmass(iavg,igto) - coreavg ) /
     $    diffmass(iavg,igto,iavg-1,igto)
       fracmass = max( ZERO, min( ONE, fracmass ) )
c
c  First the CN number concentration element
c 
       evappe(iavg-1,ieto) = evappe(iavg-1,ieto) + evdrop*fracmass

       if(evappe(iavg-1,ieto) .lt. 0.) then
         print *,'ev5',evappe(iavg-1,ieto),ieto
         stop 1
       endif

       evappe(iavg,ieto) = evappe(iavg,ieto) + evdrop*( ONE - fracmass )

       if(evappe(iavg,ieto) .lt. 0.) then
          print *,'ev6',evappe(iavg,ieto),ieto,evdrop,fracmass
          stop 1
       endif

c
c  Now the cores
c
       do ic = 1, ncore(ig)
         iecore = icorelem(ic,ig)
         ie2cn  = ievp2elem(iecore)

         if( ie2cn .ne. 0 ) then
          if( ie2cn .ne. ieto ) then
           evappe(iavg-1,ie2cn) = evappe(iavg-1,ie2cn) +
     $        rmass(iavg-1,igto)*evcore(ic)*fracmass

           if(evappe(iavg-1,ie2cn) .lt. 0.) then
             print *,'ev7',evappe(iavg-1,ie2cn),ie2cn
             stop 1
           endif

           evappe(iavg,ie2cn) = evappe(iavg,ie2cn) +
     $        rmass(iavg,igto)*evcore(ic)*( ONE - fracmass )

           if(evappe(iavg,ie2cn) .lt. 0.) then
             print *,'ev8',evappe(iavg,ie2cn),ie2cn
             stop 1
           endif

           if( if_sec_mom(ig) .and. itype(iecore).eq.I_COREMASS ) then
             iemom = imomelem(ig)
             ie2c2m = ievp2elem(iemom) 

             evappe(iavg-1,ie2c2m) = evappe(iavg-1,ie2c2m) +
     $         pc3(ixyz,iavg-1,iemom)*evaplg(iavg-1,ig)*fracmass

             evappe(iavg,ie2c2m) = evappe(iavg,ie2c2m) +
     $         pc3(ixyz,iavg,iemom)*evaplg(iavg,ig)*(ONE-fracmass)

c            write(*,*) 'iavg',evappe(iavg-1,ie2c2m),
c    $                   evappe(iavg,ie2c2m)
           endif 

          endif !(ie2cn .ne. ieto)
         else

c    Return growcore component to the gas phase for a total
c    evaporating particle

          if( itype(iecore) .eq. I_GROWCORE ) then
            igas = igrowgas(iecore)
            gprod_grow(ig,igas) = gprod_grow(ig,igas) 
     $               + pc3(ixyz,ibin,iecore) * evaplg(ibin,ig)
     $               + pc3(ixyz,ibin,iecore+1) * evaplg(ibin,ig)

c          if(ixyz.eq.iwa)
c    $      write(*,*) 'gprod_grow(mono): ',gprod_grow(ig,igas)

            gprod_mono(ig,igas) = gprod_mono(ig,igas)
     $               + pc3(ixyz,ibin,iecore) * evaplg(ibin,ig)
     $               + pc3(ixyz,ibin,iecore+1) * evaplg(ibin,ig)

          endif     
         endif 
       enddo
 
      endif

      return
      end
