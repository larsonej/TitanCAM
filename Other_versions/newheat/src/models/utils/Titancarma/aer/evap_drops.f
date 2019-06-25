      subroutine evap_drops(ibin,ig,ip)
c
c
c  @(#) evap_drops.f  Ackerman  Aug-2001
c  This routine calculates particle source terms <evappe> of droplets
c  evaporating within a particle group.
c
c  Distinct evaporation of cores has not been treated.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Define formats
c
    1 format(a,f20.15,2(i4),2x,f13.8,' days')
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter evap_drops'
c
c
c-------------------------------------------------------------------------------
c
c
c  The smallest bin cannot be a source to smaller bins in same group
c
      if( ibin .eq. 1 )then
        return
      endif
c
c
c  Evaluate evaporation source term <evappe> for all elements in group
c
      do isub = 1, nelemg(ig)
        ie = ip + isub - 1

c       if( itype(ie) .eq. I_VOLCORE ) then
c         evappe(ibin-1,ie) = evappe(ibin-1,ie) +
c    $                   diffmass(ibin,ig,ibin-1,ig) *
c    $                   pc3(ixyz,ibin,ip)*evaplg(ibin,ig)
c    $                 ( pc3(ixyz,ibin,ie)/ ( pc3(ixyz,ibin,ip)
c    $                    * rmass(ibin,ig) - pc3(ixyz,ibin,icm) )

        if( itype(ie) .eq. I_VOLCORE )
     $    vcmf = pc3(ixyz,ibin,ie)/(pc3(ixyz,ibin,ie)*rmass(ibin,ig))

        if( itype(ie) .eq. I_VOLCORE .and. vcmf .gt. FIX_COREF ) then
            evappe(ibin-1,ie) = evappe(ibin-1,ie) +
     $           pc3(ixyz,ibin,ip)*evaplg(ibin,ig)* 
     $                       diffmass(ibin,ig,ibin-1,ig) 
        else
          evappe(ibin-1,ie) = evappe(ibin-1,ie) +
     $                   pc3(ixyz,ibin,ie)*evaplg(ibin,ig)
        endif

          if(evappe(ibin-1,ie) .lt. ZERO) then
            print *,'ev13',evappe(ibin-1,ie),ie,ixyz,ibin,
     $              pc3(ixyz,ibin,ie),evaplg(ibin,ig)
              
c           stop 1
          endif

        if( itype(ie) .eq. I_GROWCORE ) then

          igas = igrowgas(ie)

c  <remove> is the total g/cm3 of growcore gas that was lost from the 
c  particle in ibin due to evaporation.  <togas> is the g/cm3 of growcore
c  gas that will be liberated to the vapor phase (from the ibin particle).
c  <remove> - <togas> is the g/cm3 of growcore gas which will be 
c  transfered with the particle to ibin-1

          remove = pc3(ixyz,ibin,ie) * evaplg(ibin,ig)

          dmdt_main = dmdt_gro(ibin-1,ig) - dmdte_gro(ibin-1,ie)
          if( dmdt_main .lt. ZERO ) then
           if( dmdte_gro(ibin-1,ie) .lt. ZERO ) then

            ! Both maingas and growcore gas evaporate from particle
            ! (keep at constant fraction for now)

c           if(ixyz.eq.iwa)
c    $        write(*,*) 'Evap both',ibin,evaplg(ibin,ig),itime,
c    $            dmdt_main,dmdt_gro(ibin,ig),dmdte_gro(ibin,ie)

cccc        remove = pc3(ixyz,ibin,ie) * evaplg(ibin,ig)
            togas = dmdte_gro(ibin-1,ie)/dmdt_gro(ibin-1,ig) *  
     $            pc3(ixyz,ibin,ip) * evaplg(ibin,ig) *
     $            diffmass(ibin,ig,ibin-1,ig) 

           else

            ! Maingas evaporates, but not growcore

c           if(ixyz.eq.iwa)
c    $        write(*,*) 'Evap main',ibin,evaplg(ibin,ig),itime

cccc        remove = ZERO
            togas = ZERO

           endif
          else
           ! Growcore evaporates, but not primary cloud volatile
           ! (may need to change this for multiple growcores)

c           if(ixyz.eq.iwa)
c    $        write(*,*) 'Evap growcore',ibin,evaplg(ibin,ig),itime

cccc       remove = pc3(ixyz,ibin,ie) * evaplg(ibin,ig)
           togas =  pc3(ixyz,ibin,ip) *diffmass(ibin,ig,ibin-1,ig) 
     $              * evaplg(ibin,ig) 

          endif

          if(togas .gt. remove) then
cc          if( max(togas,remove,togas-remove) .gt. 1d-14 )
cc   $       write(*,*) 'evap_drops: togas > remove',ibin-1,
cc   $                 remove,togas,remove-togas,ixyz,itime
            togas = remove
          endif

	  evappe(ibin-1,ie) = evappe(ibin-1,ie)  - togas

          if(evappe(ibin-1,ie) .lt. ZERO) then
            print *,'ev14',evappe(ibin-1,ie),ie,togas
            stop 1
          endif

          gprod_grow(ig,igas) = gprod_grow(ig,igas) + togas
	
          gprod_drop(ig,igas) = gprod_drop(ig,igas) + togas  !for debugging

ccc     else
        endif  !ielem = GROWCORE
 
      enddo

      return
      end
