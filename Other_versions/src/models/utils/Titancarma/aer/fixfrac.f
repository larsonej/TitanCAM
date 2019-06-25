      subroutine fixfrac
c
c
c  @(#) fixfrac.f  Barth  Jul-2003
c  (Aug-2007 compatible with volcores)
c  Check fractions of all volatile cores (VOLCORE,GROWCORE) and
c  adjust if necessary
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
c  Define formats
c
    1 format(a,i4,'  pc(gcore): ',1pe11.4,'  df: ',e11.4,
     $       '  binmass: ',e11.4,'  f: ',e11.4,'  fmax: ',
     $       e11.4,'  volmass: ',e11.4,'  cmf: ',e11.4,
     $       '  pc(cloud): ',e11.4)

c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter fixfrac'
   
      do igrp = 1,NGROUP

       if(is_grp_mixed_phase(igrp) .or. is_grp_mixed_comp(igrp) ) then

         iepart = ienconc(igrp)
         imaingas = igrowgas(iepart)
         do ibin = 1,NBIN
          if( pc3(ixyz,ibin,iepart) .gt. FEW_PC ) then
c
c  Calculate core mass fractions
c    cmf - involatile core
c   gcmf - growcore (ice+liquid)
c   vcmf - ice core corresponding to primary
c
           gcoretot = 0.d0
           do ic = 1,ncore(igrp)
            iecore = icorelem(ic,igrp)
            if(itype(iecore) .eq. I_COREMASS ) then
             coreavg = pc3(ixyz,ibin,iecore) / pc3(ixyz,ibin,iepart)
            elseif(itype(iecore) .eq. I_VOLCORE) then
             if( itype(iecore-1) .eq. I_GROWCORE ) then
              gcoretot = gcoretot + pc3(ixyz,ibin,iecore)
              if( t3(ixyz) .le. Tfreez(igrowgas(iecore-1)) ) then
                iegc = iecore
              else
                iegc = iecore - 1
              endif
             else
              vcmf = pc3(ixyz,ibin,iecore) / pc3(ixyz,ibin,iepart) 
     $                          / rmass(ibin,igrp)
              ievc = iecore
             endif
            elseif(itype(iecore) .eq. I_GROWCORE) then
              gcoretot = gcoretot + pc3(ixyz,ibin,iecore)
              if( .not. is_grp_mixed_phase(igrp)) iegc = iecore
            endif
           enddo !cores
           coreavg = min( rmass(ibin,igrp),coreavg )
           cmf(ibin,igrp) = coreavg / rmass(ibin,igrp)
           one_cmf = ONE - cmf(ibin,igrp)

           gcoreavg = gcoretot / pc3(ixyz,ibin,iepart)
           gcmf = gcoreavg / rmass(ibin,igrp)

           binmass = pc3(ixyz,ibin,iepart) * rmass(ibin,igrp)
c
c  Check to see if growcore fraction > 1
c
           if( gcmf .gt. one_cmf ) then

c  Calculate mass <df*binmass> to remove from growcore and 
c  adjust gases accordingly (calculations simplified to one
c  line but detailed below)
ccc         fmax = ONE - cmf(ibin,igrp)
ccc         f = pc3(ixyz,ibin,iesub) / binmass
ccc         df = fmax - f

            dfbinmass = one_cmf*binmass - gcoretot

c           write(*,*) itime,' fixfrac/gc:',ixyz,ibin,vcmf,
c    $        cmf(ibin,igrp),gcmf,dfbinmass,pc3(ixyz,ibin,iepart),
c    $        igrp

c           if(itime.eq.14694.and.ixyz.eq.14.and.igrp.eq.3)
c    $       write(*,*) 'cores in fixfrac:',ibin,pc3(ixyz,ibin,ievc),
c    $             pc3(ixyz,ibin,iegc),pc3(ixyz,ibin,iegc+1),
c    $             rmass(ibin,igrp),ievc,iegc
c
c   Update growcore to reflect new f and distribute gases
c   accordingly
c
            pc3(ixyz,ibin,iegc) = pc3(ixyz,ibin,iegc)
     $                                  + dfbinmass

            gc3(ixyz,igas) = gc3(ixyz,igas) - dfbinmass
            gc3(ixyz,imaingas) = gc3(ixyz,imaingas) + dfbinmass

           endif
c
c  Check to see if too much mass in ice core assoc w/ primary
c  condensate on cloud particle
c
           if( vcmf .gt. one_cmf - gcmf ) then
c           write(*,*) itime,' fixfrac:',ixyz,ibin,vcmf,cmf(ibin,igrp),
c    $         gcmf,ievc,pc3(ixyz,ibin,ievc),(one_cmf-gcmf)*binmass

             pc3(ixyz,ibin,ievc) = 
     $                      (one_cmf - gcmf)*binmass  

c          if(itime.gt.220000 .and. ixyz.eq.11)
c    $      write(*,*) itime,' fixfrac:',ixyz,ibin,vcmf,one_cmf,gcmf,
c    $        pc3(ixyz,ibin,ievc)/binmass,pc3(ixyz,ibin,iepart)

           endif

          endif !small_pc
         enddo !bins
       endif !mixed group
      enddo !group
c
c
c  Return to caller with growcore fractions fixed
c
      return
      end
