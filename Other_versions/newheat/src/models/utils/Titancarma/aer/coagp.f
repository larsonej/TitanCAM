       subroutine coagp(ibin,ielem)
c
c
c  @(#) coagp.f  Jensen  Oct-1995
c  This routine calculates coagulation production terms <coagpe>.
c  See [Jacobson, et al., Atmos. Env., 28, 1327, 1994] for details
c  on the coagulation algorithm.
c
c  Argument list input:
c    ibin, ielem
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
c  Local declarations
c
c     dimension rmi(NXYZ)
c
c
c  Define formats
c
    1 format(1pe13.6,0p,6(i4),3x,1pe13.6,i4)
    2 format(2(i4),3x,3(1pe13.6,3x))
    3 format(a,1pe13.6,2x,3(i4),2x,'i:',1pe13.6,2x,
     $       2(i4),2x,'j:',1pe13.6,2x,2(i4),2x,
     $       'm, K:',2(2x,1pe13.6))
    4 format(i4,2(3x,a,3(i4)))
c
c
c  Announce entry to this routine
c
      if (DEBUG) write(LUNOPRT,'(/,a)') 'Enter coagp'
c
c  Definition of i,j,k,n used in comments: colision between i and j bins
c  yields particle between bins k and k+1.  Production in bin n is calculated.
c
c
c  Determine group that particles are produced in
c
      igrp = igelem(ielem)
c
c
c  Particle number production
c
c  Coagulation between particle in group <ig> bin <i> with particle in
c  group <jg> bin <j> results in particle with mass between bins k and k+1.
c  First, loop over group-bin quads <ig,i,jg,j> resulting in production in
c  bin <ibin> = k+1.  The set of quads <igup,jgup,iup,jup> is
c  defined in setupcoag.
c
cc    write(*,*) 'Coag Production:',ielem,ibin,igrp
cc    write(*,*) 'npairu',npairu(igrp,ibin)

      do iquad = 1, npairu(igrp,ibin)

       ig = igup(igrp,ibin,iquad)           ! source group
       jg = jgup(igrp,ibin,iquad)           ! source group
       i  = iup(igrp,ibin,iquad)            ! source bin
       j  = jup(igrp,ibin,iquad)            ! source bin

       iefrom = icoagelem(ielem,ig)         ! source element for <i> particle
       if( if_sec_mom(igrp) )then
         iefrom_cm = icoagelem_cm(ielem,ig) ! core mass moment source element
       endif
c
c  If <iefrom> = 0 then there is no contribution to production
c
cc     if( iefrom .eq. 0 ) 
cc   $  write(*,4) iquad,'no prod ',iefrom,i,ig,
cc   $                   'source j',ienconc(jg),j,jg

       if( iefrom .ne. 0 ) then

        je = ienconc(jg)                    ! source element for <j> particle
        if( if_sec_mom(igrp) )then
          je_cm = icoagelem_cm(ielem,jg)    ! core mass moment source element
        endif

cc      write(*,4) iquad,'source i',iefrom,i,ig,'source j',je,j,jg
cc      if( itype(ielem) .eq. I_CORE2MOM )
cc   $  write(*,4) iquad,'cmsrce i',iefrom_cm,i,ig,
cc   $                   'cmsrce j',je_cm,j,jg

c
c
c  If ielem is core mass type and <ig> is a CN type and <ig> is different
c  from <igrp>, then we must multiply production by mass
c  per particle (<rmi>) of element <icoagelem>.  (this is <rmass> for all source
c  elements except particle number concentration in a multicomponent CN group).
c
        do ixyz = 1, NXYZ
c
c
c  Bypass calculation if few source particles present
c
         if( pconmax(ixyz,ig)*pconmax(ixyz,jg) .gt. FEW_PC**2 )then

          rmi = 1.
          rmj = 1.
          i_pkern = 1
 
          if( itype(ielem) .eq. I_COREMASS .or.
     $        itype(ielem) .eq. I_GROWCORE .or.
     $        itype(ielem) .eq. I_VOLCORE )then          ! core mass

            i_pkern = 3                                 ! Use different kernel for core mass prod.

            if( ( itype(ienconc(ig)) .eq. I_INVOLATILE .or.
     $            itype(ienconc(ig)) .eq. I_VOLATILE )
     $          .and. ig .ne. igrp ) then               ! CN source and ig different from igrp

              if( ncore(ig) .eq. 0 )then                ! No cores in source group
                if( icomp(ielem) .eq. icomp(iefrom) )then
                 rmi = rmass(i,ig)
                else
                 rmi = 0.
                endif
              elseif( itype(iefrom) .eq. I_INVOLATILE .or.
     $                itype(iefrom) .eq. I_VOLATILE ) then
c
c  Source element is number concentration elem of mixed CN group
c
               totmass  = pc3(ixyz,i,iefrom) * rmass(i,ig)
               rmasscore = pc3(ixyz,i,icorelem(1,ig))
               do ic = 2,ncore(ig)
                 iecore = icorelem(ic,ig)
                 rmasscore = rmasscore + pc3(ixyz,i,iecore)
               enddo
               fracmass = 1. - rmasscore/totmass
               elemass  = fracmass * rmass(i,ig)
               rmi = elemass
              endif

            endif

          elseif( itype(ielem) .eq. I_CORE2MOM )then         ! core mass^2
                                                        ! Use different kernel for
            i_pkern = 5                                 ! core mass^2 production
            if( itype(ienconc(ig)) .eq. I_INVOLATILE ) 
     $        rmi = rmass(i,ig)
            if( itype(ienconc(jg)) .eq. I_INVOLATILE ) 
     $        rmj = rmass(j,jg)

          endif
c
c  For each spatial grid point, sum up coagulation production
c  contributions from each quad.
c 
          k = (ixyz - 1)/NXY + 1                        ! vertical level
c
          if( itype(ielem) .ne. I_CORE2MOM ) then

           coagpe(ixyz,ibin,ielem) = coagpe(ixyz,ibin,ielem) 
     $         + pc3(ixyz,i,iefrom)*pcl(ixyz,j,je)*rmi
     $         * pkernel(k,i,j,ig,jg,igrp,i_pkern)
 
          else
            coagpe(ixyz,ibin,ielem) = coagpe(ixyz,ibin,ielem) +
     $           ( pc3(ixyz,i,iefrom)*pcl(ixyz,j,je)*rmi**2 +
     $             pc3(ixyz,i,iefrom_cm)*rmi*
     $             pcl(ixyz,j,je_cm)*rmj ) *
     $           pkernel(k,i,j,ig,jg,igrp,i_pkern)

cc          if(ixyz.eq.21 .and. ielem.eq.7) 
cc   $     write(*,3) 'Coagpe1:',coagpe(ixyz,ibin,ielem),ixyz,
cc   $      ibin,ielem,pc3(ixyz,i,iefrom_cm),i,iefrom_cm,
cc   $      pcl(ixyz,j,je_cm),j,je_cm,rmj,
cc   $      pkernel(k,i,j,ig,jg,igrp,i_pkern)
          endif

c         if(coagpe(ixyz,ibin,ielem) .lt. 0.) 
c    $     write(*,3) 'Coagpe:',coagpe(ixyz,ibin,ielem),ixyz,
c    $      ibin,ielem,pc3(ixyz,i,iefrom),i,iefrom,
c    $      pcl(ixyz,j,je),j,je,rmi,
c    $      pkernel(k,i,j,ig,jg,igrp,i_pkern)
  
         endif    ! end of ( pconmax .gt. FEW_PC )

        end do    ! end of (ixyz = 1, NXYZ)

       endif      ! end of (iefrom .ne. 0)

      end do      ! end of (iquad = 1, npairu(igrp,ibin))
c
c
c  Next, loop over group-bin quads for production in bin <ibin> = k from
c  bin <i> due to collision between bins <i> and <j>.
c  Production will only occur if either k != <i> or igrp != <ig>
c 
      do iquad = 1, npairl(igrp,ibin)

       ig = iglow(igrp,ibin,iquad)
       jg = jglow(igrp,ibin,iquad)
       i  = ilow(igrp,ibin,iquad)
       j  = jlow(igrp,ibin,iquad)

       iefrom = icoagelem(ielem,ig)         ! source element for <i> particle
       if( if_sec_mom(igrp) )then
         iefrom_cm = icoagelem_cm(ielem,ig) ! core mass moment source element
c        if(iefrom_cm .eq. 0)
c    $     write(*,*) 'iefrom_cm = 0',ielem,ig,iefrom,igrp
       endif

       if( iefrom .ne. 0 ) then

        je = ienconc(jg)                    ! source element for <j> particle
        if( if_sec_mom(igrp) )then
          je_cm = icoagelem_cm(ielem,jg)    ! core mass moment source element
        endif

        do ixyz = 1, NXYZ
c
c
c  Bypass calculation if few particles present
c
         if( pconmax(ixyz,ig)*pconmax(ixyz,jg) .gt. FEW_PC**2 )then

          rmi = 1.
          rmj = 1.
          i_pkern = 2
 
          if( itype(ielem) .eq. I_COREMASS .or.
     $        itype(ielem) .eq. I_GROWCORE .or.
     $        itype(ielem) .eq. I_VOLCORE )then       ! core mass

            i_pkern = 4                         ! Use different kernel for core mass production

            if( ( itype(ienconc(ig)) .eq. I_INVOLATILE .or.
     $            itype(ienconc(ig)) .eq. I_VOLATILE )
c            if( itype(ienconc(ig)) .eq. I_INVOLATILE
     $          .and. ig .ne. igrp ) then               ! CN source and ig different from igrp

              if( ncore(ig) .eq. 0 )then                ! No cores in source group
                if( icomp(ielem) .eq. icomp(iefrom) )then
                 rmi = rmass(i,ig)
                else
                 rmi = 0.
                endif
              elseif( itype(iefrom) .eq. I_INVOLATILE .or.
     $                itype(iefrom) .eq. I_VOLATILE ) then
c
c  Source element is number concentration elem of mixed CN group
c
               totmass  = pc3(ixyz,i,iefrom) * rmass(i,ig)
               rmasscore = pc3(ixyz,i,icorelem(1,ig))
               do ic = 2,ncore(ig)
                 iecore = icorelem(ic,ig)
                 rmasscore = rmasscore + pc3(ixyz,i,iecore)
               enddo
               fracmass = 1. - rmasscore/totmass
               elemass  = fracmass * rmass(i,ig)
               rmi = elemass
              endif

            endif

          elseif( itype(ielem) .eq. I_CORE2MOM )then      ! core mass^2
                                                     ! Use different kernel for 
            i_pkern = 6                              ! core mass^2 production
            if( itype(ienconc(ig)) .eq. I_INVOLATILE ) 
     $        rmi = rmass(i,ig)
            if( itype(ienconc(jg)) .eq. I_INVOLATILE )
     $        rmj = rmass(j,jg)

          endif

          k = (ixyz - 1)/NXY + 1                        ! vertical level
c
          if( itype(ielem) .ne. I_CORE2MOM )then 

           coagpe(ixyz,ibin,ielem) = coagpe(ixyz,ibin,ielem) 
     $         + pc3(ixyz,i,iefrom)*pcl(ixyz,j,je)*rmi
     $         * pkernel(k,i,j,ig,jg,igrp,i_pkern)

          else
           coagpe(ixyz,ibin,ielem) = coagpe(ixyz,ibin,ielem) +
     $          ( pc3(ixyz,i,iefrom)*pcl(ixyz,j,je)*rmi**2 +
     $            pc3(ixyz,i,iefrom_cm)*rmi*
     $            pcl(ixyz,j,je_cm)*rmj ) *
     $          pkernel(k,i,j,ig,jg,igrp,i_pkern)

c           if(ixyz.eq.36 .and. ielem.eq.4) 
c    $     write(*,3) 'Coagpe2:',coagpe(ixyz,ibin,ielem),ixyz,
c    $      ibin,ielem,pc3(ixyz,i,iefrom_cm),i,iefrom_cm,
c    $      pcl(ixyz,j,je_cm),j,je_cm,rmi,
c    $      pkernel(k,i,j,ig,jg,igrp,i_pkern)
          endif

         endif    ! end of ( pconmax .gt. FEW_PC )

        end do    ! end of (ixyz = 1, NXYZ)

       endif      ! end of (iefrom .ne. 0)

      end do
c
c
c  Return to caller with coagulation production terms evaluated.
c 
      return
      end