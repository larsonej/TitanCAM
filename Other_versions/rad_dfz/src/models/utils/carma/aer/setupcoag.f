       subroutine setupcoag
c
c
c  @(#) setupcoag.f  Jensen  Oct-1995
c  This routine sets up mapping arrays and precomputed
c  factors for coagulation.
c
c  This routine requires that <ckernel> has been defined.
c  (setupckern.f must be called before this)
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
      dimension kbin(NGROUP,NGROUP,NGROUP,NBIN,NBIN)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupcoag'
c
c-----------------------------------------------------------------
c
c
c  Initialize <icoagelem> with zeros
c
      do ielem = 1,NELEM
        do ig = 1,NGROUP
          icoagelem(ielem,ig) = 0
          icoagelem_cm(ielem,ig) = 0
        enddo
      enddo
c
c  For each element <ielem> and each group <ig>, determine which element in <ig>
c  contributes to production  in <ielem>: <icoagelem(ielem,ig)>.
c  If no elements in <ig> are transfered into element <ielem> during coagulation,
c  then set <icoagelem(ielem,ig)> to 0.
c
      do ielem = 1,NELEM

        isolto = isolelem(ielem)           ! target solute type
        icompto = icomp(ielem)            ! target element compound
        igto = igelem(ielem)               ! target group

        do ig = 1, NGROUP                 ! source group

          iepart = ienconc(ig)            ! source particle number concentration element
c
c
c  If <ielem> is particle number concentration type or <ig> is pure group, then
c  the source element is the particle number concentration element of group <ig>.
c
          if( ( itype(ielem) .eq. I_INVOLATILE .or. 
     $          itype(ielem) .eq. I_VOLATILE ) .or.
     $          ncore(ig) .eq. 0 ) then
            icoagelem(ielem,ig) = iepart
c
c
c  If <ielem> is a core mass element, use element compound names to determine
c  which source element matches target core mass element.
c
          else if( itype(ielem) .eq. I_COREMASS )then
            icompfrom = icomp(iepart)       ! source element compound
            if( icompfrom .eq. icompto )then
              icoagelem(ielem,ig) = iepart
            else
              do ic = 1,ncore(ig)
                iecore = icorelem(ic,ig)       ! absolute element number of core
                icompfrom = icomp(iecore)    ! source element compound
                if( icompfrom .eq. icompto )
     $            icoagelem(ielem,ig) = iecore
              enddo
            endif
c
c
c  If <ielem> is a core second moment element, the source element is either
c  the number concentration element or the second moment concentration.
c
          else if( itype(ielem) .eq. I_CORE2MOM )then
            icompfrom = icomp(iepart)       ! source element compound
            if( icompfrom .eq. icompto )then
              icoagelem(ielem,ig) = iepart
            else
              do ic = 1,ncore(ig)
                iecore = icorelem(ic,ig)
                icompfrom = icomp(iecore)
                if( icompfrom .eq. icompto )
     $            icoagelem(ielem,ig) = imomelem(ig)
              enddo
            endif
c
c
c  For core second moment elements, we need additional pairs of source
c  elements c  to account for core moment production due to products
c  of source particle core mass.
c
c  If <ig> is a pure group, use particle number concentration element,
c  otherwise, use core mass element that matches the core second
c  moment composition
c
            icompfrom = icomp(iepart)       ! source element compound
            if( icompfrom .eq. icompto )then
              icoagelem_cm(ielem,ig) = iepart
            else
              do ic = 1,ncore(ig)
                iecore = icorelem(ic,ig)       ! absolute element number of core
                icompfrom = icomp(iecore)    ! source element compound
                if( icompfrom .eq. icompto )
     $            icoagelem_cm(ielem,ig) = iecore
              enddo
            endif

          endif
c
c
c  If <ielem> is a core mass type and <ig> is a pure CN group and the
c  solutes don't match, then set <icoagelem> to zero to make sure no
c  coag production occurs.
c
          if( itype(ielem) .eq. I_COREMASS .and.
     $        itype(ienconc(ig)).eq. I_INVOLATILE
     $        .and. ncore(ig) .eq. 0 ) then
            isolfrom = isolelem(ienconc(ig))
            if( isolfrom .ne. isolto ) then
              icoagelem(ielem,ig) = 0
            endif
          endif

        enddo          ! end of (ig = 1, NGROUP)

      enddo            ! end of (ielem = 1,NELEM)
c
c
c  Calculate lower bin <kbin> which coagulated particle goes into
c  and make sure it is less than <NBIN>+1
c
c  Colliding particles come from group <ig>, bin <i> and group <jg>, bin <j>
c  Resulting particle lands in group <igrp>, between <ibin> and <ibin> + 1
c
      do igrp = 1, NGROUP
        do ig = 1, NGROUP
          do jg = 1, NGROUP
            do i = 1, NBIN
              do j = 1, NBIN

                rmsum = rmass(i,ig) + rmass(j,jg)

                do ibin = 1, NBIN-1
                  if( rmsum .ge. rmass(ibin,igrp) .and.
     $                rmsum .lt. rmass(ibin+1,igrp) ) then
                    kbin(igrp,ig,jg,i,j) = ibin
                  endif
                enddo

                if( rmsum .ge. rmass(NBIN,igrp) )
     $                   kbin(igrp,ig,jg,i,j) = NBIN

              enddo
            enddo
          enddo
        enddo
      enddo
c
c
c  Calculate partial loss fraction
c
c  This fraction is needed because when a particle in bin <i> collides
c  with a particle in bin <j> resulting in a particle whose mass falls
c  between <i> and <i>+1, only partial loss occurs from bin <i>.
c
c  Since different particle groups have different radius grids, this
c  fraction is a function of the colliding groups and the resulting group.
c
      do igrp = 1, NGROUP
        do ig = 1, NGROUP
          do jg = 1, NGROUP

            if( igrp .eq. icoag(ig,jg) ) then

              do i = 1, NBIN
                do j = 1,NBIN

                  volx(igrp,ig,jg,i,j) = 1.
 
                  if(kbin(igrp,ig,jg,i,j).eq.i) then
 
                    rmkbin = rmass(kbin(igrp,ig,jg,i,j),igrp)
                    volx(igrp,ig,jg,i,j) = 1. -
     $                  (rmrat(igrp)*rmkbin-rmass(i,ig)-rmass(j,jg))
     $                  /(rmrat(igrp)*rmkbin-rmkbin)*
     $                  rmass(i,ig)/(rmass(i,ig) + rmass(j,jg))

                  endif

                enddo
              enddo
            endif
          enddo
        enddo
      enddo
c
c
c  Calculate mapping functions that specify sets of quadruples
c  (group pairs and bin pairs) that contribute to production
c  in each bin. Mass transfer from <ig,i> to <igrp,ibin> occurs due to
c  collisions between particles in <ig,i> and particles in <jg,j>.
c  2 sets of quadruples must be generated:
c     low: k = ibin and (k != i or ig != igrp)  and  icoag(ig,jg) = igrp
c      up: k+1 = ibin        and  icoag(ig,jg) = igrp
c
c  npair#(igrp,ibin) is the number of pairs in each set (# = l,u)
c  i#, j#, ig#, and jg# are the bin pairs and group pairs in each
c  set (# = low, up)
c
      do igrp = 1, NGROUP
        do ibin = 1, NBIN

          npairl(igrp,ibin) = 0
          npairu(igrp,ibin) = 0

          do ig = 1, NGROUP
          do jg = 1, NGROUP
            do i = 1, NBIN
            do j = 1, NBIN

              kb = kbin(igrp,ig,jg,i,j)
              ncg = icoag(ig,jg)

              if( kb+1.eq.ibin .and. ncg.eq.igrp ) then
                npairu(igrp,ibin) = npairu(igrp,ibin) + 1
                iup(igrp,ibin,npairu(igrp,ibin)) = i
                jup(igrp,ibin,npairu(igrp,ibin)) = j
                igup(igrp,ibin,npairu(igrp,ibin)) = ig
                jgup(igrp,ibin,npairu(igrp,ibin)) = jg
              endif

              if( kb.eq.ibin .and. ncg.eq.igrp .and.
     $          (i.ne.ibin .or. ig.ne.igrp) ) then
                npairl(igrp,ibin) = npairl(igrp,ibin) + 1
                ilow(igrp,ibin,npairl(igrp,ibin)) = i
                jlow(igrp,ibin,npairl(igrp,ibin)) = j
                iglow(igrp,ibin,npairl(igrp,ibin)) = ig
                jglow(igrp,ibin,npairl(igrp,ibin)) = jg
              endif

            enddo
            enddo
          enddo
          enddo
        enddo
      enddo
c
c
c  Do some extra debugging reports  (normally commented)
c
c      write(LUNOPRT,*) ' '
c      write(LUNOPRT,*) 'Coagulation group mapping:'
c      do ig = 1, NGROUP
c        do jg = 1, NGROUP
c          print*, 'ig jg icoag = ', ig, jg, icoag(ig,jg)
c        enddo
c      enddo
c      write(LUNOPRT,*) ' '
c      write(LUNOPRT,*) 'Coagulation element mapping:'
c      do ielem = 1, NELEM
c        do ig = 1, NGROUP
c          print*, 'ielem ig icoagelem icomp(ielem) = ',
c     $     ielem, ig, icoagelem(ielem,ig), icomp(ielem)
c        enddo
c      enddo
c      write(LUNOPRT,*) ' '
c      write(LUNOPRT,*) 'Coagulation bin mapping arrays'
c      do igrp = 1, NGROUP
c        do ibin = 1,3
c          write(LUNOPRT,*) 'igrp, ibin = ',igrp, ibin
c          do ipair = 1,npairl(igrp,ibin)
c            write(LUNOPRT,*) 'low:np,ig,jg,i,j ',
c     $          ipair,iglow(igrp,ibin,ipair),
c     $      jglow(igrp,ibin,ipair), ilow(igrp,ibin,ipair),
c     $            jlow(igrp,ibin,ipair)
c          enddo
c          do ipair = 1,npairu(igrp,ibin)
c            write(LUNOPRT,*) 'up:np,ig,jg,i,j ',
c     $          ipair,igup(igrp,ibin,ipair),
c     $      jgup(igrp,ibin,ipair), iup(igrp,ibin,ipair),
c     $           jup(igrp,ibin,ipair)
c          enddo
c        enddo
c      enddo
c      call endcarma
c
c
c  Calculate variables needed in routine coagp.f
c

! -EJL using Wolf's kluge.
!      open(unit=99, file='carma_coag_pkernel.txt', status='unknown')
!      write(99,*) 'pkernel?'
	  
      do igrp = 1, NGROUP
        do ig = 1, NGROUP
        do jg = 1, NGROUP
		
          if( igrp .eq. icoag(ig,jg) ) then

            do i = 1, NBIN
            do j = 1, NBIN

              rmk = rmass(kbin(igrp,ig,jg,i,j),igrp)
              rmsum = rmass(i,ig) + rmass(j,jg)
			  
              do k = 1, NZ

                pkernl =(rmrat(igrp)*rmk - rmsum) /
     $                   (rmrat(igrp)*rmk - rmk)

                pkernu = (rmsum - rmk) /
     $                   (rmrat(igrp)*rmk - rmk)

!                write(99,*) ckernel, rmrat, pkernl

                if( kbin(igrp,ig,jg,i,j) .eq. NBIN )then
                  pkernl = rmsum / rmass(NBIN,igrp)
                  pkernu = 0.
                endif
  
                pkern0(k,i,j,ig,jg,igrp,1) = pkernu *
     $                                        rmass(i,ig)/rmsum
                pkern0(k,i,j,ig,jg,igrp,2) = pkernl *
     $                                        rmass(i,ig)/rmsum
                pkern0(k,i,j,ig,jg,igrp,3) = pkernu *
     $                                        rmk*rmrat(igrp)/rmsum
                pkern0(k,i,j,ig,jg,igrp,4) = pkernl *
     $                                        rmk/rmsum
                pkern0(k,i,j,ig,jg,igrp,5) = pkernu *
     $                                 ( rmk*rmrat(igrp)/rmsum )**2
                pkern0(k,i,j,ig,jg,igrp,6) = pkernl *
     $                                        ( rmk/rmsum )**2


c  Set constant <pkernel>s for time-independent <ckernel>s
                pkernel(k,i,j,ig,jg,igrp,1) = ckernel(k,i,j,ig,jg) *
     $                                      pkern0(k,i,j,ig,jg,igrp,1)
                pkernel(k,i,j,ig,jg,igrp,2) = ckernel(k,i,j,ig,jg) *
     $                                      pkern0(k,i,j,ig,jg,igrp,2)
                pkernel(k,i,j,ig,jg,igrp,3) = ckernel(k,i,j,ig,jg) *
     $                                      pkern0(k,i,j,ig,jg,igrp,3)
                pkernel(k,i,j,ig,jg,igrp,4) = ckernel(k,i,j,ig,jg) *
     $                                      pkern0(k,i,j,ig,jg,igrp,4)
                pkernel(k,i,j,ig,jg,igrp,5) = ckernel(k,i,j,ig,jg) *
     $                                      pkern0(k,i,j,ig,jg,igrp,5)
                pkernel(k,i,j,ig,jg,igrp,6) = ckernel(k,i,j,ig,jg) *
     $                                      pkern0(k,i,j,ig,jg,igrp,6)

              enddo

            enddo
            enddo

          endif

        enddo
        enddo
      enddo

!      close(unit=99)  !EJL
c
c
c  Return to caller with coagulation mapping arrays defined
c
      return
      end
