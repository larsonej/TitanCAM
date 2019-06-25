       subroutine setupbins
c
c
c  @(#) setupbins.f  Jensen  Oct-1995
c  This routine evaluates the derived mapping arrays and sets up
c  the particle size bins.
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
    1 format(a,':  ',12i6)
    2 format(a,':  ',i6)
    3 format(a,':  ',f12.2)
    4 format(a,':  ',12f12.2)
    5 format(/,'Particle grid structure (setupbins):')
    6 format(a,':  ',1p12e12.3)
    7 format(8(1pe13.6))
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupbins'
c
c
c  Determine which elements are particle number concentrations
c  <ienconc(igroup)> is the element corresponding to particle number 
c  concentration in group <igroup>
c
      igrp = 0
      do ielem = 1, NELEM
        if( itype(ielem) .eq. I_INVOLATILE .or. 
     $      itype(ielem) .eq. I_VOLATILE )then

          igrp = igrp + 1
          ienconc(igrp) = ielem
        endif
      enddo
      
      if( igrp .gt. NGROUP )then
        write(LUNOPRT,'(/,a)') 'bad itype array'
        stop 1
      endif

c
c  Determine which group each element belongs to
c  i.e., <igelem(ielem)> is the group to which element <ielem> belongs
c
      igrp = 0
      do ielem = 1, NELEM
        if( itype(ielem) .eq. I_INVOLATILE .or. 
     $      itype(ielem) .eq. I_VOLATILE )then
          igrp = igrp + 1
        endif
        igelem(ielem) = igrp
      enddo
c
c
c  Determine how many cores are in each group <ncore>.
c  The core elements in a group are given by <icorelem(1:ncore,igroup)>.
c
c  Also evaluate whether or not second moment is used <if_sec_mom> for each group.
c
      ielem = 0

      do igrp = 1, NGROUP

        ncore(igrp) = 0
        if_sec_mom(igrp) = .false.
        imomelem(igrp) = 0

        do j = 1, nelemg(igrp)

          ielem = ielem + 1

          if( itype(ielem) .eq. I_COREMASS .or.
     $        itype(ielem) .eq. I_GROWCORE .or.
     $        itype(ielem) .eq. I_VOLCORE )then

            ncore(igrp) = ncore(igrp) + 1
            icorelem(ncore(igrp),igrp) = ielem

          elseif( itype(ielem) .eq. I_CORE2MOM )then

            if_sec_mom(igrp) = .true.
            imomelem(igrp) = ielem

          endif

        enddo
      enddo
c
c
c  Particle mass densities (NXYZ*NBIN for each group) -- the user might want
c  to modify this (this code segment does not appear in setupaer subroutine
c  because <igelem> is not defined until this subroutine).
c
      do ig = 1,NGROUP
        ie = ienconc(ig)
        do ibin = 1,NBIN
          do ixyz = 1,NXYZ
            rhop3(ixyz,ibin,ig) = rhoelem(ie)
c  Set initial density of all hydrometeor groups to 1 such that nucleation
c  mapping arrays are calculated correctly.
c           if( itype(ie) .ne. I_INVOLATILE ) then
c             rhop3(ixyz,ibin,ig) = 1.
c           endif
          enddo
        enddo
      enddo
c
c
c  Set up the particle bins.
c  For each particle group, the mass of a particle in
c  bin j is <rmrat> times that in bin j-1
c
c    rmass(NBIN,NGROUP)     =  bin center mass [g]
c    r(NBIN,NGROUP)         =  bin mean (volume-weighted) radius [cm]
c    vol(NBIN,NGROUP)       =  bin center volume [cm^3]
c    dr(NBIN,NGROUP)        =  bin width in radius space [cm]
c    dv(NBIN,NGROUP)        =  bin width in volume space [cm^3]
c    dm(NBIN,NGROUP)        =  bin width in mass space [g]
c
      cpi = 4./3.*PI

      do igrp = 1, NGROUP

        rmassmin(igrp) = cpi*rhop3(1,1,igrp)*rmin(igrp)**3

       !For liquid methane cloud group adjust so that mass bins are
       !same as in ice methane cloud group (but radius will now be
       !different since density is different)
        if( groupname(igrp) .eq. 'Methane Liq Cloud group') then
          rmassmin(igrp) = rmassmin(igrp-1)
          rmin(igrp) = (rmassmin(igrp)/cpi/rhop3(1,1,igrp))**(ONE/3.)
          write(*,*) 'Adjusting mass and radius for liq methane clouds'
        endif

        vrfact = ( (3./2./PI/(rmrat(igrp)+1.))**(ONE/3.) )*
     $           ( rmrat(igrp)**(ONE/3.) - 1. )

        do j = 1, NBIN

          rmass(j,igrp)   = rmassmin(igrp) * rmrat(igrp)**(j-1)
          rmassup(j,igrp) = 2.*rmrat(igrp)/(rmrat(igrp)+1.)*
     $                      rmass(j,igrp)
          dm(j,igrp)      = 2.*(rmrat(igrp)-1.)/(rmrat(igrp)+1.)*
     $                      rmass(j,igrp)
          vol(j,igrp) = rmass(j,igrp) / rhop3(1,1,igrp)
          r(j,igrp)   = ( rmass(j,igrp)/rhop3(1,1,igrp)/cpi )**(ONE/3.)
          rup(j,igrp) = ( rmassup(j,igrp)/rhop3(1,1,igrp)/cpi )**
     $                                                          (ONE/3.)
          dr(j,igrp)  = vrfact*(rmass(j,igrp)/rhop3(1,1,igrp))**(ONE/3.)
          rlow(j,igrp) = rup(j,igrp) - dr(j,igrp)

c         write(*,7) rmass(j,igrp),rmassup(j,igrp),dm(j,igrp),
c    $               vol(j,igrp),r(j,igrp),rup(j,igrp),dr(j,igrp),
c    $               rlow(j,igrp)
        enddo
      enddo
c
c
c  Evaluate differences between values of <rmass> in different bins.
c
      do igrp = 1, NGROUP
       do jgrp = 1, NGROUP
        do i = 1, NBIN
         do j = 1, NBIN
           diffmass(i,igrp,j,jgrp) = rmass(i,igrp) - rmass(j,jgrp)
         enddo
        enddo
       enddo
      enddo
c
c
c  Report some initialization values
c
      if (do_print_setup) then
      write(LUNOPRT,5)
      write(LUNOPRT,2) 'NGROUP ',NGROUP
      write(LUNOPRT,2) 'NELEM  ',NELEM
      write(LUNOPRT,2) 'NBIN   ',NBIN
      write(LUNOPRT,6) 'Massmin',(rmassmin(i),i=1,NGROUP)
      write(LUNOPRT,4) 'Mrat   ',(rmrat(i),i=1,NGROUP)
      write(LUNOPRT,1) 'nelemg ',(nelemg(i),i=1,NGROUP)
      write(LUNOPRT,1) 'itype  ',(itype(i),i=1,NELEM)
      write(LUNOPRT,1) 'ienconc',(ienconc(i),i=1,NGROUP)
      write(LUNOPRT,1) 'igelem ',(igelem(i),i=1,NELEM)
      write(LUNOPRT,1) 'ncore  ',(ncore(i),i=1,NGROUP)
      endif
c
c
c  Return to caller with particle grid initialized
c
      return
      end
