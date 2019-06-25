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
    6 format(a,':  ',1p12e12.3)
    5 format(/,'Particle grid structure (setupbins):')
    7 format(/,'Fractal particles enabled (setupbins):')
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
        call endcarma
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
c  or not
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
c    rm(NBIN,NGROUP)		=  bin mobility radius [cm]
c	 rf(NBIN,NGROUP)		=  bin fractal radius [cm]
c    umon(NBIN,NGROUP)		=  number of monomers per aggregate
c	 porosity(NBIN,NGROUP)	=  porosity
c
      cpi = 4./3.*PI

      do igrp = 1, NGROUP

        rmassmin(igrp) = cpi*rhop3(1,1,igrp)*rmin(igrp)**3

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
c
c         define fractal particle radius, fractal dimension, and mobility radius
          if (ifractal(igrp) .eq. 1) then
            if (r(j,igrp) .le. rmon(igrp)) then
  
               umon(j,igrp) = 1.0
               df(j,igrp) = 3.0
  
               rf(j,igrp) = r(j,igrp)
               rm(j,igrp) = r(j,igrp)
               porosity(j,igrp) = 0.0

            elseif (r(j,igrp) .gt. rmon(igrp)) then
c  EJL - change df
               umon(j,igrp) = (r(j,igrp)/rmon(igrp))**3.
               df(j,igrp) = 3.0 - 1.0*exp(-umon(j,igrp)**(2./3.)/5000.0)
               rf(j,igrp) = r(j,igrp)**(3./df(j,igrp))* 
     $                  rmon(igrp)**(1.-3./df(j,igrp))
c  
c    calculating mobility radius for permeable aggregates us Vainshtein (2003)
c    formulation - from E. Wolf - EJL
c  
              tpor = 1.0 - (r(j,igrp)/rmon(igrp))**
     $           (3.*(1.-3./df(j,igrp)))
              porosity(j,igrp) = 1.0-(1.0 - tpor)*sqrt(df(j,igrp)/3.0)
              gamm = (1. - porosity(j,igrp))**(1./3.)
              happel = 2./(9.*(1.-porosity(j,igrp)))* 
     $               (3.-4.5*gamm+4.5*gamm**5.0-3.*gamm**6.)/
     $               (3.+2.*gamm**5.)
              perm = happel*rmon(igrp)**2.0
              brinkman = umon(j,igrp)**(1./df(j,igrp))*1./sqrt(happel)
              epsi = 1.0 - brinkman**(-1.)*tanh(brinkman)
              omega = 2.0/3.0*epsi/(2./3.+epsi/brinkman**2.0)
  
              rm(j,igrp) = rf(j,igrp) * omega
c  
            endif
          elseif (ifractal(igrp) .ne. 1) then
            umon(j,igrp) = 1.0
            df(j,igrp) = 3.0
            rf(j,igrp) = r(j,igrp)
            rm(j,igrp) = r(j,igrp)
            porosity(j,igrp) = 0.0

          endif
c
        enddo
      enddo
c
c
c  Evaluate differences between valuse of <rmass> in different bins.
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
c  Precalculate small values for use in smallconc()
c
      do ielem = 1, NELEM
        ig = igelem(ielem)
        ip = ienconc(ig)
        do ibin = 1, NBIN
          if( ielem .eq. ip )then
            small_val(ibin,ielem) = SMALL_PC
          elseif( itype(ielem) .eq. I_COREMASS .or.
     $            itype(ielem) .eq. I_VOLCORE )then
            small_val(ibin,ielem) =
     $         SMALL_PC*rmass(ibin,ig)*FIX_COREF
          elseif( itype(ielem) .eq. I_CORE2MOM )then
            small_val(ibin,ielem) =
     $         SMALL_PC*(rmass(ibin,ig)*FIX_COREF)**2
          endif
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
        if (ifractal(1) .eq. 1) then
          write(LUNOPRT,7)
        endif
      endif
c
c
c  Return to caller with particle grid initialized
c
c  Print bin structure
      open(unit=99, file='carma_bins.txt', status='unknown')

      do j = 1,NGROUP

        write(99,*) 'Group: ', groupname(j)
        write(99,*) 'bin    r(cm)     rf(cm)     rm(cm)',
     $              '      df     porosity      umon'
          do i = 1,NBIN
        write(99,'(i3,1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.4,1pe11.3)')
     $          i,r(i,j),rf(i,j),rm(i,j),df(i,j),porosity(i,j),umon(i,j)
          enddo
        write(99,*) ''
      enddo

      close(unit=99)

      return
      end
