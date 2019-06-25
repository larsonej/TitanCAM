      subroutine seinfeld86
c
c
c  @(#) seinfeld86.f  Jensen  Oct-1995
c  This routine sets up a coagulation calculation for
c  comparison with the self-preserving solution given by
c  Seinfeld [Atmos. Chem. and Phys. of Air Pollution, 1986, pgs 413-415]
c
c  Note: this routine should be called after the initialization
c  routines (at the end of init.f).
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
      include 'globals.h'
c
      dimension kbin(MAXNGRP,MAXNGRP,MAXNGRP,MAXNBINS,MAXNBINS)
c
c
c  Define formats
c
    1 format(/,3x,'Ielem',5x,'ibin',8x,'m',12x,'N(r)')
    2 format(i6,3x,i6,3x,e12.2,3x,e12.2)
    3 format(/,'Initial size distributions (seinfeld86):')
c
c
c  Set time step and nrun
c
      dtime = 1./80./pi
      ietime = 100
c
c
c  Initialize particle concentrations
c
      write(LUNOPRT,3)
      write(LUNOPRT,1)
      do ig = 1, ngroups
       do j = 1,nbins
        do ixyz=1,nxyz
         pc3(ixyz,j,ig) = 0.
        enddo
       enddo
      enddo
      do j = 1,nbins
       do ixyz=1,nxyz
        pc3(ixyz,j,1) = dm(j,1) / rmass(10,1) *
     $                exp(-rmass(j,1)/rmass(10,1))
        write(LUNOPRT,2) 1,j,rmass(j,1), c(1,1,1,j,1)
       enddo
      enddo
c
c
c  Change <icoag> mapping array
c
      do ig = 1, ngroups
       do jg = 1, ngroups
        icoag(ig,jg) = 0
       enddo
      enddo
      icoag(1,1) = 1
c
c
c  Set <ckernel> to a constant
c
      do i = 1,nbins
        do j = 1,nbins
          do k = 1, nz
            ckernel(k,i,j,1,1) = 16.*pi
          enddo
        enddo
      enddo
c
c
c  Recalculate <kbin> and <pkernel> arrays
c
      do igrp = 1, ngroups
        do ig = 1, ngroups
          do jg = 1, ngroups
            do i = 1, nbins
              do j = 1, nbins
 
                rmsum = rmass(i,ig) + rmass(j,jg)
 
                do ibin = 1, nbins-1
                  if( rmsum .ge. rmass(ibin,igrp) .and.
     $                rmsum .lt. rmass(ibin+1,igrp) ) then
                    kbin(igrp,ig,jg,i,j) = ibin
                  endif
                enddo
 
                if( rmsum .ge. rmass(nbins,igrp) )
     $                   kbin(igrp,ig,jg,i,j) = nbins
 
              enddo
            enddo
          enddo
        enddo
      enddo

      do igrp = 1, ngroups
        do ig = 1, ngroups
        do jg = 1, ngroups
 
          if( igrp .eq. icoag(ig,jg) ) then
 
            do i = 1, nbins
            do j = 1, nbins
 
              rmk = rmass(kbin(igrp,ig,jg,i,j),igrp)
              rmsum = rmass(i,ig) + rmass(j,jg)
 
              do k = 1, nz
 
                pkernl = ckernel(k,i,j,ig,jg)*
     $                   (rmrat(igrp)*rmk - rmsum) /
     $                   (rmrat(igrp)*rmk - rmk)
 
                pkernu = ckernel(k,i,j,ig,jg)*
     $                   (rmsum - rmk) /
     $                   (rmrat(igrp)*rmk - rmk)
 
                if( kbin(igrp,ig,jg,i,j) .eq. nbins )then
                  pkernl = ckernel(k,i,j,ig,jg)*
     $                     rmsum / rmass(nbins,igrp)
                  pkernu = 0.
                endif
 
                pkernel(k,i,j,ig,jg,igrp,1) = pkernu *
     $                                        rmass(i,ig)/rmsum
                pkernel(k,i,j,ig,jg,igrp,2) = pkernl *
     $                                        rmass(i,ig)/rmsum
                pkernel(k,i,j,ig,jg,igrp,3) = pkernu *
     $                                        rmk*rmrat(igrp)/rmsum
                pkernel(k,i,j,ig,jg,igrp,4) = pkernl *
     $                                        rmk/rmsum
 
              enddo
 
            enddo
            enddo
 
          endif
 
        enddo
        enddo
      enddo
c
c
c  Return to caller
c
      return
      end
