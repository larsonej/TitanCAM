      subroutine smoluchowski
c
c
c  @(#) smoluchowski.f  Jensen  Oct-1995
c  This routine sets up a coagulation calculation for
c  comparison with the Smoluchowski analytic solution.
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
c  Set time step and nrun
c
      dtime = 1./80./pi
      ietime = 100
c
c
c  Initialize particle concentrations
c
      do ig = 1, ngroups
       do j = 1,nbins
        do ic3=1,nc3
         c3(ic3,j,ig) = 0.
        enddo
       enddo
      enddo
      c3(1,1,1) = 1.
c
c
c  Change <ncoag> mapping array
c
      do ig = 1, ngroups
       do jg = 1, ngroups
        ncoag(ig,jg) = 0
       enddo
      enddo
      ncoag(1,1) = 1
c
c
c  Set <kernel> to a constant
c
      do i = 1,nbins
        do j = 1,nbins
          ckernel(1,1,i,j) = 16.*pi
        enddo
      enddo
c
c
c  Recalculate <kbin> and <pkernel> arrays
c
      do ngrp = 1, ngroups
        do ig = 1, ngroups
          do jg = 1, ngroups
            do i = 1, nbins
              do j = 1, nbins
 
                rmres = rmass(ig,i) + rmass(jg,j)
                do nbin = 1, nbins-1
                  if( rmres .ge. rmass(ngrp,nbin) .and.
     $                rmres .lt. rmass(ngrp,nbin+1) ) then
                    kbin(ngrp,ig,jg,i,j) = nbin
                  endif
                enddo
                if( rmres .ge. rmass(ngrp,nbins) )
     $                   kbin(ngrp,ig,jg,i,j) = nbins
              enddo
            enddo
          enddo
        enddo
      enddo

      do ngrp = 1, ngroups
        do ig = 1, ngroups
        do jg = 1, ngroups
          do i = 1, nbins
          do j = 1, nbins
 
            rmk = rmass(ngrp,kbin(ngrp,ig,jg,i,j))
            rmsum = rmass(ig,i) + rmass(jg,j)
 
            pkernli = ckernel(ig,jg,i,j)*theta(i,j) *
     $          (rmrat(ngrp)*rmk - rmsum) /
     $          (rmrat(ngrp)*rmk - rmk)
 
            pkernui = ckernel(ig,jg,i,j)*theta(i,j) *
     $          (rmsum - rmk) /
     $          (rmrat(ngrp)*rmk - rmk)
 
            pkernlj = ckernel(jg,ig,j,i)*theta(i,j) *
     $          (rmrat(ngrp)*rmk - rmsum) /
     $          (rmrat(ngrp)*rmk - rmk)
 
            pkernuj = ckernel(jg,ig,j,i)*theta(i,j) *
     $          (rmsum - rmk) /
     $          (rmrat(ngrp)*rmk - rmk)
 
            if(kbin(ngrp,ig,jg,i,j) .eq. nbins) then
              pkernli = ckernel(ig,jg,i,j)*theta(i,j) *
     $          rmsum / rmass(ngrp,nbins)
              pkernui = 0.
              pkernlj = ckernel(jg,ig,j,i)*theta(i,j) *
     $          rmsum / rmass(ngrp,nbins)
              pkernuj = 0.
            endif

            pkernel(ngrp,ig,jg,i,j,1) = pkernui *
     $                                  rmass(ig,i)/rmsum
            pkernel(ngrp,ig,jg,i,j,2) = pkernuj *
     $                                  rmass(jg,j)/rmsum
            pkernel(ngrp,ig,jg,i,j,3) = pkernli *
     $                                  rmass(ig,i)/rmsum
            pkernel(ngrp,ig,jg,i,j,4) = pkernlj *
     $                                  rmass(jg,j)/rmsum
 
          enddo
          enddo
        enddo
        enddo
      enddo
c
c
c  Return to caller
c
      return
      end
