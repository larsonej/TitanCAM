       subroutine coagl
c
c
c  @(#) coagl.f  Jensen  Oct-1995
c  This routine calculates coagulation loss rates <coaglgg>.
c  See [Jacobson, et al., Atmos. Env., 28, 1327, 1994] for details
c  on the coagulation algorithm.
c
c  The loss rates for all particle elements in a particle group are equal.
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
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter coagl'
c
c
c  Loop over particle groups for which coagulation loss is being
c  calculated.
c
      do ig = 1,NGROUP
c
c  Element corresponding to particle number concentration
c
        ie = ienconc(ig)
c
c  Loop over particle groups that particle in group ig might
c  collide with.
c
        do jg = 1,NGROUP

          je = ienconc(jg)
c
c  Particle resulting from coagulation between groups <ig> and <jg> goes
c  to group <igrp>
c
          igrp = icoag(ig,jg)
c
c
c  Resulting particle is in same group as particle under consideration --
c  partial loss (muliplies <volx>).
c
          if( igrp .eq. ig )then
c
c  Loop over entire spatial grid
c
           do ixyz = 1, NXYZ
c
c  Bypass calculation if few source particles present
c
            if( pconmax(ixyz,jg)*pconmax(ixyz,ig) .gt. FEW_PC**2 )then

              do i = 1, NBIN-1
               do j = 1, NBIN
                k = (ixyz - 1)/(NX*NY) + 1 
                coaglg(ixyz,i,ig) = coaglg(ixyz,i,ig) + 
     $                    ckernel(k,i,j,ig,jg) * pcl(ixyz,j,je) *
     $                    volx(igrp,ig,jg,i,j) 
              enddo
             enddo
            endif
           enddo
c
c
c  Resulting particle is in a different group -- complete loss (no <volx>).
c
          else if( igrp .ne. ig .and. igrp .ne. 0 )then
c
c  Bypass calculation if few particles present
c
           do ixyz = 1, NXYZ
            if( pconmax(ixyz,jg)*pconmax(ixyz,ig) .gt. FEW_PC**2 )then
             do i = 1, NBIN
              do j = 1, NBIN
                k = (ixyz - 1)/(NX*NY) + 1
                coaglg(ixyz,i,ig) = coaglg(ixyz,i,ig) +
     $                    ckernel(k,i,j,ig,jg)*pcl(ixyz,j,je)
              enddo
             enddo
            endif
           enddo
          endif
        enddo  
      enddo
c
c
c  Boundary condition: Particles from bin <NBIN> are only lost by
c  coagulating into other elements. (This is taken care of by <NBIN>-1
c  limit above)
c
c
c  Return to caller with particle loss rates due to coagulation evaluated.
c
      return
      end
