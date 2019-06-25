       subroutine updatecoag
c
c
c  @(#) updatecoag.f  Barth  Jul-1999
c  (Based on setupcoag.f Jensen  Oct-1995)
c  This routine updates precomputed factors for coagulation after
c  ckernel has been changed (newckern.f must be called before this)
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter updatecoag'
c
c-----------------------------------------------------------------
c
c  Calculate variables needed in routine coagp.f
c
      do igrp = 1, NGROUP
        do ig = 1, NGROUP
        do jg = 1, NGROUP

          if( igrp .eq. icoag(ig,jg) ) then

            do i = 1, NBIN
            do j = 1, NBIN

              do k = 1, NZ
c
c  Zero some ckernel values to prevent tholin in bins 36-41
c  (Since ckernel is a function of groups, need to take care of 
c   preventing growcore-tholin interactions in {coagp}).
c
                if(ig.eq.1 .or. jg.eq.1) then
                 if(i.ge.36 .or. j.ge.36 ) then
                      ckernel(k,i,j,ig,jg) = ZERO
                 endif
                endif

                pkernel(k,i,j,ig,jg,igrp,1) = ckernel(k,i,j,ig,jg) *
     $                                     pkern0(k,i,j,ig,jg,igrp,1)
                pkernel(k,i,j,ig,jg,igrp,2) = ckernel(k,i,j,ig,jg) *
     $                                     pkern0(k,i,j,ig,jg,igrp,2)
                pkernel(k,i,j,ig,jg,igrp,3) = ckernel(k,i,j,ig,jg) *
     $                                     pkern0(k,i,j,ig,jg,igrp,3)
                pkernel(k,i,j,ig,jg,igrp,4) = ckernel(k,i,j,ig,jg) *
     $                                     pkern0(k,i,j,ig,jg,igrp,4)
                pkernel(k,i,j,ig,jg,igrp,5) = ckernel(k,i,j,ig,jg) *
     $                                     pkern0(k,i,j,ig,jg,igrp,5)
                pkernel(k,i,j,ig,jg,igrp,6) = ckernel(k,i,j,ig,jg) *
     $                                     pkern0(k,i,j,ig,jg,igrp,6)

              enddo

            enddo
            enddo

          endif

        enddo
        enddo
      enddo
c
c
c  Return to caller with coagulation mapping arrays defined
c
      return
      end
