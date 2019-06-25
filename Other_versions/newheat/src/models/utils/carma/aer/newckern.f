       subroutine newckern 
c 
c 
c  @(#) newckern.f  Barth Jul-1999
c 
c  This routine evaluates the coagulation kernels, ckernel(k,j1,j2,i1,i2)
c  [cm^-3 s^-1]. Indices correspond to vertical level <k>, aerosol groups
c  <j1,j2> and bins <i1,i2> of colliding particles.
c
c  This routine requires that vertical profiles of temperature <T>,
c  air density <rhoa>, and viscosity <rmu> are defined.
c  (i.e., initatm.f must be called before this)
c  The vertical profile with ix = iy = 1 is used.
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
    1 format(4(1pe13.6)) 
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter newckern'
c
c-------------------------------------------------------------------------------
c
c  electron charge [esu]
            e = -4.8d-10
c
c  Use constant kernel if <icoagop> = 0
c
      if( icoagop .eq. 0 )then
        do j1 = 1, NGROUP
         do j2 = 1, NGROUP
          do i1 = 1, NBIN
           do i2 = 1, NBIN
            do k = 1, NZ
             ckernel(k,i1,i2,j1,j2) = ck0
            enddo
           enddo
          enddo
         enddo
        enddo
        return   ! Return to caller with coagulation kernels evaluated.
      endif

c  Use column ix = iy = 1
      ix = 1
      iy = 1
c
c  Set up charge to radius ratio array <prrat>
c
!      call charge
!	  open(unit=99, file='charge.txt',status='unknown')
!	  write(99,*) prrat
!	  write(99,*) ''
!	  write(99,*) icoag(j1,j2)
!	  close(unit=99)
c
c  Loop over groups
c
      j1 = 1
      j2 = 1

      igrp = icoag(j1,j2)
      if( igrp .ne. 0 )then
   
c
c  First particle
c
        do i1 = 1, NBIN
c
c  Second particle
c
         do i2 = 1, NBIN
c
c  Loop over vertical layers
c
          do k = 1, NZ
c
c  First calculate thermal coagulation kernel
c  
c  <cstick> is the probability that two particles that collide
c  through thermal coagulation will stick to each other.
c
c  Sticking coefficient for charged particles below 300km from
c  Toon et al. (Icarus 95, 24-53, 1992)

!            if(zl(1,1,k) .le. 300.d5) then
             cstick = exp(-r(i1,j1)*rf(i2,j2)*prrat(k)*prrat(k)*e*e /
     $                          (BK*t(1,1,k)*(rf(i1,j1)+r(i2,j2))))
!            else
!             cstick = 1.
!            endif
			

c
c
c   <cbr> is thermal (brownian) coagulation coefficient
c
            ckernel(k,i1,i2,j1,j2) = cbr_term0(k,i1,i2,j1,j2) /
     $                      ( cbr_term1(k,i1,i2,j1,j2)
     $                          + cbr_term2(k,i1,i2,j1,j2)/cstick )
c
c  Update precomputed factors for coagulation <pkernel> with new ckernels
c
c
            pkernel(k,i1,i2,j1,j2,igrp,1) = ckernel(k,i1,i2,j1,j2) *
     $                                     pkern0(k,i1,i2,j1,j2,igrp,1)
            pkernel(k,i1,i2,j1,j2,igrp,2) = ckernel(k,i1,i2,j1,j2) *
     $                                     pkern0(k,i1,i2,j1,j2,igrp,2)
            pkernel(k,i1,i2,j1,j2,igrp,3) = ckernel(k,i1,i2,j1,j2) *
     $                                     pkern0(k,i1,i2,j1,j2,igrp,3)
            pkernel(k,i1,i2,j1,j2,igrp,4) = ckernel(k,i1,i2,j1,j2) *
     $                                     pkern0(k,i1,i2,j1,j2,igrp,4)
            pkernel(k,i1,i2,j1,j2,igrp,5) = ckernel(k,i1,i2,j1,j2) *
     $                                     pkern0(k,i1,i2,j1,j2,igrp,5)
            pkernel(k,i1,i2,j1,j2,igrp,6) = ckernel(k,i1,i2,j1,j2) *
     $                                     pkern0(k,i1,i2,j1,j2,igrp,6)
c
c
           enddo     ! vertical level
          enddo    ! second particle bin
         enddo    ! first particle bin
        endif     ! icoag ne 0 

c
c
c  return to caller with coagulation kernels evaluated.
c
      return
      end
