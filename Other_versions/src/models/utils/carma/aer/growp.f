      subroutine growp(ibin,ielem)
c
c
c  @(#) growp.f  Ackerman  Dec-1995
c  This routine calculates particle source terms due to growth <growpe>
c  for one particle size bin at one spatial grid point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c
c  Argument list input:
c    ibin, ielem
c
c  Argument list output:
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter growp'

      xyzmet = xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)
      rhoa_cgs = rhoa3(ixyz) /  xyzmet
c
c
c  Define group & particle # concentration indices for current element
c
      igroup = igelem(ielem)      ! target particle group 
      iepart = ienconc(igroup)	  ! target particle number concentration element
c
c
c  Calculate production terms due to condensational growth <growpe>
c  only if group to which element belongs grows.
c
      if( igrowgas(iepart) .ne. 0 .and. ibin .ne. 1 )then
c
c
c  Bypass calculation if few droplets are present 
c
         if( pconmax(ixyz,igroup) .gt. FEW_PC )then

          growpe(ibin,ielem) = pc3(ixyz,ibin-1,ielem)
     $                            * growlg(ibin-1,igroup) 
         endif
      endif
c
c
c  Return to caller with growth production terms evaluated.
c
      return
      end
