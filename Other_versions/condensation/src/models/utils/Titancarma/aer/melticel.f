      subroutine melticel
c
c
c  @(#) melticel.f  Jensen  Jan-2000
c  This routine evaluates particle loss rates due to nucleation <rnuclg>:
c  Ice crystal melting (to mixed phase particles) only.
c  
c  The loss rates for all particle elements in a particle group are equal.
c
c  Calculations are done at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c
c  Argument list input:
c
c  Argument list output:
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter melticel'
c
c
c  Loop over particle groups.
c
      do igroup = 1,NGROUP

       iepart = ienconc( igroup )            ! particle number density element
c
c
c  Calculate nucleation loss rates.
c
c  <ienucto> is index of target nucleation element;
c  <ignucto> is index of target nucleation group.
c
       do inuc = 1,nnuc2elem(iepart)

        ienucto = inuc2elem(inuc,iepart)
        if( ienucto .ne. 0 )then
          ignucto = igelem( ienucto )
        else
          ignucto = 0
        endif
c
c
c  Only compute nucleation rate for ice crystal melting
c
        if( inucproc(iepart,ienucto) .eq. I_ICEMELT ) then
c
c
c  Loop over particle bins.  Loop from largest to smallest for 
c  evaluation of index of smallest bin nucleated during time step <inucstep>.
c
         do ibin =NBIN,1,-1
c
c  Bypass calculation if few particles are present 
c
cc        if( pconmax(ixyz,igroup) .gt. FEW_PC )then
          if( pc3(ixyz,ibin,iepart) .gt. FEW_PC )then
c
c
c  Temporary simple kludge: Set <rnuclg> to 1.e2 if T > 0C
c  
           if( t3(ixyz) .gt. T0(igrowgas(iepart)) ) then

            rnuclg(ibin,igroup,ignucto) = 1.e2

cc          if( ixyz.eq.9 .and. rnuclg(ibin,igroup,ignucto) .gt. 0. ) 
cc   $         write(*,*) 'Melting:',ixyz,ibin,igroup,ignucto,
cc   $           rnuclg(ibin,igroup,ignucto),itime,time/8.64d4,
cc   $           pc3(ixyz,ibin,iepart)

            inucstep(igroup) = ibin

           endif

          endif   ! pconmax(ixyz,igroup) .gt. FEW_PC
         enddo      ! ibin = 1,NBIN
        endif       ! inucproc(iepart,ienucto) .eq. I_DROPFREEZE
       enddo       ! inuc = 1,nnuc2elem(iepart)
      enddo         ! igroup = 1,NGROUP
c
c
c  Return to caller with particle loss rates due to nucleation evaluated.
c
      return
      end
