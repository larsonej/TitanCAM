      subroutine freezdropl
c
c
c  @(#) freezdropl.f  Jensen  Jan-2000
c  This routine evaluates particle loss rates due to nucleation <rnuclg>:
c  droplet freezing only.
c  
c  The loss rates for all particle elements in a particle group are equal.
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
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter freezdropl'
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
c  Only compute nucleation rate for droplet freezing
c
        if( inucproc(iepart,ienucto) .eq. I_DROPFREEZE ) then
c
c
c  Loop over particle bins.  
c
         do ibin = NBIN,1,-1
c
c
c  <inucto> is index of target nucleation bin.
c
          if( ignucto .ne. 0 )then
            inucto = inuc2bin(ibin,igroup,ignucto)
          else
            inucto = 0
          endif
c
c
c  Bypass calculation if few particles are present 
c
          if( pc3(ixyz,ibin,iepart) .gt. FEW_PC )then
c
c
c  Temporary simple kludge: Set <rnuclg> to 1.e2 if T < -40C
c  
           if( t3(ixyz) .lt. Tfreez(igrowgas(iepart)) ) then
             rnuclg(ibin,igroup,ignucto) = 1.e1
             inucstep(igroup) = ibin
           endif

          endif     ! pc(source particles) .gt. FEW_PC
         enddo      ! ibin = 1,NBIN
        endif       ! inucproc(iepart,ienucto) .eq. I_DROPFREEZE
       enddo        ! inuc = 1,nnuc2elem(iepart)
      enddo         ! igroup = 1,NGROUP
c
c
c  Return to caller with particle loss rates due to nucleation evaluated.
c
      return
      end
