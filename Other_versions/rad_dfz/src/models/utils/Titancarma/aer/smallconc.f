      subroutine smallconc(i,ie)
c
c
c  @(#) smallconc.f  Ackerman  Oct-1997
c  This routine ensures limits all particle concentrations in a grid box
c  to SMALL_PC.  In bins where this limitation results in the particle 
c  concentration changing, the core mass fraction and second moment fraction 
c  are set to <FIX_COREF>. 
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter smallconc'
   
        ig = igelem(ie)
        ip = ienconc(ig)


        if( ie .eq. ip )then
c
c
c  Element is particle concentration
c 
          pc3(ixyz,i,ie) = max( pc3(ixyz,i,ie), SMALL_PC )

        else
c
c
c  Element is core mass or core second moment
c 
            if( itype(ie) .eq. I_COREMASS ) then

              small_val = SMALL_PC*rmass(i,ig)*FIX_COREF

            elseif( itype(ie) .eq. I_GROWCORE .or.
     $          itype(ie) .eq. I_VOLCORE ) then

              small_val = SMALL_PC*rmass(i,ig)*FIX_COREF**2

            elseif( itype(ie) .eq. I_CORE2MOM )then

              small_val = SMALL_PC*(rmass(i,ig)*FIX_COREF)**2

            endif

            if( pc3(ixyz,i,ip) .le. SMALL_PC .or.
     $          pc3(ixyz,i,ie) .lt. small_val )then

              pc3(ixyz,i,ie) = small_val
c             pc3(ixyz,i,ie) = 0.

            endif

        endif
c
c
c  Return to caller with particle concentrations limited to SMALL_PC
c
      return
      end
