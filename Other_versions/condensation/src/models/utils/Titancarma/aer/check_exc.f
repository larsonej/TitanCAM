      subroutine check_exc
c
c
c  @(#) check_exc.f  Barth  May-2003
c  Check that applying exchange term to gas production won't cause
c  <gc> to go negative
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables.
c
      include 'globaer.h'

c
c  Define formats
c
    1 format(a,2(i3,1pe11.4))
    2 format(a,i3,2x,'gas 1:',1pe13.6,3x,'gas 2:',1pe13.6)
    3 format(a,2x,'f:',f20.15,2x,'fold:',f20.15,2x,'df:',
     $       f20.15,2x,'fmax:',f20.15)
    4 format(a,i4,'   cloud: ',1pe13.6,0p,'  f: ',f17.15,
     $       '   gas 1: ',1pe13.6,'   df: ',1pe13.6,
     $       '  mass: ',1pe13.6,'  dt: ',0p,f10.5)
c
c
c  Announce entry into this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter check_exc'
c
c
c  Check if <gasprod> will cause <gc> to go negative.  If so, increment
c  df_factor by 1% and call {isotope_exc} again

      do igas = 1,NGAS
        do igrp = 1,NGROUP

          if( -gasprod(igas) .ge. gc3(ixyz,igas) / dtime ) then
            df_factor = df_factor - 0.01d0
            if(df_factor .lt. 0.d0) then
              print *, 'ISOTOPE_EXC: df reached 0!'
            else
              print *,itime,df_factor
              call isotope_exc
            endif
          endif

        enddo !igrp
      enddo   !igas
c
c
c  Return to caller with .
c
      return
      end
