       subroutine gsolve
c
c
c  @(#) gsolve.f  Ackerman  Dec-1995
c  This routine calculates new gas concentrations.
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter gsolve'
c
c
c  Update each gas concentration using gas production rates
c 
      do igas = 1,NGAS
        gc3(ixyz,igas) = gc3(ixyz,igas) + dtime * gasprod(igas)
        
        if (gc3(ixyz, igas) .lt. 0.0) then
          write(LUNOPRT, *) 'FATAL: Negative gas concentration for ',
     $      gassname(igas), '...'
          write(LUNOPRT, *) '      ixyz  = ', ixyz
          write(LUNOPRT, '(a,e)') ' (new) gc3   = ', gc3(ixyz, igas)
          write(LUNOPRT, '(a,e)') ' (old) gc3   = ', gc3(ixyz, igas) -
     $      dtime * gasprod(igas)
          write(LUNOPRT, '(a,e)') '   gasprod   = ', gasprod(igas)
          write(LUNOPRT, *) '     dtime  = ', dtime
          write(LUNOPRT, *) 'ntsubsteps  = ', ntsubsteps
          write(LUNOPRT, *) 'maxsubsteps = ', maxsubsteps
          write(LUNOPRT, *) '   supsati  = ', supsati3(ixyz,igas)
          write(LUNOPRT, *) '   supsatl  = ', supsatl3(ixyz,igas)
          write(LUNOPRT, *) '   pconmax  = ', pconmax(ixyz,:)
          write(LUNOPRT, *) '    conmax  = ', conmax
          
          call endcarma
        end if
      enddo
c
c
c  Return to caller with new gas concentrations.
c
      return
      end
