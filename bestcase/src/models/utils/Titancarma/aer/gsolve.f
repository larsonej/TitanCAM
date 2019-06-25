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
c  Define formats
c
    1 format(a,f5.2,' km ',i3,3x,'gc: ',1pe19.12,'  prod: ',1pe13.6,
     $       3x,0p,'ss: ',f13.8,3x,f16.8,' days ',i8,' itime')
    2 format(a,2(i4),'  gc: ',1pe11.4,'  gcl: ',1pe11.4,' prod: ',
     $       1pe11.4,'  flux: ',1pe11.4,'  grow: ',1pe11.4,a,
     $       '  drop: ',1pe11.4,'  mono: ',1pe11.4,'  poly: ',
     $       1pe11.4,'  dt: ',1pe11.4,0p,f15.8,' days ',i8,' itime')
    3 format(a,2(i4),'  orig: ',1pe11.4,'  indiv: ',1pe11.4,
     $       '  ss: ',0p,f13.8,3x,f13.8,' days')
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
        gc_old = gc3(ixyz,igas)

        gc3(ixyz,igas) = gc3(ixyz,igas) + dtime * gasprod(igas)
      
c         iwa = 14
c        if(ixyz .eq. iwa .and. itime.gt.14690) then 
c        if( itime.ge.25680 .and. itime.le.25800 ) then
c         write(*,1) 'GSOLVE ',zl3(ixyz)/1.d5,igas,gc3(ixyz,igas),
c    $           gasprod(igas),supsati3(ixyz,igas),time/60.**2/24.,
c    $           itime
c        endif

c       if( gasname(igas) .eq. 'ethane' .and. ixyz.eq.11 )
c    $     write(*,*) 'Ethane Vapor:',gc3(ixyz,igas)-gc_old,
c    $                 gasprod(igas)

c       if(ixyz.le.30)
c    $   write(*,3) 'gasprod',ixyz,igas,gasprod(igas),gasprod(igas),
c    $              supsati3(ixyz,igas),time/60.**2/24.

        igroup=2

        if(itime.gt.1 .and. dtime.le.dtmin .and. 
     $                                gc3(ixyz,igas) .lt. 0.) then
c !       if( abs(gc3(ixyz,igas)) .lt. ALMOST_ZERO ) then
c !         gc3(ixyz,igas) = ALMOST_ZERO * pvapi3(ixyz,igas) *
c !  $                       gwtmol(igas) / RGAS / t3(ixyz)
c !       else
           write(*,2) 'gsolve',ixyz,igas,gc3(ixyz,igas),
     $        gcl(ixyz,igas),gasprod(igas),fbotevap(igas),
     $        gprod_grow(igroup,igas),' Growcore: ',
     $        gprod_drop(igroup,2),
     $        gprod_mono(igroup,2),gprod_poly(igroup,2),
     $        dtime,time/60.**2/24.,itime
          stop 1
c !      endif
        endif

      enddo
c
c
c  Return to caller with new gas concentrations.
c
      return
      end
