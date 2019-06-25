       subroutine downgevapply
c
c
c  @(#) downgevapply.f  Ackerman  Dec-1995
c  This routine applies evaporation and nucleation production terms to
c  particle concentrations.
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
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
c  Define formats
c
    1 format(a,i2,',',1pe13.6,2x,0p,f5.2,' km, bin ',i2,
     $       0p,f20.10)
    2 format('Production   Nucl: ',1pe11.4,' Evap: ',1pe11.4,/)
    3 format('EvapProd',3(i4),3x,'s: ',0p,f11.5,3x,
     $       'new pc: ',1pe11.4,3x,
     $       'old pc: ',1pe11.4,3x,'evappe: ',1pe11.4,3x,
     $       'rnucpe: ',1pe11.4,2x,0p,f15.8,' days',i10) 
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter downgevapply'
c
c  Visit each radius bin for each element to compute particle production 
c  due to evaporation and element transfer processes for which the source
c  element number is greater than the target element number
c
      pcsmall = 1000.d0
      ibinsmall = 0
      ixyzsmall = 0

      do ielem = 1,NELEM
        do ibin = 1,NBIN

            pc_old = pc3(ixyz,ibin,ielem)

            pc3(ixyz,ibin,ielem) = pc3(ixyz,ibin,ielem) +
     $                         dtime * ( evappe(ibin,ielem) +
     $                                   rnucpe(ibin,ielem) )

c           if( itype(ielem) .eq. I_GROWCORE ) then
c           if( itype(ielem) .eq. I_VOLATILE ) then
c           if( itype(ielem).eq.I_INVOLATILE .and. itime.gt.52535) then
c           if( itype(ielem).eq.I_VOLCORE .and. itime.gt.240709) then
c       if(itype(ielem).eq.I_VOLCORE.and.evappe(ibin,ielem).gt.0.) then
c            if( t3(ixzy).gt.Tfreez(1) ) then
c              ss = supsatl3(ixyz,1)
c            else
c              ss = supsati3(ixyz,1)
c            endif
c            if(ss.lt.ZERO) 
c            if( ielem.eq.7 .and. itime.gt.3604.and.itime.lt.3607
c    $           .and. ibin.ge.20.and.ibin.le.22 .and.ixyz.eq.14)
c            if(ixyz.eq.11.and.ibin.ge.36)
c    $        write(*,3) ixyz,ielem,ibin,ss,
c    $                 pc3(ixyz,ibin,ielem),pc_old,
c    $                 evappe(ibin,ielem),rnucpe(ibin,ielem),
c    $       pc3(ixyz,ibin,ielem)/pc3(ixyz,ibin,ielem-1)/rmass(ibin,2),
c    $       itime
c    $                 time/60.**2/24.,itime
c           endif !write
c
        enddo
      enddo

c
c  Return to caller with evaporation and down-grid element transfer
c  production terms applied to particle concentrations.
c
      return
      end
