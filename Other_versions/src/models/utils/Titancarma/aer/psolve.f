       subroutine psolve(ibin,ielem)
c
c
c  @(#) psolve.f  Jensen  Oct-1995
c  This routine calculates new particle concentrations.
c
c   The basic form from which the solution is derived is
c   ( new_value - old_value ) / dtime = source_term - loss_rate*new_value
c
c  Modified  Sep-1997  (McKie)
c    New particle concentrations due to coagulation processes
c    were moved to the csolve routine.  Csolve is called to
c    update particle concentrations due to coagulation.
c    This new psolve now updates particle concentrations due
c    to the faster calcs of the non-coag microphysical processes.
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c  Argument list input:
c    ielem, ibin
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
    2 format('Production   Nucl: ',1pe11.4,' Grow: ',1pe11.4)
    3 format('Loss         Nucl: ',1pe11.4,' Grow: ',
     $        1pe11.4,' Evap: ',1pe11.4,/)
    4 format(2(i4),2x,3(1pe13.6,2x))
    5 format(2(i4),2x,'CH4 cloud: ',1pe13.6,'  core: ',
     $       1pe13.6,'  C2H6 core: ',1pe13.6,
     $       '  CH4 only: ',1pe13.6)
    6 format(a,i3,2x,1pe11.4,'  prod: (tot/nuc/grow)  ',3(1pe11.4,2x),
     $       '  loss: (tot/nuc/grow/evap)  ',4(1pe11.4,2x))
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter psolve'

c
c
c  Define current group & particle number concentration element indices
c
      igroup = igelem(ielem)
      iepart = ienconc(igroup)
c
c
c  Compute total production rate
c
        ppd = rnucpe(ibin,ielem) + rprod(ixyz)   !should be changed so that 
     $      + growpe(ibin,ielem)                 !rprod is only added for 13 A
     $      + evappe(ibin,ielem)                 !tholin particles (and moved)
c
c
c  Sum up nucleation loss rates
c
        rnuclgtot = 0.
        do igto = 1,NGROUP
          rnuclgtot = rnuclgtot + rnuclg(ibin,igroup,igto)
        enddo

c
c  Compute total loss rate
c
        pls = rnuclgtot
     $      + growlg(ibin,igroup) 
     $      + evaplg(ibin,igroup) 
c
c
c  Update net particle number concentration during current timestep
c  due to production and loss rates
c
        pc_old = pc3(ixyz,ibin,ielem)

        pc3(ixyz,ibin,ielem) = ( pc3(ixyz,ibin,ielem) + dtime*ppd ) /
     $                         ( ONE + pls*dtime )

c       if(ixyz.eq.7 .and. ibin.le.6 .and. ielem.eq.8)
c       if( zl3(ixyz).eq.40.d5 .and. ielem.eq.2 .and. ibin.le.15 )
cc      if( ixyz.eq.9 ) then
c         if( ielem .eq. 5 ) 
c    $   write(*,6) 'psolve',ielem,pc3(ixyz,ibin,ielem),ppd,
c    $    rnucpe(ibin,ielem),growpe(ibin,ielem),pls,rnuclgtot,
c    $    growlg(ibin,igroup),evaplg(ibin,igroup)
c         if( ielem .eq. 9 ) 
c    $   write(*,6) 'psolve',ielem,pc3(ixyz,ibin,ielem),ppd,
c    $    rnucpe(ibin,ielem),growpe(ibin,ielem),pls,rnuclgtot,
c    $    growlg(ibin,igroup),evaplg(ibin,igroup)
cc      endif
c       if( itype(ielem) .eq. I_GROWCORE .and. 
c    $      rnucpe(ibin,ielem) .ne. 0          ) 
c    $    write(*,*) 'Prod: ',ibin,dtime*ppd

c       if( ielem .eq. 2 .and. rnuclgtot .ne. 0. ) 
c    $   write(*,*) 'Loss: ',ibin,(rmass(ibin,igroup)*pc_old -
c    $                       pc3(ixyz,ibin,3) )*
c    $                         (ONE/(ONE + pls*dtime) - ONE)

c       if( itype(ielem) .eq. I_GROWCORE .and. 
c    $      rnucpe(ibin,ielem) .ne. 0          ) then
c         write(*,5) ixyz,ibin,rnucpe(ibin,5)*rmass(ibin,igroup),
c    $               rnucpe(ibin,6),rnucpe(ibin,ielem),
c    $               rnucpe(ibin,5)*rmass(ibin,igroup)
c    $                 -rnucpe(ibin,6)-rnucpe(ibin,ielem)
c       endif
       
c
c       if(igroup.eq.2 .and. ixyz.eq.iwa)
c    $   write(*,*) 'rnuclg:',ielem,ibin,rnuclgtot,
c    $              pc_old,pc3(ixyz,ibin,ielem)
c
c       if( do_write(3) ) then
c        if( ixyz .eq. iwa .and. ibin .eq. iwb ) then
c        if( ixyz .eq. iwa .and. dtmax .lt. 4.d3 ) then
c         if(itype(ielem) .eq. I_VOLATILE .or.
c    $       itype(ielem) .eq. I_GROWCORE ) then
c
c          igas = igrowgas(ielem)
c          if(itype(ielem) .eq. I_VOLATILE) then
c           write(*,1) 'New g/cm3 for elem',ielem,
c    $             rmass(ibin,2)*pc3(ixyz,ibin,ielem),
c    $             zl3(ixyz)/1.d5,ibin,supsati3(ixyz,igas)
c          else
c           write(*,1) 'New g/cm3 for elem',ielem,pc3(ixyz,ibin,ielem),
c    $             zl3(ixyz)/1.d5,ibin,supsati3(ixyz,igas)
c          endif
c
c         iwa = 14
c         if(ixyz.eq.iwa .and. icomp(ielem).ne.I_CxHyNx .and. 
c    $       (itime.gt.14689.and.itime.lt.14696) .and. 
c    $        rnucpe(ibin,6).gt.ZERO .and.
c    $         itype(ielem).eq.I_VOLCORE .and.
c    $          ielem.eq.7 .and.
c    $      (ibin.ge.35.and.ibin.le.48) ) then
c          write(*,*) 'Element: ',ielem,ibin,itime,pc3(ixyz,ibin,ielem)
c          write(*,2) rnucpe(ibin,ielem),growpe(ibin,ielem)
c          write(*,3) rnuclgtot,growlg(ibin,igroup),evaplg(ibin,igroup)
c          write(*,*) pc_old
c         endif
c
c         endif
c        endif
c       endif
c
c  Return to caller with new particle number concentrations.
c
      return
      end
