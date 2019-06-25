       subroutine csolve(ibin,ielem)
c
c
c  @(#) csolve.f  McKie  Sep-1997
c  This routine calculates new particle concentrations from coagulation
c  microphysical processes.
c
c   The basic form from which the solution is derived is
c   ( new_value - old_value ) / dtime = source_term - loss_rate*new_value
c
c  This routine derived from psolve.f code, in which particle concentrations
c  due to coagulation were formerly included, before the relatively slow
c  coagulation calcs were separated from the other microphysical processes
c  so that time splitting could be applied to these fast & slow calcs.
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
    1 format('Days:',f13.8,2x,i3,2x,'bin:'i3,2x,'prod:',
     $       1pe13.6,2x,'loss:',1pe13.6,2x,0p,f20.15)
    2 format('Time:',f13.8,2x,i3,2x,'bin:',i3,2x,'cld :',
     $       1pe13.6,2x,'core:',1pe13.6,2x,'gcor:',1pe13.6)
    5 format(a,2(i4),2x,'gc/vol:',f15.10,2x,'gc/cld:',
     $       f15.10,2x,'cmf:',f15.10,2x)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter csolve'
c
c
c  Define current group & particle number concentration element indices
c
      igroup = igelem(ielem)         ! particle group
      iepart = ienconc(igroup)       ! particle number concentration element
c
c
c  Visit each spatial point
c
      do ixyz = 1,NXYZ
c
c
c  Metric scaling factor
c
        xyzmet = xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)
c
c
c  Compute total production rate due to coagulation
c
        ppd = coagpe(ixyz,ibin,ielem) / xyzmet
c
c
c  Compute total loss rate due to coagulation
c
        pls = coaglg(ixyz,ibin,igroup) / xyzmet
c
c
c  Update net particle number concentration during current timestep
c  due to production and loss rates for coagulation
c
        pc3(ixyz,ibin,ielem) = ( pc3(ixyz,ibin,ielem) + dtime*ppd ) /
     $                         ( ONE + pls*dtime )
c
c
c  Prevent particle concentrations from dropping below SMALL_PC
c
c       if(ixyz.eq.36 .and. ielem.eq.4) 
c    $    write(*,*) 'csolve',pc3(ixyz,ibin,ielem),ibin,ppd,pls

ccc     call smallconc(ibin,ielem)
c
c
c     if(do_write(12) .and. itime.gt.1) then
c      if(itype(ielem).eq.I_GROWCORE) then
c       iecore = icorelem(1,igroup)
c       volpart = pc3(ixyz,ibin,iepart)*rmass(ibin,igroup)
c    $               -  pc3(ixyz,ibin,iecore)
c       if(volpart .lt. 0.) 
c    $    write(*,2) time/60.**2/24.,ixyz,ibin,
c    $    pc3(ixyz,ibin,iepart)*rmass(ibin,igroup),
c    $    pc3(ixyz,ibin,iecore),pc3(ixyz,ibin,ielem)
c       if( pc3(ixyz,ibin,ielem) .ge. volpart )
c    $    write(*,1) time/60.**2/24.,ixyz,ibin,ppd,pls,
c    $    pc3(ixyz,ibin,ielem)/volpart
c      endif
c     endif

c      if(ixyz.eq.iwa) then
c      if(itype(ielem).eq.I_GROWCORE) then
c       iecore = icorelem(1,igroup)
c       volpart = pc3(ixyz,ibin,iepart)*rmass(ibin,igroup)
c    $               -  pc3(ixyz,ibin,iecore)
cc      if( pc3(ixyz,ibin,ielem) .ge. volpart )then
c         write(*,5) 'csolve fraction:',ixyz,ibin,
c    $       pc3(ixyz,ibin,ielem)/volpart,
c    $       pc3(ixyz,ibin,ielem)/
c    $         (pc3(ixyz,ibin,iepart)*rmass(ibin,igroup)),
c    $       pc3(ixyz,ibin,iecore)/
c    $         (pc3(ixyz,ibin,iepart)*rmass(ibin,igroup))
cc        write(*,*) 'Csolve decreasing growcore',ixyz,ibin
cc        pc3(ixyz,ibin,ielem) = 0.75d0*volpart
cc      endif
c      endif
c      endif

      enddo

c
c
c  Return to caller with new particle number concentrations.
c
      return
      end
