      subroutine prestep
c
c
c  @(#) prestep.f  McKie  Oct-1995
c  This routine handles all preliminary setup at the beginning
c  of every timestep.  Things that would appropriately be done
c  here include:
c    Input or otherwise define interface quantities from other submodels.
c    Save any model state that is needed to compute tendencies.
c    Save any model state that might be needed for comparison at end of step.
c    Update timestep counter and simulation time.
c
c  NOTE: Within the CAM framework, allthe last vaules are maintained by CAM
c and have been set before calling this routine. CAM needs to keep this,
c since CARMA is only working on a column and the last timestep is
c probably from a diferent column. HAving CAm store the last values
c is also needed for restarts.
c
c NOTE: The variables: pold, told, rhoaold, zlod, zcold onl seemd to be
c used with varstep, so they have not been implemented. This saveson the
c the nukber of things that need to be stored by CAM.
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter prestep'
c
c
c  Don't allow particle concentrations to get too small.
c  
      do ixyz = 1,NXYZ
        do ibin = 1,NBIN
          do ielem = 1,NELEM
            call smallconc(ibin,ielem)
          enddo
        enddo
      enddo
c
c
c  Set <pcl> to <pc> from previous time step
c  and find maximum for each element.
c
      do ielem = 1,NELEM
        pcmax(ielem) = 0.
        do i = 1,NBIN
          do ixyz = 1,NXYZ
c            pcl(ixyz,i,ielem) = pc3(ixyz,i,ielem)
            pcmax(ielem) = max( pcmax(ielem), pc3(ixyz,i,ielem) )
          enddo
        enddo
      enddo
c
c
c  Find maximum particle concentration for each spatial grid box
c  (in units of cm^-3)
c
      do igroup = 1,NGROUP
        iep = ienconc(igroup)
        do ixyz = 1,NXYZ
          pconmax(ixyz,igroup) = 0.
          do i = 1,NBIN
            pconmax(ixyz,igroup) = max( pconmax(ixyz,igroup),
     $                                  pc3(ixyz,i,iep) )
          enddo
          pconmax(ixyz,igroup) = pconmax(ixyz,igroup) /
     $         ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )
        enddo
      enddo
c
c
c  Set <gcl> to <gc> from previous time step
c  and store old ice and liquid supersaturations for variable time stepping
c
c      do igas = 1,NGAS
c        do ixyz = 1,NXYZ
c          gcl(ixyz,igas) = gc3(ixyz,igas)
c          supsatlold(ixyz,igas) = supsatl3(ixyz,igas)
c          supsatiold(ixyz,igas) = supsati3(ixyz,igas)
c        enddo
c      enddo
c
c
c  Store old values of atmospheric properties for variable time stepping.
c
c      do ixyz = 1,NXYZ
c        zmetl3(ixyz) = zmet3(ixyz)
c      enddo

c      do ixyz = 1,NXYZ
c        ptcl(ixyz) = ptc3(ixyz)
c        told(ixyz) = t3(ixyz)
c        pold(ixyz) = p3(ixyz)
c        rhoaold(ixyz) = rhoa3(ixyz)
c      enddo

c      if( do_parcel )then
c        do ixyz = 1,NXYZP1
c          zlold(ixyz) = zl3(ixyz)
c        enddo
c        do ixyz = 1,NXYZ
c          zcold(ixyz) = zc3(ixyz)
c        enddo
c      endif
c
c
c  Set production terms and loss rates due to slow microphysics
c  processes (coagulation) to zero.
c
      if( ibtime .eq. 0 .or. do_coag )then

        do i = 1,NBIN

          do ielem = 1,NELEM
            do ixyz = 1,NXYZ
              coagpe(ixyz,i,ielem) = 0.
            enddo
          enddo

          do igroup = 1,NGROUP
            do ixyz = 1,NXYZ
              coaglg(ixyz,i,igroup) = 0.
            enddo
          enddo

        enddo

      endif
c
c
c  Set production terms and loss rates due to fast microphysics
c  processes to zero if just starting a simulation.
c
      call zeromicro
c
c
c  Temperature oscillation (hack for a test simulation)
c
c      if( itime .eq. 0 )then
c        print*,'warning: ugly T oscillation in prestep()'
c      endif
c      t_amp = 2.0
c      t_period = 1200.
c      t_omega = 2.*PI/t_period
c      cos_arg = t_omega*(time+dtime)
c      if( time .gt. 0.75*t_period )then                  ! ugliness
c        dtemp = 0.
c      else
c        dtemp = dtime * t_amp * t_omega * cos( cos_arg )
c      endif
c      do ixyz = 1,NXYZ
c        t3(ixyz) = t3(ixyz) - dtemp
c        pt = t3(ixyz) * ( PREF / p3(ixyz) )**RKAPPA
c        ptc3(ixyz) = pt * rhoa3(ixyz)
c      enddo
c
c
c  Return to caller with preliminary timestep things completed.
c
      return
      end
