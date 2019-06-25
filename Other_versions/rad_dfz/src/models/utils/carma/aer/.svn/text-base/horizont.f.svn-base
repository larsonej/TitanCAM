      subroutine horizont
c
c
c  @(#) horizont.f  Jensen  Jan-1997
c  This routine drives the horizontal transport calculations.
c  Following Lin and Rood [Mon. Wea. Rev., 124, 2046, 1996], we
c  are using a flux-form scheme for the outer advection operator
c  and an advective-form scheme for the inner operator.  For
c  the flux-form scheme, we are using the Piecewise Parabolic
c  method [Colela and Woodard, J. Comp. Phys., 54, 174-201, 1984]
c  For the advective-form scheme, we are using a first-order
c  semi-Lagrangian technique [Bates and McDonald, Mon. Wea. Rev.,
c  110, 1831-1842, 1982].
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
c
c  Local declarations
c
      dimension doxn(NX,NY),doyn(NX,NY),
     $          dixn(NX,NY),diyn(NX,NY)
c
c
c  Define formats
c
    3 format(i3,3x,6e13.4)
    4 format(i3,3x,i3,3x,i3,3x,6f9.2)
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter horizont'
c
c
c  Set flag to make sure coefficients are calculated in
c  initial call to htranglk.f
c
      do_ccoef = .true.
c
c
c  Initially, set concentration changes to 0.
c
      do iy = 1,NY
        do ix = 1,NX
          dixn(ix,iy) = 0.
          diyn(ix,iy) = 0.
          doxn(ix,iy) = 0.
          doyn(ix,iy) = 0.
        enddo
      enddo

      do k = 1,NZ     ! Loop over vertical levels
c
c
c  If using Galerkin scheme, calculate coefficients.
c
      if( ihoradv .eq. I_GALERKIN )then
        call glkcoef(k)
      endif
c
c=======================================================================
c
c  Transport of particles
c
        do ielem = 1,NELEM           ! Loop over particle elements

          ig = igelem(ielem)         ! particle group

          do ibin = 1,NBIN          ! Loop over particle mass bins
c
c  East-west transport by inner advective-form operator.
c
            if( do_ew .and. do_ns .and. ihoradv .eq. I_PPM ) then

              idir = IDIRX

              do iy = 1,NY       ! Loop over latitudes

                do ix = 1,NX
                  htrans(ix) = u(ix,iy,k)
                  dhor(ix) = dx(ix,iy,k)
                  chor(ix) = pc(ix,iy,k,ibin,ielem) / rhoa(ix,iy,k)
                enddo

                call htranfosl(idir)

                do ix = 1,NX
                  dixn(ix,iy) = chor(ix)*rhoa(ix,iy,k) -
     $                          pc(ix,iy,k,ibin,ielem)
                enddo

              enddo       ! <iy=1,NY>

            endif         ! end of if <do_ew>
c
c  North-south transport by inner advective-form operator.
c
            if( do_ns .and. do_ew .and. ihoradv .eq. I_PPM ) then

              idir = IDIRY

              do ix = 1,NX       ! Loop over latitudes

                do iy = 1,NY
                  htrans(iy) = v(ix,iy,k)
                  dhor(iy) = dy(ix,iy,k)
                  chor(iy) = pc(ix,iy,k,ibin,ielem) / rhoa(ix,iy,k)
                enddo

                call htranfosl(idir)

                do iy = 1,NY
                  diyn(ix,iy) = chor(iy)*rhoa(ix,iy,k) -
     $                          pc(ix,iy,k,ibin,ielem)
                enddo

              enddo       ! <ix=1,NX>

            endif         ! end of if <do_ns>
c
c  East-west transport by outer flux-form operator.
c
            if( do_ew ) then

              idir = IDIRX

              do iy = 1,NY       ! Loop over latitudes

                do ix = 1,NX
                  htrans(ix) = u(ix,iy,k)
                  dhor(ix) = dx(ix,iy,k)
                  chor(ix) = pc(ix,iy,k,ibin,ielem) + 0.5*diyn(ix,iy)
                enddo

                if( ihoradv .eq. I_PPM )then
                  call htranppm(idir)
                else if( ihoradv .eq. I_GALERKIN )then
                  call htranglk(idir,iy)
                endif

                do ix = 1,NX
                  doxn(ix,iy) = chor(ix) - pc(ix,iy,k,ibin,ielem) -
     $                          0.5*diyn(ix,iy)
                enddo

              enddo       ! <iy=1,NY>

            endif         ! end of if <do_ew>
c
c  North-south transport by outer flux-form operator.
c
            if( do_ns ) then

              idir = IDIRY

              do ix = 1,NX       ! Loop over longitudes

                do iy = 1,NY
                  htrans(iy) = v(ix,iy,k)
                  dhor(iy) = dy(ix,iy,k)
                  chor(iy) = pc(ix,iy,k,ibin,ielem) + 0.5*dixn(ix,iy)
                  if( ihoradv .eq. I_GALERKIN )
     $              chor(iy) = chor(iy) + doxn(ix,iy)
                enddo

                if( ihoradv .eq. I_PPM )then
                  call htranppm(idir)
                else if( ihoradv .eq. I_GALERKIN )then
                  call htranglk(idir,ix)
                endif

                do iy = 1,NY
                  doyn(ix,iy) = chor(iy) - pc(ix,iy,k,ibin,ielem) -
     $                          0.5*dixn(ix,iy)
                  if( ihoradv .eq. I_GALERKIN )
     $              doyn(ix,iy) = chor(iy) - pc(ix,iy,k,ibin,ielem) -
     $                            doxn(ix,iy)
                enddo

              enddo       ! <ix=1,NX>

            endif         ! end of if <do_ns>
c
c
c  Update particle concentrations
c
            do ix = 1,NX
              do iy = 1,NY
                pc(ix,iy,k,ibin,ielem) = pc(ix,iy,k,ibin,ielem) + 
     $                                   doxn(ix,iy) + doyn(ix,iy)
              enddo
            enddo
c
c
c  Horizontal diffusion of particles
c
c
c  East-west diffusion.
c
            if( do_ew ) then

              idir = IDIRX

              do iy = 1,NY       ! Loop over latitudes

                do ix = 1,NX
                  hdiff(ix) = dkx(ix,iy,k)
                  dhor(ix) = dx(ix,iy,k)
                  chor(ix) = pc(ix,iy,k,ibin,ielem)
                enddo

                call hordif(idir)

                do ix = 1,NX
                  pc(ix,iy,k,ibin,ielem) = chor(ix)
                enddo

              enddo       ! <iy=1,NY>

            endif         ! end of if <do_ew>
c
c  North-south diffusion.
c
            if( do_ns ) then

              idir = IDIRY

              do ix = 1,NX       ! Loop over longitudes

                do iy = 1,NY
                  hdiff(iy) = dky(ix,iy,k)
                  dhor(iy) = dy(ix,iy,k)
                  chor(iy) = pc(ix,iy,k,ibin,ielem)
                enddo

                call hordif(idir)

                do iy = 1,NY
                  pc(ix,iy,k,ibin,ielem) = chor(iy)
                enddo

              enddo       ! <ix=1,NX>

            endif         ! end of if <do_ns>

          enddo           ! <ibin=1,NBIN>

        enddo             ! <ielem=1,NELEM>
c
c=======================================================================
c
c  Transport of gases
c
        do igas = 1,NGAS           ! Loop over gases
c
c  East-west transport by inner advective-form operator.
c
          if( do_ew .and. do_ns .and. ihoradv .eq. I_PPM ) then

            idir = IDIRX

            do iy = 1,NY       ! Loop over latitudes

              do ix = 1,NX
                htrans(ix) = u(ix,iy,k)
                dhor(ix) = dx(ix,iy,k)
                chor(ix) = gc(ix,iy,k,igas) / rhoa(ix,iy,k)
              enddo

              call htranfosl(idir)

              do ix = 1,NX
                dixn(ix,iy) = chor(ix)*rhoa(ix,iy,k) -
     $                        gc(ix,iy,k,igas)
              enddo

            enddo       ! <iy=1,NY>

          endif         ! end of if <do_ew>
c
c  North-south transport by inner advective-form operator.
c
          if( do_ns .and. do_ew .and. ihoradv .eq. I_PPM ) then

            idir = IDIRY

            do ix = 1,NX       ! Loop over longitudes

              do iy = 1,NY
                htrans(iy) = v(ix,iy,k)
                dhor(iy) = dy(ix,iy,k)
                chor(iy) = gc(ix,iy,k,igas) / rhoa(ix,iy,k)
              enddo

              call htranfosl(idir)

              do iy = 1,NY
                diyn(ix,iy) = chor(iy)*rhoa(ix,iy,k) -
     $                        gc(ix,iy,k,igas)
              enddo

            enddo       ! <ix=1,NX>

          endif         ! end of if <do_ns>
c
c  East-west transport by outer flux-form operator.
c
          if( do_ew ) then

            idir = IDIRX

            do iy = 1,NY       ! Loop over latitudes

              do ix = 1,NX
                htrans(ix) = u(ix,iy,k)
                dhor(ix) = dx(ix,iy,k)
                chor(ix) = gc(ix,iy,k,igas) + 0.5*diyn(ix,iy)
              enddo

              if( ihoradv .eq. I_PPM )then
                call htranppm(idir)
              else if( ihoradv .eq. I_GALERKIN )then
                call htranglk(idir,iy)
              endif

              do ix = 1,NX
                doxn(ix,iy) = chor(ix) - gc(ix,iy,k,igas) -
     $                        0.5*diyn(ix,iy)
              enddo

            enddo       ! <iy=1,NY>

          endif         ! end of if <do_ew>
c
c  North-south transport by outer flux-form operator.
c
          if( do_ns ) then

            idir = IDIRY

            do ix = 1,NX       ! Loop over longitudes

              do iy = 1,NY
                htrans(iy) = v(ix,iy,k)
                dhor(iy) = dy(ix,iy,k)
                chor(iy) = gc(ix,iy,k,igas) + 0.5*dixn(ix,iy)
                if( ihoradv .eq. I_GALERKIN )
     $              chor(iy) = chor(iy) + doxn(ix,iy)
              enddo

              if( ihoradv .eq. I_PPM )then
                call htranppm(idir)
              else if( ihoradv .eq. I_GALERKIN )then
                call htranglk(idir,ix)
              endif

              do iy = 1,NY
                doyn(ix,iy) = chor(iy) - gc(ix,iy,k,igas) -
     $                        0.5*dixn(ix,iy)
                if( ihoradv .eq. I_GALERKIN )
     $            doyn(ix,iy) = chor(iy) - gc(ix,iy,k,igas) -
     $                          doxn(ix,iy)
              enddo

            enddo       ! <ix=1,NX>

          endif         ! end of if <do_ns>
c
c
c  Update gas concentrations
c
          do ix = 1,NX
            do iy = 1,NY
              gc(ix,iy,k,igas) = gc(ix,iy,k,igas) +
     $                           doxn(ix,iy) + doyn(ix,iy)
            enddo
          enddo
c
c
c  Horizontal diffusion of gases
c
          if( do_ew ) then

            idir = IDIRX

            do iy = 1,NY       ! Loop over latitudes

              do ix = 1,NX
                hdiff(ix) = dkx(ix,iy,k)
                dhor(ix) = dx(ix,iy,k)
                chor(ix) = gc(ix,iy,k,igas)
              enddo

              call hordif(idir)

              do ix = 1,NX
                gc(ix,iy,k,igas) = chor(ix)
              enddo

            enddo       ! <iy=1,NY>

          endif         ! end of if <do_ew>

          if( do_ns ) then

            idir = IDIRY

            do ix = 1,NX       ! Loop over longitudes

              do iy = 1,NY
                hdiff(iy) = dky(ix,iy,k)
                dhor(iy) = dy(ix,iy,k)
                chor(iy) = gc(ix,iy,k,igas)
              enddo

              call hordif(idir)

              do iy = 1,NY
                gc(ix,iy,k,igas) = chor(iy)
              enddo

            enddo       ! <ix=1,NX>

          endif         ! end of if <do_ns>

        enddo             ! <igas=1,NGAS>
c
c=======================================================================
c
c  Transport of potential temperature
c
c
c  East-west transport by inner advective-form operator.
c
        if( do_ew .and. do_ns .and. ihoradv .eq. I_PPM ) then

          idir = IDIRX

          do iy = 1,NY       ! Loop over latitudes

            do ix = 1,NX
              htrans(ix) = u(ix,iy,k)
              dhor(ix) = dx(ix,iy,k)
              chor(ix) = ptc(ix,iy,k) / rhoa(ix,iy,k)
            enddo

            call htranfosl(idir)

            do ix = 1,NX
              dixn(ix,iy) = chor(ix)*rhoa(ix,iy,k) -
     $                      ptc(ix,iy,k)
            enddo

          enddo       ! <iy=1,NY>

        endif         ! end of if <do_ew>
c
c  North-south transport by inner advective-form operator.
c
        if( do_ns .and. do_ew .and. ihoradv .eq. I_PPM ) then

          idir = IDIRY

          do ix = 1,NX       ! Loop over longitudes

            do iy = 1,NY
              htrans(iy) = v(ix,iy,k)
              dhor(iy) = dy(ix,iy,k)
              chor(iy) = ptc(ix,iy,k) / rhoa(ix,iy,k)
            enddo

            call htranfosl(idir)

            do iy = 1,NY
              diyn(ix,iy) = chor(iy)*rhoa(ix,iy,k) - 
     $                      ptc(ix,iy,k)
            enddo

          enddo       ! <ix=1,NX>

        endif         ! end of if <do_ns>
c
c  East-west transport by outer flux-form operator.
c
        if( do_ew ) then

          idir = IDIRX

          do iy = 1,NY       ! Loop over latitudes

            do ix = 1,NX
              htrans(ix) = u(ix,iy,k)
              dhor(ix) = dx(ix,iy,k)
              chor(ix) = ptc(ix,iy,k) + 0.5*diyn(ix,iy)
            enddo

            if( ihoradv .eq. I_PPM )then
              call htranppm(idir)
            else if( ihoradv .eq. I_GALERKIN )then
              call htranglk(idir,iy)
            endif

            do ix = 1,NX
              doxn(ix,iy) = chor(ix) - ptc(ix,iy,k) -
     $                      0.5*diyn(ix,iy)
            enddo

          enddo       ! <iy=1,NY>

        endif         ! end of if <do_ew>
c
c  North-south transport by outer flux-form operator.
c
        if( do_ns ) then

          idir = IDIRY

          do ix = 1,NX       ! Loop over longitudes

            do iy = 1,NY
              htrans(iy) = v(ix,iy,k)
              dhor(iy) = dy(ix,iy,k)
              chor(iy) = ptc(ix,iy,k) + 0.5*dixn(ix,iy) 
              if( ihoradv .eq. I_GALERKIN )
     $          chor(iy) = chor(iy) + doxn(ix,iy)
            enddo

            if( ihoradv .eq. I_PPM )then
              call htranppm(idir)
            else if( ihoradv .eq. I_GALERKIN )then
              call htranglk(idir,ix)
            endif

            do iy = 1,NY
              doyn(ix,iy) = chor(iy) - ptc(ix,iy,k) -
     $                      0.5*dixn(ix,iy)
              if( ihoradv .eq. I_GALERKIN )
     $          doyn(ix,iy) = chor(iy) - ptc(ix,iy,k) -
     $                        doxn(ix,iy)
            enddo

          enddo       ! <ix=1,NX>

        endif         ! end of if <do_ns>
c
c
c  Update potential temperatures
c
        do ix = 1,NX
          do iy = 1,NY
            ptc(ix,iy,k) = ptc(ix,iy,k) + 
     $                     doxn(ix,iy) + doyn(ix,iy)
          enddo
        enddo
c
c
c  Horizontal diffusion of potential temperatures
c
        if( do_ew ) then

          idir = IDIRX

          do iy = 1,NY       ! Loop over latitudes

            do ix = 1,NX
              hdiff(ix) = dkx(ix,iy,k)
              dhor(ix) = dx(ix,iy,k)
              chor(ix) = ptc(ix,iy,k)
            enddo

            call hordif(idir)

            do ix = 1,NX
              ptc(ix,iy,k) = chor(ix)
            enddo

          enddo       ! <iy=1,NY>

        endif         ! end of if <do_ew>

        if( do_ns ) then

          idir = IDIRX

          do ix = 1,NX       ! Loop over longitudes

            do iy = 1,NY
              hdiff(iy) = dky(ix,iy,k)
              dhor(iy) = dy(ix,iy,k)
              chor(iy) = ptc(ix,iy,k)
            enddo

            call hordif(idir)

            do iy = 1,NY
              ptc(ix,iy,k) = chor(iy)
            enddo

          enddo       ! <ix=1,NX>

        endif         ! end of if <do_ns>
c
c=======================================================================
c
c  Transport of air density
c
c
c  East-west transport by outer flux-form operator
c
        if( do_ew ) then

          idir = IDIRX

          do iy = 1,NY       ! Loop over latitudes

            do ix = 1,NX
              htrans(ix) = u(ix,iy,k)
              dhor(ix) = dx(ix,iy,k)
              chor(ix) = rhoa(ix,iy,k)
            enddo

            if( ihoradv .eq. I_PPM )then
              call htranppm(idir)
            else if( ihoradv .eq. I_GALERKIN )then
              call htranglk(idir,iy)
            endif

            do ix = 1,NX
              doxn(ix,iy) = chor(ix) - rhoa(ix,iy,k)
            enddo

          enddo       ! <iy=1,NY>

        endif         ! end of if <do_ew>
c
c  North-south transport by outer flux-form operator.
c
        if( do_ns ) then

          idir = IDIRY

          do ix = 1,NX       ! Loop over longitudes

            do iy = 1,NY
              htrans(iy) = v(ix,iy,k)
              dhor(iy) = dy(ix,iy,k)
              chor(iy) = rhoa(ix,iy,k)
            enddo

            if( ihoradv .eq. I_PPM )then
              call htranppm(idir)
            else if( ihoradv .eq. I_GALERKIN )then
              call htranglk(idir,ix)
            endif

            do iy = 1,NY
              doyn(ix,iy) = chor(iy) - rhoa(ix,iy,k)
            enddo

          enddo       ! <ix=1,NX>

        endif         ! end of if <do_ns>
c
c
c  Update air density
c
        do ix = 1,NX
          do iy = 1,NY
            rhoa(ix,iy,k) = rhoa(ix,iy,k) + 
     $                      doxn(ix,iy) + doyn(ix,iy)
          enddo
        enddo

      enddo               ! <k=1,NZ>
c
c
c  Return to caller with new particle and gas concentrations.
c
      return
      end
