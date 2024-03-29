#include <misc.h>
#include <params.h>


subroutine debug_srfflx_state2d(str)

   use physics_types,  only: physics_state, physics_tend
   use comsrf
   use mpishorthand
   use ppgrid,        only: pcols, pver
   use pmgrid      , only : iam

   implicit none

   ! input
   character(*) str

   ! local var
   integer m, i,c, ierr, ncol, counter

   save counter

   do m=0,1
#if (defined SPMD)
      call mpi_barrier(mpicom,ierr) 
#endif
!      call mpi_barrier(mpicom,ierr)
      if (m .eq. iam) then
         print *, 'DBUG: <proc,counter> ', str, iam, counter
         do c=begchunk, endchunk
            ncol = get_ncols_p(c)
            do i=1,ncol
               write(6,'(a12,i5, i5, e12.4, e12.4)') &
                    'DBUG: (ts)', c, i, srfflx_parm2d(c)%ts(i), srfflx_state2d(c)%ts(i)
            enddo
         enddo
      endif
   enddo

   counter = counter + 1

end subroutine debug_srfflx_state2d


subroutine physpkg(phys_state, gw, ztodt, phys_tend, pbuf)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Loop over time, calling driving routines for physics
! 
! Method: 
! COUP_CSM and must be checked in order to invoke the proper calling
! sequence for running the CSM model
! 
! Author: 
! Original version:  CCM3
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plon, plat, masterproc
   use tracers,    only: set_state_pdry
   use ppgrid       , only: pcols, pver
   use buffer, only: pblht, tpert, qpert, qrs, qrl
   use check_energy, only: check_energy_gmean
   use comsrf
   use comsrfdiag
   use phys_gmean, only: diur_av
#ifdef DIUR_AV
   use tmp_qrs_mod, only: tmp_qrs
#endif
#ifdef COUP_CSM
   use ccsm_msg, only: ccsmave, dorecv, dosend, ccsmsnd, ccsmrcv
#else
   use atm_lndMod, only: atmlnd_drv
#endif
#ifdef SPMD
   use mpishorthand
#endif
   use phys_buffer,    only: pbuf_size_max, pbuf_fld, pbuf_allocate, pbuf_deallocate, &
                             pbuf_update_tim_idx
   use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
   use physics_types,  only: physics_state, physics_tend
   use diagnostics,    only: diag_surf
   use time_manager,   only: get_nstep, is_first_step, is_first_restart_step, &
                             is_end_curr_month, get_curr_date, is_end_curr_day
   use physconst, only: stebol,zvir, rair, gravit
   use dycore, only: dycore_is
! AJF 3-28-08:  Add taulw_cam module containing interpck_cam, interpcia_cam, etc.
   use taulw_cam,      only: interpck_cam, interpcia_cam, get_taulw_cam
   use radcnst,        only: taui_lw  !Big lw-optical depth array 


#if (!defined COUP_CSM)
   use ice_constants, only: TfrezK
#endif
   use history, only: outfld
#if (!defined COUP_CSM)
#if (!defined COUP_SOM)
   use sst_data, only: sstint
   use ice_data, only: iceint
#endif
#endif

!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: gw(plat)                    ! Gaussian weights
   real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
!
! Input/Output arguments
!
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf
!
!---------------------------Local workspace-----------------------------
!
   integer :: i,m,lat,c,lchnk                   ! indices
   integer :: lats(pcols)                       ! array of latitude indices
   integer :: lons(pcols)                       ! array of longitude indices
   integer :: ncol                              ! number of columns
   integer :: nstep                             ! current timestep number
   integer :: ncdate                            ! current date in integer format [yyyymmdd]
   integer :: ncsec                             ! current time of day [seconds]
   integer :: yr, mon, day                      ! year, month, and day components of a date

   real(r8) fsds(pcols,begchunk:endchunk)        ! Surface solar down flux
   real(r8) :: tmp(pcols,begchunk:endchunk)
!AJF, 9/07/07 - for diurnal averaging:
   real(r8) :: tmpz(pcols,begchunk:endchunk,pver)
   real(r8), dimension(pcols,begchunk:endchunk) :: tmp_sols,tmp_solsd
!AJF, 4/24/08 - for debugging tau_lw
   logical, parameter :: debug_tau=.false.
                                                   
!-----------------------------------------------------------------------

   call t_startf ('physpkg_st')
   nstep = get_nstep()

   call pbuf_allocate('physpkg')

! Compute total energy of input state and previous output state
   call t_startf ('chk_en_gmean')
   call check_energy_gmean(phys_state, pbuf, ztodt, nstep)
   call t_stopf ('chk_en_gmean')

!-----------------------------------------------------------------------
! Advance time information
!-----------------------------------------------------------------------

   call advnce()
   call t_stopf ('physpkg_st')
!
!  set the state vector dry quantities
!

   if ( .not. dycore_is('LR') )then
      ! for LR, this is done in d_p_coupling since dynamics is called first

!$OMP PARALLEL DO PRIVATE (C,NCOL)

      do c=begchunk, endchunk
         call set_state_pdry( phys_state(c)) 
      end do

   endif ! .not. dycore_is('LR')

#ifdef TRACER_CHECK
   call gavglook ('before tphysbc DRY', phys_state, gw)
#endif


!-----------------------------------------------------------------------
! Tendency physics before flux coupler invokation
!-----------------------------------------------------------------------
!

#if (defined BFB_CAM_SCAM_IOP )
   do c=begchunk, endchunk
      call outfld('Tg',srfflx_state2d(c)%ts,pcols   ,c     )
   end do
#endif

#define KLUDGE_TREF
#ifdef KLUDGE_TREF
!!!  AJF 7/11/07:  Setting tref to ts on first step or first restart step
      if ( is_first_step() .or. is_first_restart_step() ) then
         do c=begchunk, endchunk
            srfflx_state2d(c)%tref=srfflx_state2d(c)%ts
         enddo
      endif
#endif

   call t_startf ('bc_physics')
   
! AJF, 3-28-08:  Interpolate reference longwave Ck, Cia coefficients to model
!                global-mean pressure-temperature profile.  Can do this 
!                periodically by setting "iradae" in namelist

   if (doabsems) then
   
     call t_startf ('get_taulw')

     call interpck_cam(phys_state)
     call interpcia_cam
     
     do c=begchunk, endchunk
      ncol=phys_state(c)%ncol
      do i=1,ncol      
       call get_taulw_cam(phys_state(c)%q(i,:,:), &
                          phys_state(c)%pint(i,:), & 
                          phys_state(c)%t(i,:),  &
                          taui_lw(c,i,:,:,:), & 
                          debug_tau)
      enddo
     enddo
     
     call t_stopf ('get_taulw')

   endif
   
 

!$OMP PARALLEL DO PRIVATE (C,NCOL)

   do c=begchunk, endchunk

!!!#define FAO_TEST_CONST  !ajf, 8-23-07: tref, ts look ok after restart
#ifdef FAO_TEST_CONST      !but still need kludge to initialize tref upon restart
      print*, 'FAO_UUU', c, srfflx_state2d(c)%tref(4), srfflx_state2d(c)%ts(4)
!!!      do i=1,pcols
!!!         srfflx_state2d(c)%tref(i) = srfflx_state2d(c)%ts(i)
!!!      enddo
#endif

      call t_startf ('tphysbc')

      call tphysbc(ztodt, pblht(1,c), tpert(1,c),             &
	           srfflx_state2d(c)%ts, srfflx_state2d(c)%sst,        &
                   qpert(1,1,c), surface_state2d(c)%precl,               &
	   	   surface_state2d(c)%precc, surface_state2d(c)%precsl,&
                   surface_state2d(c)%precsc,                          &
                   srfflx_state2d(c)%asdir, srfflx_state2d(c)%asdif,    &
                   srfflx_state2d(c)%aldir, srfflx_state2d(c)%aldif,  &
                   snowhland(1,c),                                    &
                   qrs(1,1,c), qrl(1,1,c), surface_state2d(c)%flwds,     &
                   fsns(1,c), fsnt(1,c),                               &
                   flns(1,c),    flnt(1,c), srfflx_state2d(c)%lwup,       &
                   surface_state2d(c)%srfrad, surface_state2d(c)%sols,    &
                   surface_state2d(c)%soll, surface_state2d(c)%solsd,   &
                   surface_state2d(c)%solld,                           &
                   phys_state(c), phys_tend(c),                        &
                   pbuf, prcsnw(1,c), fsds(1,c), landm(1,c), landfrac(1,c), &
                   ocnfrac(1,c),icefrac(1,c), &
                   srfflx_state2d(c)%tref, &
                   srfflx_state2d(c)%wsx, srfflx_state2d(c)%wsy)
                        
      call t_stopf ('tphysbc')
      if (dosw .or. dolw) then
	call output_flns_fsns_fluxes(surface_state2d(c),c)
      end if	

#ifdef DIUR_AV  !If using diurnal-averaging, store chunked data in proper format

      tmpz(:,c,:)=tmp_qrs(:,:)
      tmp_sols(:,c)=surface_state2d(c)%sols(:)
      tmp_solsd(:,c)=surface_state2d(c)%solsd(:)

#endif



#if ( ! defined COUP_CSM )
!
! zero surface fluxes at beginning of each time step.  Land Ocean and Ice
! processes will will write into process specific flux variables
! at the end of the time step these separate fluxes will be combined over the
! entire grid
!
!!!    call srfflx_state_reset (srfflx_state2d(c)) !ajf 7/12- don't call
                                                    !when land driver off!
    call srfflx_state_reset (srfflx_state2d(c))
#endif

   end do   !end loop over chunks
   call t_stopf ('bc_physics')

#ifdef DIUR_AV  !AJF, 9/7/07: Return diurnal-average fields in *_da, pass thru comsrf

   call diur_av(tmpz(:,begchunk:endchunk,:),qrs_da(:,begchunk:endchunk,:),pver)
   call diur_av(fsns,fsns_da,1)
   call diur_av(fsds,fsds_da,1)
   call diur_av(fsnt,fsnt_da,1)
   call diur_av(tmp_sols,sols_da,1)
   call diur_av(tmp_solsd,solsd_da,1)
   do c=begchunk,endchunk
    surface_state2d(c)%sols(:)=sols_da(:,c)
    surface_state2d(c)%solsd(:)=solsd_da(:,c)
   enddo

#endif

#ifdef TRACER_CHECK
   call gavglook ('between DRY', phys_state, gw)
#endif



#if ( ! defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities - no flux coupler
!-----------------------------------------------------------------------
!
   if (.not. aqua_planet) then
!
! Call land model driving routine
!
#ifdef TIMING_BARRIERS
      call t_startf ('sync_tphysbc_lnd')
      call mpibarrier (mpicom)
      call t_stopf ('sync_tphysbc_lnd')
#endif
      call t_startf ('atmlnd_drv')

#if ( defined SCAM )
       if (landfrac(1,begchunk).gt.0) &
#endif

! print out diagnostic info to understand atm--> land coupling...

!!!!      if (masterproc) then
!!!!       write(6,*) 'Just before atmlnd_drv in physpkg'
!!!!       c=begchunk
!!!!       write(6,*) ' THBOT (1st chunk)= ',surface_state2d(c)%thbot
!!!!       write(6,*) ' UBOT (1st chunk)= ',surface_state2d(c)%ubot
!!!!       write(6,*) ' TAUX (1st chunk)= ',srfflx_state2d(c)%wsx
!!!!      endif


      call atmlnd_drv(nstep, iradsw, eccen, obliqr, lambm0,&
                      mvelpp,surface_state2d,srfflx_parm2d)

!!!          print*,' ts,tref after atmlnd_drv '
!!!          do c =begchunk,endchunk
!!!              ncol = get_ncols_p(c)         
!!!                 print*,' chunk= ',c,' ncol= ',ncol
!!!                 do i=1,ncol
!!!                   print*,'i, TS= ',i,srfflx_state2d(c)%ts(i)
!!!                 enddo
!!!                 do i=1,ncol
!!!                   print*,'i, TREF= ',i,srfflx_state2d(c)%tref(i)
!!!                 enddo
!!!          enddo



!!!!      if (masterproc) then
!!!!       write(6,*) 'Just after atmlnd_drv in physpkg'
!!!!       c=begchunk
!!!!       write(6,*) ' THBOT (1st chunk)= ',surface_state2d(c)%thbot
!!!!       write(6,*) ' UBOT (1st chunk)= ',surface_state2d(c)%ubot
!!!!       write(6,*) ' TAUX (1st chunk)= ',srfflx_parm2d(c)%wsx
!!!!       write(6,*) ' SHF (1st chunk)= ',srfflx_parm2d(c)%shf
!!!!       write(6,*) ' LHF (1st chunk)= ',srfflx_parm2d(c)%lhf
!!!!      endif

      call t_stopf ('atmlnd_drv')
#ifdef TIMING_BARRIERS
      call t_startf ('sync_after_lnd')
      call mpibarrier (mpicom)
      call t_stopf ('sync_after_lnd')
#endif


!
! save off albedos and longwave for som offline vars
!
!$OMP PARALLEL DO PRIVATE (C,NCOL,I)

      do c=begchunk,endchunk

         !print*, 'FAO_VVV', c, srfflx_parm2d(c)%tref(1), srfflx_parm2d(c)%ts(1)

         ncol = get_ncols_p(c)
         do i=1,ncol
            if (landfrac(i,c) > 0.) then
               asdirlnd(i,c) = srfflx_parm2d(c)%asdir(i)
               asdiflnd(i,c) = srfflx_parm2d(c)%asdif(i)
               aldirlnd(i,c) = srfflx_parm2d(c)%aldir(i)
               aldiflnd(i,c) = srfflx_parm2d(c)%aldif(i)
               lwuplnd(i,c)  = srfflx_parm2d(c)%lwup(i)
            else
               asdirlnd(i,c) = 0. 
               asdiflnd(i,c) = 0. 
               aldirlnd(i,c) = 0. 
               aldiflnd(i,c) = 0. 
               lwuplnd(i,c)  = 0. 
            end if
         end do
!
!output shf/lhf fluxes for land model
!
         call output_shf_lhf_fluxes(srfflx_parm2d(c), c, ncol, landfrac(1,c), 'LND')
         call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), landfrac(1,c), ncol)

      end do

   end if                    ! end of not aqua_planet if block


#if (defined COUP_SOM)
!
! Set ocean surface quantities - ocn model internal to atm
!
   if (is_end_curr_day ()) then
      call print_coverage ('icefrac', ' million km^2', icefrac, 1.d-12)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = icefrac(i,c)*sicthk(i,c)
         end do
      end do
      call print_coverage ('icevol ', ' 10^13m^3', tmp, 1.d-13)

      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = icefrac(i,c)*snowhice(i,c)
         end do
      end do
      call print_coverage ('snowvol', ' 10^13m^3', tmp, 1.d-13)
   end if

   call t_startf ('somint')
   call somint ()
   call t_stopf ('somint')

   call t_startf ('somoce')
   call somoce (surface_state2d, srfflx_parm2d_ocn)
   call t_stopf ('somoce')

#else

!!!!!!!! TURN OF COUPLING TO OCEAN, ICE "LANDFORMS"
!!!!   call t_startf ('sstint')
!!!!   call sstint ()
!!!!   call t_stopf ('sstint')
!
! iceint may change ocean fraction, so call it before camoce
!
!!!!   call t_startf ('iceint')
!!!!   call iceint ()
!!!!   call t_stopf ('iceint')

!!!!   call t_startf ('camoce')
!!!!   call camoce (surface_state2d, srfflx_parm2d_ocn)
!!!!   call t_stopf ('camoce')
#endif
!
! Set ice surface quantities - icn model internal to atm
!
!!!!!!!! TURN OFF ICE COUPLINGS
!!!!   call t_startf('camice')
!!!!   call camice (surface_state2d, srfflx_parm2d)
!!!!   call t_stopf('camice')
!
! output shf/lhf fluxes for ice/ocn/som_offline 
!
!$OMP PARALLEL DO PRIVATE (C, NCOL, I)
   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
      do i=1,ncol
         if(icefrac(i,c) > 0.) then
            tsice_rad(i,c) = sqrt(sqrt(srfflx_parm2d(c)%lwup(i)/stebol))
         else
            tsice_rad(i,c) = TfrezK
         endif
      end do
      call output_shf_lhf_fluxes (srfflx_parm2d(c), c, ncol, icefrac(1,c), 'ICE')
      call output_shf_lhf_fluxes (srfflx_parm2d_ocn(c), c, ncol, ocnfrac(1,c), 'OCN')
      call output_shfoi_lhfoi_fluxes (srfflx_parm2d_ocn(c), srfflx_parm2d(c), c)

!JR SOM case: Have to wait to call update routine till after both ocean and ice have
!JR operated, since the fractions can change internal to the parameterization
      do i = 1, ncol
         srfflx_state2d(c)%sst(i) = srfflx_parm2d_ocn(c)%ts(i)
      enddo
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d_ocn(c), ocnfrac(1,c), ncol)
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), icefrac(1,c), ncol)
   end do
#endif
!!!!!!!END NOT COUP_CSM BLOCK

#if ( defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities using csm flux coupler
!-----------------------------------------------------------------------
!
! If send data to flux coupler only on radiation time steps:
!
   if (flxave) then
!
! Average the precipitation input to lsm between radiation calls.
!
      call ccsmave(iradsw, nstep, dosw)
!
! Use solar radiation flag to determine data exchange steps 
! with flux coupler. This processes are not independent since 
! instantaneous radiative fluxes are passed, valid over the 
! interval to the next radiation calculation. The same 
! considerations apply to the long and shortwave fluxes, so 
! the intervals must be the same. Data is received from the 
! coupler one step after it is sent.
!
      if (nstep == 0) then
         dorecv = .true.
         dosend = .true.
      else if (nstep == 1) then
         dorecv = .false.
         dosend = .false.
      else if ( (nstep == 2) .and. (iradsw == 1) ) then
         dorecv = .true.
         dosend = dosw
      else
         dorecv = dosend
         dosend = dosw
      end if
   endif
!
! If send data to flux coupler on every time step
!
   if (.not. flxave) then
      if (nstep /= 1) then
         dorecv = .true.
         dosend = .true.
      else 
         dorecv = .false.
         dosend = .false.
      endif
   endif
!
! Send/recv data to/from the csm flux coupler.
!
   if (dosend) call ccsmsnd ( )
   if (dorecv) call ccsmrcv ( )
#endif   
!!!!!!!END COUP_CSM BLOCK
!
!-----------------------------------------------------------------------
! Tendency physics after coupler 
! Not necessary at terminal timestep.
!-----------------------------------------------------------------------
!

   call t_startf ('ac_physics')

!$OMP PARALLEL DO PRIVATE (C, NCOL)

   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
!
! surface diagnostics for history files
!

      call diag_surf (srfflx_state2d(c), surface_state2d(c), icefrac(1,c), ocnfrac(1,c), landfrac(1,c), &
                      sicthk(1,c), snowhland(1,c), snowhice(1,c), tsice(1,c), trefmxav(1,c), &
                      trefmnav(1,c) )

      call t_startf ('tphysac')

      call tphysac (ztodt, pblht(1,c), qpert(1,1,c), tpert(1,c), srfflx_state2d(c)%shf,        &
                    srfflx_state2d(c)%wsx,srfflx_state2d(c)%wsy, srfflx_state2d(c)%cflx, sgh(1,c), srfflx_state2d(c)%lhf,        &
                    landfrac(1,c), snowhland(1,c),srfflx_state2d(c)%tref, surface_state2d(c)%precc, surface_state2d(c)%precl,    &
                    surface_state2d(c)%precsc, surface_state2d(c)%precsl, phys_state(c), phys_tend(c), pbuf, &
                    ocnfrac(1,c), fsds(1,c), icefrac(1,c), fv(1,c), ram1(1,c))
      call t_stopf ('tphysac')
   end do                    ! Chunk loop

   call t_stopf('ac_physics')

#ifdef TRACER_CHECK
   call gavglook ('after tphysac FV:WET)', phys_state, gw )
#endif

   call pbuf_deallocate('physpkg')
   call pbuf_update_tim_idx()

end subroutine physpkg



subroutine gavglook (title, state, gw)

  !
  ! process info from state vectors. Only useful when data in all chunks are in sync
  ! e.g. before and after tphysac and tphysbc
  !

  use physics_types,  only: physics_state, physics_tend
!  use comsrf,         only: srfflx_state
  use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
  use ppgrid,         only: begchunk, endchunk
  use ppgrid,         only: pcols, pver
  use constituents,   only: ppcnst, cnst_name
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use pmgrid,         only: plon, plat, masterproc
  use physconst,      only: gravit
  use time_manager,   only: dtime
#ifdef SPMD
  use mpishorthand
#endif

  implicit none

  ! arguments
  character(len=*), intent(in) :: title
  type(physics_state), intent(in), dimension(begchunk:endchunk) :: state
  real(r8), intent(in) :: gw(plat)                    ! Gaussian weights
  
  ! local
  integer i, lat, c, lon, k
  integer :: lats(pcols)                       ! array of latitude indices
  integer :: lons(pcols)                       ! array of longitude indices
  integer m
  integer :: ncol                              ! number of columns
  real(r8) twodfld(plon,plat,ppcnst)          ! summed at each grid point
  real(r8) twodfle(plon,plat,ppcnst)          ! summed at each grid point
  real(r8) twodflx(plon,plat,ppcnst)          ! summed at each grid point
  real(r8) twodfly(plon,plat,ppcnst)          ! summed at each grid point
#ifdef SPMD                                     
  real(r8) :: twodfld_glob(plon,plat,ppcnst)   ! global summed at each grid point
  real(r8) :: twodfle_glob(plon,plat,ppcnst)   ! global summed at each grid point
  real(r8) :: twodflx_glob(plon,plat,ppcnst)   ! global summed at each grid point
  real(r8) :: twodfly_glob(plon,plat,ppcnst)   ! global summed at each grid point
#endif                                          
  real(r8) :: zonal(plat), zonalw(plat)              !  summed along each latitude
  real(r8) gavg, gavgw
  real(r8) col, wmin, wmax, colw

  !!--------------------------------------------------------


  ! operations on each processor
  twodfld(:,:,:) = 0.
  twodfle(:,:,:) = 0.
  twodflx(:,:,:) = 0.
  twodfly(:,:,:) = 0.
  do c=begchunk, endchunk
     ncol = get_ncols_p(c)
     call get_lat_all_p(c, ncol, lats)
     call get_lon_all_p(c, ncol, lons)
     do m = 1,ppcnst
        do i=1,ncol
           lat = lats(i)
           lon = lons(i)
           col = 0.
           colw = 0.
!           fluxcol = 0.
           wmax = -1.e36
           wmin = 1.e36
           do k = 1,pver
              col  = col + state(c)%pdeldry(i,k)*state(c)%q(i,k,m)*gw(lats(i))
              colw = colw + state(c)%pdel(i,k)  *state(c)%q(i,k,m)*gw(lats(i))
              wmax = max(wmax,state(c)%q(i,k,m))
              wmin = min(wmin,state(c)%q(i,k,m))
           end do ! k
           col = col/gravit
           colw = colw/gravit
           twodfld(lons(i),lats(i),m) = twodfld(lons(i),lats(i),m) + col
           twodfle(lons(i),lats(i),m) = twodfle(lons(i),lats(i),m) + colw
           twodflx(lons(i),lats(i),m) = twodflx(lons(i),lats(i),m) + wmin
           twodfly(lons(i),lats(i),m) = twodfly(lons(i),lats(i),m) + wmax
        enddo ! i
     enddo ! m
  end do ! c

  ! move data to masterproc
#ifdef SPMD
#ifdef TIMING_BARRIERS
  call mpibarrier (mpicom)
#endif
  call mpisum(twodfld, twodfld_glob, plon*plat*ppcnst, mpir8, 0, mpicom)
  call mpisum(twodfle, twodfle_glob, plon*plat*ppcnst, mpir8, 0, mpicom)
  call mpisum(twodflx, twodflx_glob, plon*plat*ppcnst, mpir8, 0, mpicom)
  call mpisum(twodfly, twodfly_glob, plon*plat*ppcnst, mpir8, 0, mpicom)
  if (masterproc) then
     twodfld(:,:,:) = twodfld_glob(:,:,:) 
     twodfle(:,:,:) = twodfle_glob(:,:,:) 
     twodflx(:,:,:) = twodflx_glob(:,:,:) 
     twodfly(:,:,:) = twodfly_glob(:,:,:) 
  endif
#endif

  ! process the data
  if (masterproc) then
     do m = 1,ppcnst
        wmax = -1.e36
        wmin = 1.e36
        do lat=1,plat
           zonal(lat) = 0.
           zonalw(lat) = 0.
           do i=1,plon
              zonal(lat) = zonal(lat) + twodfld(i,lat,m)
              zonalw(lat) = zonalw(lat) + twodfle(i,lat,m)
              wmax = max(wmax,twodfly(i,lat,m))
              wmin = min(wmin,twodflx(i,lat,m))
           end do
        end do
        gavg = 0.
        gavgw = 0.
        do lat=1,plat
           gavg = gavg + zonal(lat)
           gavgw = gavgw + zonalw(lat)
        end do
        gavg = gavg/(2.*plon)
        gavgw = gavgw/(2.*plon)

        write (6,66) trim(title)//' m=',m,'name='//trim(cnst_name(m))//' gavg dry, wet, min, max ' &
             , gavg, gavgw,wmin,wmax
66      format (a24,i2,a36,1p,4g25.14)

     end do
  endif

end subroutine gavglook
