#include <misc.h>

module check_energy

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to check 
!   1. vertically integrated total energy and water conservation for each
!      column within the physical parameterizations
!
!   2. global mean total energy conservation between the physics output state
!      and the input state on the next time step.
!
!   3. add a globally uniform heating term to account for any change of total energy in 2.
!
! Author: Byron Boville  Oct 31, 2002
!         
! Modifications:
!   03.03.29  Boville  Add global energy check and fixer.        
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
  use pmgrid,          only: masterproc
  use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_times
  use phys_gmean,      only: gmean
  use physconst,       only: gravit, latvap, latice
  use physics_types,   only: physics_state, physics_tend, physics_ptend
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind
  use time_manager,    only: is_first_step

  implicit none
  private

! Public types:
  public check_tracers_data

! Public methods
  public :: check_energy_register  ! register fields in physics buffer
  public :: check_energy_init      ! initialization of module
  public :: check_energy_timestep_init  ! timestep initialization of energy integrals and cumulative boundary fluxes
  public :: check_energy_chng      ! check changes in integrals against cumulative boundary fluxes
  public :: check_energy_gmean     ! global means of physics input and output total energy
  public :: check_energy_fix       ! add global mean energy difference as a heating
  public :: check_tracers_init      ! initialize tracer integrals and cumulative boundary fluxes
  public :: check_tracers_chng      ! check changes in integrals against cumulative boundary fluxes


! Private module data

  integer  :: teout_idx            ! teout index in physics buffer

  real(r8) :: teout_glob           ! global mean energy of output state
  real(r8) :: teinp_glob           ! global mean energy of input state
  real(r8) :: tedif_glob           ! global mean energy difference
  real(r8) :: psurf_glob           ! global mean surface pressure
  real(r8) :: ptopb_glob           ! global mean top boundary pressure
  real(r8) :: heat_glob            ! global mean heating rate


  type check_tracers_data
     real(r8) :: tracer(pcols,pcnst)       ! initial vertically integrated total (kinetic + static) energy
     real(r8) :: tracer_tnd(pcols,pcnst)   ! cumulative boundary flux of total energy
     integer :: count(pcnst)               ! count of values with significant imbalances
  end type check_tracers_data

!===============================================================================
contains
!===============================================================================

  subroutine check_energy_register()
!
! Register fields in the physics buffer.
! 
!-----------------------------------------------------------------------
    use phys_buffer,  only: pbuf_times, pbuf_add

    implicit none
!-----------------------------------------------------------------------

! Request physics buffer space for fields that persist across timesteps.
    call pbuf_add('TEOUT', 'global' , 1, 1, pbuf_times, teout_idx)

  end subroutine check_energy_register

!===============================================================================

  subroutine check_energy_init()
!
! Initialize the energy conservation module
! 
!-----------------------------------------------------------------------
    use history,       only: addfld, phys_decomp

    implicit none
!-----------------------------------------------------------------------

! register history variables
    call addfld('TEINP   ', 'W/m2', 1, 'A', 'Total energy of physics input',  phys_decomp)
    call addfld('TEOUT   ', 'W/m2', 1, 'A', 'Total energy of physics output', phys_decomp)
    call addfld('TEFIX   ', 'W/m2', 1, 'A', 'Total energy after fixer',       phys_decomp)

  end subroutine check_energy_init

!===============================================================================

  subroutine check_energy_timestep_init(state, tend, pbuf)

!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state),   intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    type(pbuf_fld),        intent(inout)    :: pbuf(pbuf_size_max)

!---------------------------Local storage-------------------------------

    real(r8) :: ke(pcols)                          ! vertical integral of kinetic energy
    real(r8) :: se(pcols)                          ! vertical integral of static energy
    real(r8) :: wv(pcols)                          ! vertical integral of water (vapor)
    real(r8) :: wl(pcols)                          ! vertical integral of water (liquid)
    real(r8) :: wi(pcols)                          ! vertical integral of water (ice)
    real(r8) :: tr(pcols)                          ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                ! pdel for tracer

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: itim                                ! pbuf time index
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0.
    se = 0.
    wv = 0.
    wl = 0.
    wi = 0.
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
          wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
          wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
       end do
    end do

! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
       state%te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       state%tw_ini(i) = wv(i) + wl(i) + wi(i)

       state%te_cur(i) = state%te_ini(i)
       state%tw_cur(i) = state%tw_ini(i)
    end do

! zero cummulative boundary fluxes 
    tend%te_tnd(:ncol) = 0.
    tend%tw_tnd(:ncol) = 0.

    state%count = 0

! initialize teout in physics buffer
    if (is_first_step()) then
       do itim = 1, pbuf_times
          pbuf(teout_idx)%fld_ptr(1,:ncol,1,lchnk,itim) = state%te_ini(:ncol)
       end do
    end if

  end subroutine check_energy_timestep_init

!===============================================================================
  subroutine check_energy_chng(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)

!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(inout) :: state
    type(physics_tend )    , intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: flx_vap(pcols)      ! boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in   ) :: flx_cnd(pcols)      ! boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in   ) :: flx_ice(pcols)      ! boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in   ) :: flx_sen(pcols)      ! boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8) :: te_xpd(pcols)                   ! expected value (f0 + dt*boundary_flux)
    real(r8) :: te_dif(pcols)                   ! energy of input state - original energy
    real(r8) :: te_tnd(pcols)                   ! tendency from last process
    real(r8) :: te_rer(pcols)                   ! relative error in energy column

    real(r8) :: tw_xpd(pcols)                   ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tw_dif(pcols)                   ! tw_inp - original water
    real(r8) :: tw_tnd(pcols)                   ! tendency from last process
    real(r8) :: tw_rer(pcols)                   ! relative error in water column

    real(r8) :: ke(pcols)                          ! vertical integral of kinetic energy
    real(r8) :: se(pcols)                          ! vertical integral of static energy
    real(r8) :: wv(pcols)                          ! vertical integral of water (vapor)
    real(r8) :: wl(pcols)                          ! vertical integral of water (liquid)
    real(r8) :: wi(pcols)                          ! vertical integral of water (ice)

    real(r8) :: te(pcols)                          ! vertical integral of total energy
    real(r8) :: tw(pcols)                          ! vertical integral of total water

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    logical :: flag
!-----------------------------------------------------------------------
!!$    if (.true.) return

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0.
    se = 0.
    wv = 0.
    wl = 0.
    wi = 0.
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
          wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
          wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
       end do
    end do

! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
       te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       tw(i) = wv(i) + wl(i) + wi(i)
    end do

!!$    print *, "chk_chng", state%lchnk,ncol,name
! compute expected values and tendencies
    do i = 1, ncol
! change in static energy and total water
       te_dif(i) = te(i) - state%te_cur(i)
       tw_dif(i) = tw(i) - state%tw_cur(i)

! expected tendencies from boundary fluxes for last process
       te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000.*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000.

! cummulative tendencies from boundary fluxes
       tend%te_tnd(i) = tend%te_tnd(i) + te_tnd(i)
       tend%tw_tnd(i) = tend%tw_tnd(i) + tw_tnd(i)

! expected new values from previous state plus boundary fluxes
       te_xpd(i) = state%te_cur(i) + te_tnd(i)*ztodt
       tw_xpd(i) = state%tw_cur(i) + tw_tnd(i)*ztodt

! relative error, expected value - input state / previous state 
       te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i)
       tw_rer(i) = (tw_xpd(i) - tw(i)) / state%tw_cur(i)
    end do

#ifndef PERGRO
! error checking
   !if (any(abs(te_rer(1:ncol)) > 1.E-10 .or. abs(tw_rer(1:ncol)) > 1.E-10)) then
   if (any(abs(te_rer(1:ncol)) > 1.E-7 .or. abs(tw_rer(1:ncol)) > 1.E-7)) then
    print *, 'FAO_CHK1: ', te_rer(1:ncol)
    print *, 'FAO_CHK2: ', tw_rer(1:ncol)

       do i = 1, ncol
! the relative error threshold for the water budget has been reduced to 1.e-10
! to avoid messages generated by QNEG3 calls
! PJR- change to identify if error in energy or water 
          if (abs(te_rer(i)) > 1.E-10 ) then 
             state%count = state%count + 1
             write(6,*) "significant conservations error energy after ", name,        &
                " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
             write(6,*) te(i),te_xpd(i),te_dif(i),tend%te_tnd(i)*ztodt,  &
                te_tnd(i)*ztodt,te_rer(i)
          endif
          if ( abs(tw_rer(i)) > 1.E-10) then
             state%count = state%count + 1
             write(6,*) "significant conservations error water after ", name,        &
                " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
             write(6,*) tw(i),tw_xpd(i),tw_dif(i),tend%tw_tnd(i)*ztodt,  &
                tw_tnd(i)*ztodt,tw_rer(i)
          end if
       end do
    end if
#endif

! copy new value to state
    do i = 1, ncol
       state%te_cur(i) = te(i)
       state%tw_cur(i) = tw(i)
    end do

    return
  end subroutine check_energy_chng


!===============================================================================
  subroutine check_energy_gmean(state, pbuf, dtime, nstep)

!-----------------------------------------------------------------------
! Compute global mean total energy of physics input and output states
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
    type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf

    real(r8), intent(in) :: dtime        ! physics time step
    integer , intent(in) :: nstep        ! current timestep number

!---------------------------Local storage-------------------------------
    integer :: ncol                      ! number of active columns
    integer :: itim                      ! time  index in pbuf
    integer :: lchnk                     ! chunk index

    real(r8) :: te(pcols,begchunk:endchunk,3)   
                                         ! total energy of input/output states (copy)
    real(r8) :: te_glob(3)               ! global means of total energy
!-----------------------------------------------------------------------

! Find previous total energy in buffer
    itim = pbuf_old_tim_idx()

! Copy total energy out of input and output states
!DIR$ CONCURRENT
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
! input energy
       te(:ncol,lchnk,1) = state(lchnk)%te_ini(:ncol)
! output energy
       te(:ncol,lchnk,2) = pbuf(teout_idx)%fld_ptr(1,:ncol,1,lchnk,itim)
! surface pressure for heating rate
       te(:ncol,lchnk,3) = state(lchnk)%pint(:ncol,pver+1)
    end do

! Compute global means of input and output energies and of
! surface pressure for heating rate (assume uniform ptop)
    call gmean(te, te_glob, 3)
    teinp_glob = te_glob(1)
!!$    print *, 'teinp_glob',teinp_glob
    teout_glob = te_glob(2)
!!$    print *, 'teout_glob',teout_glob
    psurf_glob = te_glob(3)
    ptopb_glob = state(begchunk)%pint(1,1)
!!$    print *, 'psurf_glob',psurf_glob
!!$    print *, 'ptopb_glob',ptopb_glob

! Global mean total energy difference
    tedif_glob =  teinp_glob - teout_glob
    heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)
!!$    print *, 'tedif_glob',tedif_glob
!!$    print *, 'heat_glob',heat_glob

#ifdef VERBOSE
    if (masterproc) then
       print *   , "nstep, te", nstep, teinp_glob, tedif_glob/dtime, heat_glob, psurf_glob
    end if
#endif
    
  end subroutine check_energy_gmean

!===============================================================================
  subroutine check_energy_fix(state, ptend, nstep, eshflx)

!-----------------------------------------------------------------------
! Add heating rate required for global mean total energy conservation
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ) :: state
    type(physics_ptend), intent(inout) :: ptend

    integer , intent(in   ) :: nstep          ! time step number
    real(r8), intent(out  ) :: eshflx(pcols)  ! effective sensible heat flux

!---------------------------Local storage-------------------------------
    integer  :: i,k                      ! column, level indexes
    integer  :: ncol                     ! number of atmospheric columns in chunk
!-----------------------------------------------------------------------
    ncol = state%ncol
    ptend%name         = 'chkenergyfix'
    ptend%ls           = .TRUE.

! add (-) global mean total energy difference as heating
    ptend%s(:ncol,:pver) = heat_glob
!!$    print *, "chk_fix: heat", state%lchnk, ncol, heat_glob

! compute effective sensible heat flux
    do i = 1, ncol
       eshflx(i) = heat_glob * (state%pint(i,pver+1) - state%pint(i,1)) / gravit
    end do
!!!    if (nstep > 0) print *, "heat", heat_glob, eshflx(1)

    return
  end subroutine check_energy_fix


!===============================================================================
  subroutine check_tracers_init(state, tracerint)

!-----------------------------------------------------------------------
! Compute initial values of tracers integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    type(physics_state),   intent(in)    :: state
    type(check_tracers_data), intent(out)   :: tracerint

!---------------------------Local storage-------------------------------

    real(r8) :: tr(pcols)                          ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                ! pdel for tracer

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k,m                                 ! column, level,constituent indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ and tracer indices

!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    do m = 1,pcnst

       if (m.eq.1.or.m.eq.ixcldice.or.m.eq.ixcldliq) exit   ! dont process water substances
                                                            ! they are checked in check_energy
       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals of tracer
       tr = 0.
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
          end do
       end do

       ! Compute vertical integrals of frozen static tracers and total water.
       do i = 1, ncol
          tracerint%tracer(i,m) = tr(i)
       end do

       ! zero cummulative boundary fluxes 
       tracerint%tracer_tnd(:ncol,m) = 0.

       tracerint%count(m) = 0

    end do

    return
  end subroutine check_tracers_init

!===============================================================================
  subroutine check_tracers_chng(state, tracerint, name, nstep, ztodt, cflx)

!-----------------------------------------------------------------------
! Check that the tracers and water change matches the boundary fluxes
! these checks are not save when there are tracers transformations, as 
! they only check to see whether a mass change in the column is
! associated with a flux
!-----------------------------------------------------------------------

    use abortutils, only: endrun 


    implicit none

!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(in   ) :: state
    type(check_tracers_data), intent(inout) :: tracerint! tracers integrals and boundary fluxes
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: cflx(pcols,pcnst)       ! boundary flux of tracers       (kg/m2/s)

!---------------------------Local storage-------------------------------

    real(r8) :: tracer_inp(pcols,pcnst)                   ! total tracer of new (input) state
    real(r8) :: tracer_xpd(pcols,pcnst)                   ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tracer_dif(pcols,pcnst)                   ! tracer_inp - original tracer
    real(r8) :: tracer_tnd(pcols,pcnst)                   ! tendency from last process
    real(r8) :: tracer_rer(pcols,pcnst)                   ! relative error in tracer column

    real(r8) :: tr(pcols)                           ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                       ! pdel for tracer

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: m                            ! tracer index
    character(len=8) :: tracname   ! tracername
!-----------------------------------------------------------------------
!!$    if (.true.) return

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    do m = 1,pcnst
       
       if (m.eq.1.or.m.eq.ixcldice.or.m.eq.ixcldliq) exit   ! dont process water substances
                                                            ! they are checked in check_energy

       tracname = cnst_name(m)
       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals tracers
       tr = 0.
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
          end do
       end do

       ! Compute vertical integrals of tracer
       do i = 1, ncol
          tracer_inp(i,m) = tr(i)
       end do

       ! compute expected values and tendencies
       do i = 1, ncol
          ! change in tracers 
          tracer_dif(i,m) = tracer_inp(i,m) - tracerint%tracer(i,m)

          ! expected tendencies from boundary fluxes for last process
          tracer_tnd(i,m) = cflx(i,m)

          ! cummulative tendencies from boundary fluxes
          tracerint%tracer_tnd(i,m) = tracerint%tracer_tnd(i,m) + tracer_tnd(i,m)

          ! expected new values from original values plus boundary fluxes
          tracer_xpd(i,m) = tracerint%tracer(i,m) + tracerint%tracer_tnd(i,m)*ztodt

          ! relative error, expected value - input value / original 
          tracer_rer(i,m) = (tracer_xpd(i,m) - tracer_inp(i,m)) / tracerint%tracer(i,m)
       end do

!! final loop for error checking
!    do i = 1, ncol

!! error messages
!       if (abs(enrgy_rer(i)) > 1.E-14 .or. abs(water_rer(i)) > 1.E-14) then
!          tracerint%count = tracerint%count + 1
!          write(6,*) "significant conservations error after ", name,        &
!               " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col", i
!          write(6,*) enrgy_inp(i),enrgy_xpd(i),enrgy_dif(i),tracerint%enrgy_tnd(i)*ztodt,  &
!               enrgy_tnd(i)*ztodt,enrgy_rer(i)
!          write(6,*) water_inp(i),water_xpd(i),water_dif(i),tracerint%water_tnd(i)*ztodt,  &
!               water_tnd(i)*ztodt,water_rer(i)
!       end if
!    end do


       ! final loop for error checking
       if ( maxval(tracer_rer) > 1.E-14 ) then
          write(6,*) "CHECK_TRACERS TRACER large rel error"
          write(6,*) tracer_rer
       endif

       do i = 1, ncol
          ! error messages
          if (abs(tracer_rer(i,m)) > 1.E-14 ) then
             tracerint%count = tracerint%count + 1
             write(6,*) "CHECK_TRACERS TRACER significant conservation error after ", name,        &
                  " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col",i
             write(6,*)' process name, tracname, index ',  name, tracname, m
             write(6,*)" input integral              ",tracer_inp(i,m)
             write(6,*)" expected integral           ", tracer_xpd(i,m)
             write(6,*)" input - inital integral     ",tracer_dif(i,m)
             write(6,*)" cumulative tend      ",tracerint%tracer_tnd(i,m)*ztodt
             write(6,*)" process tend         ",tracer_tnd(i,m)*ztodt
             write(6,*)" relative error       ",tracer_rer(i,m)
             call endrun()
          end if
       end do
    end do

    return
  end subroutine check_tracers_chng


end module check_energy
