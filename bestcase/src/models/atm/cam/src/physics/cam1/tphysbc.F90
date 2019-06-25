#include <misc.h>
#include <params.h>


subroutine tphysbc                                               &
                   (ztodt,   pblht,   tpert,   ts,      sst,     &
                    qpert,   precl,   precc,   precsl,  precsc,  &
                    asdir,   asdif,   aldir,   aldif,   snowh,   &
                    qrs,     qrl,     flwds,   fsns,    fsnt,    &
                    flns,    flnt,    lwup,    srfrad,  sols,    &
                    soll,    solsd,   solld,   state,   tend,    &
                    pbuf,    prcsnw,  fsds ,   landm,   landfrac,&
		    ocnfrac, icefrac,                            &
                    tref,    taux,    tauy)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Tendency physics BEFORE coupling to land, sea, and ice models.
! 
! Method: 
! Call physics subroutines and compute the following:
!     o cloud calculations (cloud fraction, emissivity, etc.)
!     o radiation calculations
! Pass surface fields for separate surface flux calculations
! Dump appropriate fields to history file.
! 
! Author: CCM1, CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p, get_lon_all_p, get_lat_all_p
   use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use cldcond,         only: cldcond_tend, cldcond_zmconv_detrain, cldcond_sediment
   use param_cldoptics, only: param_cldoptics_calc
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use diagnostics,     only: diag_dynvar
   use history,         only: outfld
   use physconst,       only: gravit, latvap, cpair, tmelt, cappa, zvir, rair, rga
   use radheat,         only: radheat_net
   use constituents,    only: pcnst, pnats, ppcnst, qmin
   use constituents,    only: dcconnam, cnst_get_ind
   use zm_conv,         only: zm_conv_evap, zm_convr
   use time_manager,    only: is_first_step, get_nstep, get_curr_calday, get_curr_titan_calday
   use moistconvection, only: cmfmca
   use check_energy,    only: check_energy_chng, check_energy_fix
   use check_energy,    only: check_tracers_data, check_tracers_init, check_tracers_chng
   use dycore,          only: dycore_is
   use cloudsimulator,  only: doisccp, cloudsimulator_run
   use aerosol_intr, only: aerosol_wet_intr
   use carma,		only: carma_is_active, carma_wetdep_tend

   ! FAO
   use geopotential , only: geopotential_t
   use pmgrid,        only:  masterproc   !AJF 6/9/06
   use abortutils,   only: endrun
! AJF,9/7/07 :  Add comsrf to access diurnally-averaged sw fields
   use comsrf,       only: qrs_da,fsds_da,fsns_da,fsnt_da,sols_da,solsd_da

   implicit none

#include <comctl.h>
!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: ts(pcols)                      ! surface temperature
   real(r8), intent(in) :: sst(pcols)                     ! sea surface temperature
   real(r8), intent(inout) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(inout) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(inout) :: qpert(pcols,ppcnst)         ! Thermal humidity & constituent excess
   real(r8), intent(in) :: asdir(pcols)                  ! Albedo: shortwave, direct
   real(r8), intent(in) :: asdif(pcols)                  ! Albedo: shortwave, diffuse
   real(r8), intent(in) :: aldir(pcols)                  ! Albedo: longwave, direct
   real(r8), intent(in) :: aldif(pcols)                  ! Albedo: longwave, diffuse
   real(r8), intent(in) :: snowh(pcols)                  ! Snow depth (liquid water equivalent)
   real(r8), intent(inout) :: qrs(pcols,pver)            ! Shortwave heating rate
   real(r8), intent(inout) :: qrl(pcols,pver)            ! Longwave  heating rate
   real(r8), intent(inout) :: flwds(pcols)               ! Surface longwave down flux
   real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: lwup(pcols)                    ! Surface longwave up flux
   real(r8), intent(out) :: srfrad(pcols)                 ! Net surface radiative flux (watts/m**2)
   real(r8), intent(inout) :: sols(pcols)                   ! Direct beam solar rad. onto srf (sw)
   real(r8), intent(inout) :: soll(pcols)                   ! Direct beam solar rad. onto srf (lw)
   real(r8), intent(inout) :: solsd(pcols)                  ! Diffuse solar radiation onto srf (sw)
   real(r8), intent(inout) :: solld(pcols)                  ! Diffuse solar radiation onto srf (lw)
   real(r8), intent(out) :: precl(pcols)                  ! Large-scale precipitation rate
   real(r8), intent(out) :: precc(pcols)                  ! Convective-scale preciptn rate
   real(r8), intent(out) :: precsl(pcols)                 ! L.S. snowfall rate
   real(r8), intent(out) :: precsc(pcols)                 ! C.S. snowfall rate
   real(r8), intent(out) :: prcsnw(pcols)                 ! snowfall rate (precsl + precsc)
   real(r8), intent(out) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
   real(r8), intent(in) :: landfrac(pcols)                ! land fraction
   real(r8), intent(in) :: ocnfrac(pcols)                ! land fraction
   real(r8), intent(in) :: icefrac(pcols)                ! land fraction

   ! FAO
   real(r8), intent(in)  :: tref(pcols)            ! surface reference temperature
   real(r8), intent(out) :: taux(pcols)            ! X surface stress (zonal)
   real(r8), intent(out) :: tauy(pcols)            ! Y surface stress (meridional)

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: rhdfda(pcols,pver)            ! dRh/dcloud, old 
   real(r8) :: rhu00 (pcols,pver)            ! Rh threshold for cloud, old

   type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies

   integer :: nstep                          ! current timestep number
   integer      lat(pcols)                   ! current latitudes(indices)
   integer      lon(pcols)                   ! current longtitudes(indices)

   real(r8) :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)

   real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: zmrprd(pcols,pver)            ! rain production in ZM convection
   real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c
   real(r8) :: cmfsl(pcols,pver)             ! Moist convection lw stat energy flux
   real(r8) :: cmflq(pcols,pver)             ! Moist convection total water flux
   real(r8) :: dtcond(pcols,pver)            ! dT/dt due to moist processes
   real(r8) :: dqcond(pcols,pver,ppcnst)     ! dq/dt due to moist processes

   real(r8) cldst(pcols,pver)
   real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
   real(r8) cllow(pcols)                      !       "     low  cloud cover
   real(r8) clmed(pcols)                      !       "     mid  cloud cover
   real(r8) clhgh(pcols)                      !       "     hgh  cloud cover
   real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
   real(r8) cmfdqr2(pcols,pver)               ! dq/dt due to moist convective rainout
   real(r8) cmfmc2(pcols,pverp)               ! Moist convection cloud mass flux
   real(r8) cmfsl2(pcols,pver)                ! Moist convection lw stat energy flux
   real(r8) cmflq2(pcols,pver)                ! Moist convection total water flux
   real(r8) cnt(pcols)                        ! Top level of convective activity
   real(r8) cnb(pcols)                        ! Lowest level of convective activity
   real(r8) cnt2(pcols)                       ! Top level of convective activity
   real(r8) cnb2(pcols)                       ! Bottom level of convective activity
   real(r8) concld(pcols,pver)             
   real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
   real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from convection
   real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
   real(r8) prect(pcols)                      ! total (conv+large scale) precip rate
   real(r8) dlf2(pcols,pver)                   ! dq/dt due to rainout terms
   real(r8) qpert2(pcols,ppcnst)              ! Perturbation q
   real(r8) rtdt                              ! 1./ztodt
   real(r8) tpert2(pcols)                     ! Perturbation T
   real(r8) pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
!                                             !    maximally overlapped region.
!                                             !    0->pmxrgn(i,1) is range of pressure for
!                                             !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                             !    2nd region, etc
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   integer nmxrgn(pcols)                      ! Number of maximally overlapped regions
   integer  i,k,m                             ! Longitude, level, constituent indices
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
                                           
!  real(r8) engt                              ! Thermal   energy integral
!  real(r8) engk                              ! Kinetic   energy integral
!  real(r8) engp                              ! Potential energy integral
   real(r8) rel(pcols,pver)                   ! Liquid cloud particle effective radius
   real(r8) rei(pcols,pver)                   ! Ice effective drop size (microns)
   real(r8) emis(pcols,pver)                  ! Cloud longwave emissivity
   real(r8) clc(pcols)                        ! Total convective cloud (cloud scheme)
   real(r8) :: cicewp(pcols,pver)             ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,pver)             ! in-cloud cloud liquid water path
!
   real(r8) dellow(pcols)                     ! delta p for bottom three levels of model
   real(r8) tavg(pcols)                       ! mass weighted average temperature for 

! physics buffer fields to compute tendencies for cloud condensation package
   integer itim, ifld
   real(r8), pointer, dimension(:,:) :: qcwat, tcwat, lcwat, cld

! physics buffer fields for total energy and mass adjustment
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: tini
!                                          
! Used for OUTFLD only                     
!                                          
   real(r8) icwmr1(pcols,pver)                ! in cloud water mixing ration for zhang scheme
   real(r8) icwmr2(pcols,pver)                ! in cloud water mixing ration for hack scheme
   real(r8) fracis(pcols,pver,ppcnst)         ! fraction of transported species that are insoluble
   real(r8) timestep(pcols)
!
!     Variables for doing deep convective transport outside of zm_convr
!
   real(r8) mu2(pcols,pver)
   real(r8) eu2(pcols,pver)
   real(r8) du2(pcols,pver)
   real(r8) md2(pcols,pver)
   real(r8) ed2(pcols,pver)
   real(r8) dp(pcols,pver)
   real(r8) dpdry(pcols,pver)
   real(r8) dsubcld(pcols)
   real(r8) conicw(pcols,pver)
   real(r8) cmfdqrt(pcols,pver)               ! dq/dt due to moist convective rainout

! stratiform precipitation variables
   real(r8) :: prec_pcw(pcols)                ! total precip from prognostic cloud scheme
   real(r8) :: snow_pcw(pcols)                ! snow from prognostic cloud scheme
   real(r8) :: prec_sed(pcols)                ! total precip from cloud sedimentation
   real(r8) :: snow_sed(pcols)                ! snow from cloud ice sedimentation

! convective precipitation variables
   real(r8) :: prec_zmc(pcols)                ! total precipitation from ZM convection
   real(r8) :: snow_zmc(pcols)                ! snow from ZM convection
   real(r8) :: prec_cmf(pcols)                ! total precipitation from Hack convection
   real(r8) :: snow_cmf(pcols)                ! snow from Hack convection

! energy checking variables
   real(r8) :: zero(pcols)                    ! array of zeros
   real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
   real(r8) :: flx_cnd(pcols)
   real(r8) :: flx_heat(pcols)
   logical  :: conserve_energy = .true.       ! flag to carry (QRS,QRL)*dp across time steps
   type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes

   integer jt(pcols)
   integer maxg(pcols)
   integer ideep(pcols)
   integer lengath
   real(r8) cldc(pcols,pver)
   real(r8) nevapr(pcols,pver)
   real(r8) qme(pcols,pver)
   real(r8) prain(pcols,pver)
   real(r8) cflx(pcols,ppcnst)

   ! FAO
   real(r8) frac_day, day_in_year
   real(r8) cld_zero(pcols,pver)   ! 2D array of zeros to disable clouds
   ! AJF
   real(r8) tmp(pcols)             ! workspace for outfld (diurnal averaging)


!-----------------------------------------------------------------------

   ! FAO: initialize -- probably unnecessary
   precl = 0.0
   precc = 0.0
   precsl = 0.0
   precsc = 0.0
   srfrad = 0.0
   prcsnw = 0.0

   zero = 0.
   cld_zero = 0.     ! FAO -- kludge to disable clouds

   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1./ztodt

   nstep = get_nstep()
   calday = get_curr_calday()
   call get_curr_titan_calday(frac_day, day_in_year)

!
! Output NSTEP for debugging
!
   timestep(:ncol) = nstep
   call outfld ('NSTEP   ',timestep, pcols, lchnk)

!  dry surface pressure
   call outfld ('PSDRY',  state%psdry, pcols, lchnk)

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('QCWAT')
   qcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('TCWAT')
   tcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('LCWAT')
   lcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('TEOUT')
   teout => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,itim)
   ifld = pbuf_get_fld_idx('QINI')
   qini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('TINI')
   tini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

!
! Set physics tendencies to 0
   tend %dTdt(:ncol,:pver)  = 0.
   tend %dudt(:ncol,:pver)  = 0.
   tend %dvdt(:ncol,:pver)  = 0.
     do i=1,ncol
        tend%flx_net(i) = 0.
     end do

   call physics_ptend_init (ptend) ! Initialize parameterization tendency structure

!
! Make sure that input tracers are all positive (probably unnecessary)
!
   call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
              ppcnst,qmin  ,state%q )
!
! Setup q and t accumulation fields
!
   dqcond(:ncol,:,:) = state%q(:ncol,:,:)
   dtcond(:ncol,:)   = state%s(:ncol,:)

   fracis (:ncol,:,1:ppcnst) = 1.

! compute mass integrals of input tracers state
   call check_tracers_init(state, tracerint)

!===================================================
! Global mean total energy fixer
!===================================================
   !*** BAB's FV heating kludge *** save the initial temperature
   tini(:ncol,:pver) = state%t(:ncol,:pver)
   if (dycore_is('LR')) then
      call check_energy_fix(state, ptend, nstep, flx_heat)
      call physics_update(state, tend, ptend, ztodt)
      call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
   end if
   qini(:ncol,:pver) = state%q(:ncol,:pver,1)

   call outfld('TEOUT', teout       , pcols, lchnk   )
   call outfld('TEINP', state%te_ini, pcols, lchnk   )
   call outfld('TEFIX', state%te_cur, pcols, lchnk   )
!
!===================================================
! Dry adjustment
!===================================================

! Copy state info for input to dadadj
! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

   ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
   ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)


   call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
                ptend%s, ptend%q(1,1,1))
   ptend%name  = 'dadadj'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
   ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
   call physics_update (state, tend, ptend, ztodt)


!!! skip moist convection calc when TITAN_FAO turned on, ajf 2/24/06
#if (! defined TITAN_FAO)
!
!===================================================
! Moist convection
!===================================================
!
! Since the PBL doesn't pass constituent perturbations, they
! are zeroed here for input to the moist convection routine
!
   qpert(:ncol,2:ppcnst) = 0.0
!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr')
   call zm_convr( lchnk,    ncol, &
                  state%t,   state%q,    prec_zmc,   cnt,     cnb,      &
                  pblht,   state%zm, state%phis,    state%zi,   ptend%q(:,:,1),     &
                  ptend%s, state%pmid,   state%pint,  state%pdel,       &
                  .5*ztodt,cmfmc,    cmfcme,             &
                  tpert,   dlf,      pflx,    zdu,     zmrprd,   &
                  mu2,      md2,     du2,     eu2,     ed2,      &
                  dp,       dsubcld, jt,      maxg,    ideep,    &
                  lengath, icwmr1,   rliq    )
   ptend%name  = 'zm_convr'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   cmfsl (:ncol,:) = 0. ! This is not returned from zm, hence it is zeroed.
   cmflq (:ncol,:) = 0. ! This is not returned from zm, hence it is zeroed.

   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend%q(1,1,1) ,pcols   ,lchnk   )
   call t_stopf('zm_convr')

   call physics_update(state, tend, ptend, ztodt)
!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
   call zm_conv_evap(state, ptend, zmrprd, cld, ztodt, prec_zmc, snow_zmc, .false.)
   call physics_update(state, tend, ptend, ztodt)
! Check energy integrals, including "reserved liquid"
   flx_cnd(:ncol) = prec_zmc(:ncol) + rliq(:ncol)
   call check_energy_chng(state, tend, "zm_evap", nstep, ztodt, zero, flx_cnd, snow_zmc, zero)

! Transport cloud water and ice only
!
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   ptend%name = 'convtran1'
   ptend%lq(ixcldice) = .true.
   ptend%lq(ixcldliq) = .true.
   call convtran (lchnk,                                        &
                  ptend%lq,state%q, ppcnst,  mu2,     md2,   &
                  du2,     eu2,     ed2,     dp,      dsubcld,  &
                  jt,      maxg,    ideep,   1,       lengath,  &
                  nstep,   fracis,  ptend%q, dpdry  )
   call physics_update (state, tend, ptend, ztodt)
!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   cmfmc(:ncol,:pver) = cmfmc(:ncol,:pver) * 100./gravit
!
! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
!
   call t_startf('cmfmca')
   tpert2(:ncol  ) =0.
   qpert2(:ncol,:) = qpert(:ncol,:)  ! BAB Why is this not zero, if tpert2=0???
   call cmfmca (lchnk,   ncol, &
                nstep,   ztodt,   state%pmid,  state%pdel,   &
                state%rpdel,   state%zm,      tpert2,  qpert2,  state%phis,     &
                pblht,   state%t,   state%q,   ptend%s,   ptend%q,      &
                cmfmc2,  cmfdqr2, cmfsl2,  cmflq2,  prec_cmf,   &
                dlf2,     cnt2,    cnb2,    icwmr2   , rliq2, & 
                state%pmiddry, state%pdeldry, state%rpdeldry)
   ptend%name  = 'cmfmca'
   ptend%ls    = .TRUE.
   ptend%lq(:) = .TRUE.

! Add shallow cloud water detrainment to cloud water detrained from ZM
   dlf(:ncol,:pver) = dlf(:ncol,:pver) + dlf2(:ncol,:pver)
   rliq(:ncol) = rliq(:ncol) + rliq2(:ncol)
   
   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('CMFDT   ',ftem          ,pcols   ,lchnk   )
   call outfld('CMFDQ   ',ptend%q(1,1,1),pcols   ,lchnk   )
   call t_stopf('cmfmca')
   call physics_update (state, tend, ptend, ztodt)
!
! Determine the phase of the precipitation produced and add latent heat of fusion
   call zm_conv_evap(state, ptend, cmfdqr2, cld, ztodt, prec_cmf, snow_cmf, .true.)
   call physics_update(state, tend, ptend, ztodt)
   flx_cnd(:ncol) = prec_cmf(:ncol) + rliq2(:ncol)
   call check_energy_chng(state, tend, "hk_evap", nstep, ztodt, zero, flx_cnd, snow_cmf, zero)
!
! Merge shallow/mid-level output with prior results from Zhang-McFarlane
!
   do i=1,ncol
      if (cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
      if (cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
   end do
!
   cmfmc (:ncol,:pver) = cmfmc (:ncol,:pver) + cmfmc2 (:ncol,:pver)
   cmfsl (:ncol,:pver) = cmfsl (:ncol,:pver) + cmfsl2 (:ncol,:pver)
   cmflq (:ncol,:pver) = cmflq (:ncol,:pver) + cmflq2 (:ncol,:pver)
   call outfld('CMFMC' , cmfmc  , pcols, lchnk)
!  output new partition of cloud condensate variables, as well as precipitation 
   call outfld('QC      ',dlf2           ,pcols   ,lchnk   )
   call outfld('PRECSH  ',prec_cmf      ,pcols   ,lchnk       )
   call outfld('CMFDQR', cmfdqr2, pcols, lchnk)
   call outfld('CMFSL' , cmfsl  , pcols, lchnk)
   call outfld('CMFLQ' , cmflq  , pcols, lchnk)
   call outfld('DQP'   , dlf2    , pcols, lchnk)

! Allow the cloud liquid drops and ice particles to sediment
! Occurs before adding convectively detrained cloud water, because the phase of the
! of the detrained water is unknown.
   call t_startf('cldwat_sediment')
   call cldcond_sediment(state, ptend, ztodt,cld, icefrac, landfrac, ocnfrac, prec_sed, &
                         snow_sed, landm, snowh)
   call physics_update(state, tend, ptend, ztodt)
   call t_stopf('cldwat_sediment')

! check energy integrals
   call check_energy_chng(state, tend, "cldwat_sediment", nstep, ztodt, zero, prec_sed, snow_sed, zero)

! Put the detraining cloud water from convection into the cloud and environment. 

   call cldcond_zmconv_detrain(dlf, cld, state, ptend)
   call physics_update(state, tend, ptend, ztodt)

! check energy integrals, reserved liquid has now been used
   flx_cnd(:ncol) = -rliq(:ncol)
   call check_energy_chng(state, tend, "cldwat_detrain", nstep, ztodt, zero, flx_cnd, zero, zero)
!
! cloud fraction after transport and convection,
! derive the relationship between rh and cld from 
! the employed cloud scheme
!
   call cldnrh(lchnk,   ncol,                                &
               state%pmid,    state%t,   state%q(1,1,1),   state%omega, &
               cnt,     cnb,     cld,    clc,     state%pdel,   &
               cmfmc,   cmfmc2, landfrac,snowh,   concld,  cldst,    &
               ts,      sst, state%pint(1,pverp),       zdu,  ocnfrac, &
               rhdfda,   rhu00 , state%phis)
   call outfld('CONCLD  ',concld, pcols,lchnk)
   call outfld('CLDST   ',cldst,  pcols,lchnk)
   call outfld('CNVCLD  ',clc,    pcols,lchnk)

! cloud water and ice parameterizations
   call t_startf('cldwat_tend')
   call cldcond_tend(state, ptend, ztodt, &
       tcwat, qcwat, lcwat, prec_pcw, snow_pcw, icefrac, rhdfda, rhu00, cld, nevapr, prain, qme, snowh)
   call physics_update (state, tend, ptend, ztodt)
   call t_stopf('cldwat_tend')

! check energy integrals
   call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_pcw, snow_pcw, zero)

! Save off q and t after cloud water
   do k=1,pver
      qcwat(:ncol,k) = state%q(:ncol,k,1)
      tcwat(:ncol,k) = state%t(:ncol,k)
      lcwat(:ncol,k) = state%q(:ncol,k,ixcldice) + state%q(:ncol,k,ixcldliq)
   end do
!
!  aerosol wet chemistry determines scavenging fractions, and transformations
   call get_lat_all_p(lchnk, ncol, lat)
   call get_lon_all_p(lchnk, ncol, lon)
   call get_rlat_all_p(lchnk, ncol, clat)
   conicw(:ncol,:) = icwmr1(:ncol,:) + icwmr2(:ncol,:)
   cmfdqrt(:ncol,:) = zmrprd(:ncol,:) + cmfdqr2(:ncol,:)
   call aerosol_wet_intr (state, ptend, cflx, nstep, ztodt, lat, clat, lon,&
        qme, prain, &
       nevapr, cldc, cld, fracis, calday, cmfdqrt, conicw)
   call physics_update (state, tend, ptend, ztodt)

! Wet scavenging for carma particles (and maybe gases).
   if (carma_is_active()) then
     call t_startf('carma_wetdep_tend')
     call carma_wetdep_tend(state, ptend, ztodt, pbuf)
     call t_stopf('carma_wetdep_tend')
     call physics_update (state, tend, ptend, ztodt)
   end if

!
!     Convective transport of all trace species except water vapor and
!     cloud liquid and ice done here because we need to do the scavenging first
!     to determine the interstitial fraction.
!
   ptend%name  = 'convtran2'
   ptend%lq(:) = .true.
   ptend%lq(ixcldice) = .false.
   ptend%lq(ixcldliq) = .false.
   dpdry = 0
   do i = 1,lengath
      dpdry(i,:) = state%pdeldry(ideep(i),:)/100.
   end do
   call convtran (lchnk,                                           &
                  ptend%lq,state%q, ppcnst,     mu2,     md2,      &
                  du2,     eu2,     ed2,        dp,      dsubcld,  &
                  jt,      maxg,    ideep,      1,       lengath,  &
                  nstep,   fracis,  ptend%q,    dpdry)

   call physics_update (state, tend, ptend, ztodt)

! check tracer integrals
   call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt, cflx)

!
! Compute rates of temperature and constituent change due to moist processes
!
   dtcond(:ncol,:) = (state%s(:ncol,:) - dtcond(:ncol,:))*rtdt / cpair
   dqcond(:ncol,:,:) = (state%q(:ncol,:,:) - dqcond(:ncol,:,:))*rtdt
   call outfld('DTCOND  ',dtcond, pcols   ,lchnk   )
   do m=1,ppcnst
      call outfld(dcconnam(m),dqcond(1,1,m),pcols   ,lchnk )
   end do

! Compute total convective and stratiform precipitation and snow rates
   do i=1,ncol
      precc (i) = prec_zmc(i) + prec_cmf(i)
      precl (i) = prec_sed(i) + prec_pcw(i)
      precsc(i) = snow_zmc(i) + snow_cmf(i)
      precsl(i) = snow_sed(i) + snow_pcw(i)
! jrm These checks should not be necessary if they exist in the parameterizations
      if(precc(i).lt.0.) precc(i)=0.
      if(precl(i).lt.0.) precl(i)=0.
      if(precsc(i).lt.0.) precsc(i)=0.
      if(precsl(i).lt.0.) precsl(i)=0.
      if(precsc(i).gt.precc(i)) precsc(i)=precc(i)
      if(precsl(i).gt.precl(i)) precsl(i)=precl(i)
! end jrm
   end do
   prcsnw(:ncol) = precsc(:ncol) + precsl(:ncol)   ! total snowfall rate: needed by slab ocean model
!
!===================================================
! Moist physical parameteriztions complete: 
! send dynamical variables, and derived variables to history file
!===================================================
!
!!!turn routine back on
#endif

   call diag_dynvar (lchnk, ncol, state)

! ============================================================================
! AJF:  Following involves inlining of tides in same manner as other phys 
!       routines.  Calculates tidal accelerations
!   -ajf 7-13-07
! ============================================================================
!
!
!EJL 2-7-13 adding EQTORQ flag
#ifdef EQTORQ
!!     write(6,*) 'defined EQTORQ'
     call t_startf('tides_titan_eqtorq')
      call tides_titan_eqtorq (ztodt, taux, tauy, state, tend, ptend, flx_heat)
     call t_stopf('tides_titan_eqtorq')

     call physics_update(state,tend,ptend,ztodt)

     call check_energy_chng(state, tend,"tides_titan_eqtorq", nstep, ztodt,  &
       zero, zero, zero, tend%flx_net)
#else
!!     write(6,*) 'EQTORQ not defined'
     call t_startf('tides_titan')
      call tides_titan (ztodt, taux, tauy, state, tend, ptend, flx_heat)
     call t_stopf('tides_titan')

     call physics_update(state,tend,ptend,ztodt)

     call check_energy_chng(state, tend,"tides_titan", nstep, ztodt,  &
       zero, zero, zero, tend%flx_net)

#endif

! ============================================================================

!
!===================================================
! Radiation computations
!===================================================
!
! Cosine solar zenith angle for current time step
!
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)

!#ifdef TITAN_FAO
   call zenith_titan (frac_day, day_in_year, clat, clon, coszrs, ncol)
!#else
!   call zenith (calday, clat, clon, coszrs, ncol)
!#endif

   if (dosw .or. dolw) then

! Compute cloud water/ice paths and optical properties for input to radiation
      call t_startf('cldoptics')
      call param_cldoptics_calc(state, cld, landfrac, landm,icefrac, &
                                cicewp, cliqwp, emis, rel, rei, pmxrgn, nmxrgn, snowh)
      call t_stopf('cldoptics')

!
! Complete radiation calculations
!

      call t_startf ('radctl')
      call radctl (lchnk, ncol, lwup, emis, state%pmid,             &
                   state%pint, state%lnpmid, state%lnpint, state%t, state%q,   &
                   cld_zero, cicewp, cliqwp, coszrs, asdir, asdif,               &
                   aldir, aldif, pmxrgn, nmxrgn, fsns, fsnt    ,flns    ,flnt    , &
                   qrs, qrl, flwds, rel, rei,                       &
                   sols, soll, solsd, solld,                  &
                   landfrac, state%zm, state, fsds, ts, tref)
      call t_stopf ('radctl')

!
! Cloud cover diagnostics
! radctl can change pmxrgn and nmxrgn so cldsav needs to follow 
! radctl.
!
      call cldsav (lchnk, ncol, cld, state%pmid, cltot, &
                   cllow, clmed, clhgh, nmxrgn, pmxrgn)
!
! Dump cloud field information to history tape buffer (diagnostics)
!
      call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
      call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
      call outfld('CLDMED  ',clmed  ,pcols,lchnk)
      call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
      call outfld('CLOUD   ',cld    ,pcols,lchnk)
      if (doisccp) then
         call cloudsimulator_run(state, ts, concld, cld, cliqwp, &
                                 cicewp, rel, rei, emis, coszrs  )
      end if
   else

! convert radiative heating rates from Q*dp to Q for energy conservation
      if (conserve_energy) then
         do k =1 , pver
            do i = 1, ncol
               qrs(i,k) = qrs(i,k)/state%pdel(i,k)
               qrl(i,k) = qrl(i,k)/state%pdel(i,k)
            end do
         end do
      end if
!	  open(unit=99,file='qrl.txt', status='unknown') !EJL
!	  write(99,*) qrl
!	  close(unit=99)
         
   end if

! FAO:  temporarily disable radiation update
#ifndef RAD_NOUPDATE

!
! Compute net flux
! Since fsns, fsnt, flns, and flnt are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
#ifdef DIUR_AV
   do i=1,ncol
      tend%flx_net(i) = fsnt_da(i,lchnk) - fsns_da(i,lchnk) - flnt(i) + flns(i)
   end do
   ftem(:ncol,:pver)= qrs_da(:ncol,lchnk,:pver)/cpair 
   call outfld('QRS_DA  ',ftem    ,pcols,lchnk)
   tmp(:ncol)=fsns_da(:ncol,lchnk)
   call outfld('FSNS_DA ',tmp,pcols,lchnk)
   tmp(:ncol)=fsds_da(:ncol,lchnk)
   call outfld('FSDS_DA ',tmp,pcols,lchnk)
   tmp(:ncol)=fsnt_da(:ncol,lchnk)
   call outfld('FSNT_DA ',tmp,pcols,lchnk)
   tmp(:ncol)=sols_da(:ncol,lchnk)
   call outfld('SOLS_DA ',tmp,pcols,lchnk)
   tmp(:ncol)=solsd_da(:ncol,lchnk)
   call outfld('SOLSD_DA',tmp,pcols,lchnk)
#else
   do i=1,ncol
      tend%flx_net(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
   end do
#endif   
!
! Compute net radiative heating
!
#ifdef DIUR_AV
   call radheat_net (state, ptend, qrl, qrs_da(:,lchnk,:))
#else
   call radheat_net (state, ptend, qrl, qrs)
#endif
!
! Add radiation tendencies to cummulative model tendencies and update profiles
!
   call physics_update(state, tend, ptend, ztodt)

! check energy integrals
   call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, tend%flx_net)
!
! Compute net surface radiative flux for use by surface temperature code.
! Note that units have already been converted to mks in RADCTL.  Since
! fsns and flwds are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
#ifdef DIUR_AV
   srfrad(:ncol) = fsns_da(:ncol,lchnk) + flwds(:ncol)
#else
   srfrad(:ncol) = fsns(:ncol) + flwds(:ncol)
#endif
   call outfld('SRFRAD  ',srfrad,pcols,lchnk)

!!! Turn routine back on to allow surface forcing and output ajf,3/30/06:
#endif
!
! Save atmospheric fields to force surface models
!

   call srfxfer (lchnk, ncol, state%ps, state%u(1,pver), state%v(1,pver),    &
                 state%t(1,pver), state%q(1,pver,1), state%exner(1,pver), state%zm(1,pver), &
                    state%pmid,      &
                 state%rpdel(1,pver))

!---------------------------------------------------------------------------------------
! Save history variables. These should move to the appropriate parameterization interface
!---------------------------------------------------------------------------------------

   call outfld('PRECL   ',precl   ,pcols   ,lchnk       )
   call outfld('PRECC   ',precc   ,pcols   ,lchnk       )
   call outfld('PRECSL  ',precsl  ,pcols   ,lchnk       )
   call outfld('PRECSC  ',precsc  ,pcols   ,lchnk       )
   
   prect(:ncol) = precc(:ncol) + precl(:ncol)
   call outfld('PRECT   ',prect   ,pcols   ,lchnk       )
   call outfld('PRECTMX ',prect   ,pcols   ,lchnk       )

#if ( defined COUP_CSM )
   call outfld('PRECLav ',precl   ,pcols   ,lchnk   )
   call outfld('PRECCav ',precc   ,pcols   ,lchnk   )
#endif

#if ( defined BFB_CAM_SCAM_IOP )
   call outfld('Prec   ',prect   ,pcols   ,lchnk       )
#endif
!     
! Compute heating rate for dtheta/dt 
!
   do k=1,pver
      do i=1,ncol
         ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5/state%pmid(i,k))**cappa
      end do
   end do
   call outfld('HR      ',ftem    ,pcols   ,lchnk   )

! convert radiative heating rates to Q*dp for energy conservation
   if (conserve_energy) then
      do k =1 , pver
         do i = 1, ncol
            qrs(i,k) = qrs(i,k)*state%pdel(i,k)
            qrl(i,k) = qrl(i,k)*state%pdel(i,k)
         end do
      end do
   end if

   return
end subroutine tphysbc
