#include <misc.h>
#include <params.h>

subroutine tides_titan_eqtorq (ztodt   ,taux    ,tauy    , state   , &
                     tend, ptend, flx_bndry)  !ajf 6/26/07 :  tend changed to ptend for output
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  algorithm 1: Held/Suarez IDEALIZED physics
!  algorithm 2: Held/Suarez IDEALIZED physics (Williamson modified stratosphere
!  algorithm 3: Held/Suarez IDEALIZED physics (Lin/Williamson modified strato/meso-sphere
!  algorithm 4: Boer/Denis  IDEALIZED physics
!  algorithm 5: For Titan 9-1-05 (ajf)
!
! Author: J. Olson
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid            , only: plev,plat,plevp,masterproc
   use ppgrid
   use phys_grid         , only: get_lat_all_p, get_rlat_all_p, get_rlon_all_p
   use vertical_diffusion, only: vd_intr
   use physics_types     , only: physics_state, physics_tend, physics_ptend
   use geopotential      , only: geopotential_t
   use history,            only: outfld
   use physconst,          only: gravit, cappa, rair, cpair
   use constituents,       only: pcnst, pnats
   use abortutils,         only: endrun
   use time_manager,       only: get_curr_date, get_start_date
   use dycore,             only: dycore_is

   implicit none

#include <comhyb.h>
!
! Input arguments
!
   real(r8), intent(in) :: ztodt   ! Two times model timestep (2 delta-t)
!                                  ! except when dyn core='LR', where =delta-t
!
! Output arguments
!
   real(r8), intent(out) :: taux(pcols)            ! X surface stress (zonal)
   real(r8), intent(out) :: tauy(pcols)            ! Y surface stress (meridional)
   real(r8), intent(out) :: flx_bndry(pcols)
   
   type(physics_state), intent(inout) :: state
   type(physics_tend), intent(inout)  :: tend
   type(physics_ptend), intent(out)   :: ptend !indivdual param tendencies

!
!---------------------------Local workspace-----------------------------
!

   integer :: lchnk                                ! chunk identifier
   integer :: ncol                                 ! number of atmospheric columns

   real(r8) clat(pcols)                        ! latitudes(radians) for columns
   real(r8) clon(pcols)                        ! longitudes(radians) for columns
   real(r8) pmid(pcols,pver)                   ! mid-point pressure
   integer  i,k                                ! Longitude, level indices
   integer  nstp                               ! nstep  ,ajf 11-2-05

   real(r8) tmp                                ! temporary
   real(r8) kf                                 ! 1./efolding_time for wind dissipation
   real(r8) ka                                 ! 1./efolding_time for temperature diss.
   real(r8) kaa                                ! 1./efolding_time for temperature diss.
   real(r8) ks                                 ! 1./efolding_time for temperature diss.
   real(r8) kv                                 ! 1./efolding_time (normalized) for wind
   real(r8) kt                                 ! 1./efolding_time for temperature diss.
   real(r8) trefa                              ! "radiative equilibrium" T
   real(r8) trefc                              ! used in calc of "radiative equilibrium" T
   real(r8) cossq(pcols)                       ! coslat**2
   real(r8) cossqsq(pcols)                     ! coslat**4
   real(r8) sinsq(pcols)                       ! sinlat**2
   real(r8) onemsig                            ! 1. - sigma_reference
   real(r8) efoldf                             ! efolding time for wind dissipation
   real(r8) efolda                             ! efolding time for T dissipation
   real(r8) efoldaa                            ! efolding time for T dissipation
   real(r8) efolds                             ! efolding time for T dissipation
   real(r8) efold_strat                        ! efolding time for T dissipation in Strat
   real(r8) efold_meso                         ! efolding time for T dissipation in Meso
   real(r8) efoldv                             ! efolding time for wind dissipation
   real(r8) p_infint                           ! effective top of model
   real(r8) constw                             ! constant
   real(r8) lapsew                             ! lapse rate
   real(r8) p0strat                            ! threshold pressure
   real(r8) phi0                               ! threshold latitude
   real(r8) dphi0                              ! del-latitude
   real(r8) a0                                 ! coefficient
   real(r8) aeq                                ! 100 mb
   real(r8) apole                              ! 2   mb
   real(r8) pi                                 ! 3.14159...
   real(r8) coslat(pcols)                      ! cosine(latitude)
   real(r8) acoslat                            ! abs(acos(coslat))
   real(r8) constc                             ! constant
   real(r8) lapsec                             ! lapse rate
   real(r8) lapse                              ! lapse rate
   real(r8) h0                                 ! scale height (7 km)
   real(r8) sigmab                             ! threshold sigma level
   real(r8) pressmb                            ! model pressure in mb
   real(r8) t00                                ! minimum reference temperature
   real(r8) etamid(pver)                       ! midpoint values of eta (a+b)

! FOR TITAN: (ajf)
   real(r8) ds0i,ds1i,ds2i,ds3i,ds4i
   real(r8) ds0,ds1,ds2,ds3,ds4
   real(r8) thu0,thu1,thu2,thu3,thu4
   real(r8) yu0,yu1,yu2,yu3,yu4,yu,theta,y_cut,t_cut,tbar,kapa
   real(r8) dt0,y_0,dt_max,sig_max,sig_u,trefa1,trefa2,trefa3,trefa4
   real(r8) n_sat,L_s,c_L_s,s_L_s
! INCL GRAV TIDES: (ajf)
   real(r8) cs2lon(pcols),sn2lon(pcols),sn2lat(pcols),csqlon(pcols)
   real(r8) utid,vtid,A_tide,c_n_t,s_n_t,n_t,n_titan
   real(r8) tim_fac                            ! =0.5 (eul core), 1.0 LR core
   real(r8) time_curr                          ! current time
   real(r8) time_s                             ! Start time of initial run, s
   integer  yr, mon, day, ncsec               ! Date parms at beginning of time step
   integer  yr_s,mon_s,day_s,ncsec_s          ! Start date at beginning of initial run
   real(r8) rho_utid(pcols,pver),rho_vtid(pcols,pver) !Store dens*tidal acc 
   integer  idlflag                            ! Flag to choose which idealized physics
!  ADD A FICTITIOUS EQUATORIAL TORQUE TO TEST CYCLOSTROPHIC REGIME BEHAVIOR
   real(r8) u_eq(pcols,pver),t_eq,ut_eqvp,g_ex(pcols), u_eq_0



!
!-----------------------------------------------------------------------
!
   idlflag = 5

   lchnk = state%lchnk
   ncol  = state%ncol
!
! Copy pressures into local array
!
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)
   do i=1,ncol
      coslat (i) = cos(clat(i))
      sn2lat (i) = sin(2.*clat(i))
!      sinsq  (i) = sin(clat(i))*sin(clat(i))
      sinsq(i) = sin(clat(i))**2
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq(i)
      cs2lon (i) = cos(2.*clon(i))
      sn2lon (i) = sin(2.*clon(i))
!      csqlon (i) = cos(clon(i))*cos(clon(i))
      csqlon(i) = cos(clon(i))**2
   end do
   do k=1,pver
      do i=1,ncol
         pmid(i,k) = state%pmid(i,k)
      end do
   end do

!  Fictitious equatorial torque:

   pi=3.1415926
   t_eq=1.e9 !s, experiment 5f:  Test separate acceleration source for strat.
!!!   t_eq=1.8e12 !s, ~570 years, for small, steady acceleration experiment
!   t_eq=2.33e8 !s, 1/4 Titan year
   phi0=30./180.*pi 
!!   phi0=15./180.*pi
!!   u_eq_0=1.56  !Default is 1.56, extra factor to enhance acceleration
   u_eq_0 = 1.56*8 !EJL 2-18-13 increasing torq factor
!! EJL 3-23-13 Changed scale factor above to 8 from 4 and 
!!   am running under name lav-eqtorq2
   do i=1,ncol
    g_ex(i)=exp(-(clat(i)*clat(i))/(phi0*phi0))
    do k=1,pver
     ut_eqvp=u_eq_0*(100.0*(pmid(i,k)/500.0)**(-0.577) - 3.771)
     if (pmid(i,k) < 500.0) ut_eqvp=u_eq_0*(100.0 - 3.771)
     u_eq(i,k) = ut_eqvp
    enddo
   enddo

!!   if (masterproc) then
!!    i=1
!!    print*,'i, clat(i), g_ex= ',i,clat(i),g_ex(i)
!!    print*,' EQ_TORQ: k, pmid, u_eq '
!!    do k=1,pver
!!     print*,k,pmid(i,k),u_eq(i,k)
!!    enddo
!!   endif



   if (idlflag == 1) then
!
!-----------------------------------------------------------------------
!
! Held/Suarez IDEALIZED physics algorithm:
!
!   Held, I. M., and M. J. Suarez, 1994: A proposal for the
!   intercomparison of the dynamical cores of atmospheric general
!   circulation models.
!   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
!
!-----------------------------------------------------------------------
!
! Add idealized radiative heating rates to temperature tendency
!
      efoldf =  1.
      efolda = 40.
      efolds =  4.
      sigmab =  0.7
      t00    = 200.
!
      onemsig = 1. - sigmab
!
      ka = 1./(86400.*efolda)
      ks = 1./(86400.*efolds)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            do i=1,ncol
               kt = ka + (ks - ka)*cossqsq(i)*(etamid(k) - sigmab)/onemsig
               tmp   = kt/(1.+ ztodt*kt)
               trefc   = 315. - 60.*sinsq(i)
               trefa = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**cappa
               trefa    = max(t00,trefa)
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         else
            tmp   = ka/(1.+ ztodt*ka)
            do i=1,ncol
               trefc   = 315. - 60.*sinsq(i)
               trefa = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**cappa
               trefa    = max(t00,trefa)
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0.
            ptend%v(i,k) = 0.
         end do
      end do
      do i=1,pcols
         taux(i) = 0.
         tauy(i) = 0.
      end do
!
      kf = 1./(86400.*efoldf)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            kv  = kf*(etamid(k) - sigmab)/onemsig
            tmp = -kv/(1.+ ztodt*kv)
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
               tend%dudt(i,k)  = tend%dudt(i,k) + ptend%u(i,k)
               tend%dvdt(i,k)  = tend%dvdt(i,k) + ptend%v(i,k)
            end do
         endif
      end do

   elseif (idlflag == 2) then
!
!-----------------------------------------------------------------------
!
! Modified Held/Suarez IDEALIZED physics algorithm
! (modified with Williamson stratosphere):
!
!   Williamson, D. L., J. G. Olson and B. A. Boville, 1998: A comparison
!   of semi--Lagrangian and Eulerian tropical climate simulations.
!   Mon. Wea. Rev., vol 126, pp. 1001-1012.
!
!-----------------------------------------------------------------------
!
! Add idealized radiative heating rates to temperature tendency
!
      efoldf  =  1.
      efolda  = 40.
      efoldaa = 40.
      efolds  =  4.
      sigmab  =  0.7
      t00     = 200.
!
      onemsig = 1. - sigmab
!
      ka  = 1./(86400.*efolda)
      kaa = 1./(86400.*efoldaa)
      ks  = 1./(86400.*efolds)
!
      pi     = 4.*atan(1.)
      phi0   = 60.*pi/180.
      dphi0  = 15.*pi/180.
      a0     = 2.65/dphi0
      aeq    = 10000.
      apole  = 200.
      lapsew = -3.345e-03
      constw = rair*lapsew/gravit
      lapsec =  2.00e-03
      constc = rair*lapsec/gravit
      do k=1,pver
         if (etamid(k) > sigmab) then
            do i=1,ncol
               kt = ka + (ks - ka)*cossqsq(i)*(etamid(k) - sigmab)/onemsig
               acoslat = abs(acos(coslat(i)))
               p0strat = aeq - (aeq - apole)*0.5*(1. + tanh(a0*(acoslat - phi0)))
               tmp     = kt/(1.+ ztodt*kt)
               trefc   = 315. - 60.*sinsq(i)
               trefa   = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)** &
                                                                                     cappa
               trefa   = max(t00,trefa)
               if (pmid(i,k) < 10000.) then
                  trefa = t00*((pmid(i,k)/10000.))**constc
                  tmp   = kaa/(1.+ ztodt*kaa)
               endif
               if (pmid(i,k) < p0strat) then
                  trefa = trefa + t00*( ((pmid(i,k)/p0strat))**constw - 1. )
                  tmp   = kaa/(1.+ ztodt*kaa)
               endif
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         else
            do i=1,ncol
               acoslat = abs(acos(coslat(i)))
               p0strat = aeq - (aeq - apole)*0.5*(1. + tanh(a0*(acoslat - phi0)))
               tmp     = ka/(1.+ ztodt*ka)
               trefc   = 315. - 60.*sinsq(i)
               trefa   = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)** &
                                                                                     cappa
               trefa   = max(t00,trefa)
               if (pmid(i,k) < 10000.) then
                  trefa = t00*((pmid(i,k)/10000.))**constc
                  tmp   = kaa/(1.+ ztodt*kaa)
               endif
               if (pmid(i,k) < p0strat) then
                  trefa = trefa + t00*( ((pmid(i,k)/p0strat))**constw - 1. )
                  tmp   = kaa/(1.+ ztodt*kaa)
               endif
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0.
            ptend%v(i,k) = 0.
         end do
      end do
      do i=1,pcols
         taux(i) = 0.
         tauy(i) = 0.
      end do
!
      kf = 1./(86400.*efoldf)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            kv  = kf*(etamid(k) - sigmab)/onemsig
            tmp = -kv/(1.+ ztodt*kv)
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
               tend%dudt(i,k)  = tend%dudt(i,k) + ptend%u(i,k)
               tend%dvdt(i,k)  = tend%dvdt(i,k) + ptend%v(i,k)
            end do
         endif
      end do

   elseif (idlflag == 3) then
!
!-----------------------------------------------------------------------
!
! Held/Suarez IDEALIZED physics algorithm:
! (modified with Lin/Williamson stratosphere/mesosphere):
!
!   Held, I. M., and M. J. Suarez, 1994: A proposal for the
!   intercomparison of the dynamical cores of atmospheric general
!   circulation models.
!   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
!
!-----------------------------------------------------------------------
!
! Add idealized radiative heating rates to temperature tendency
!
      efoldf      =  1.
      efolda      = 40.
      efolds      =  4.
      efold_strat = 40.
      efold_meso  = 10.
      efoldv      = 0.5
      sigmab      = 0.7
      lapse       = 0.00225
      h0          = 7000.
      t00         = 200.
      p_infint    = 0.01
!
      onemsig = 1. - sigmab
!
      ka = 1./(86400.*efolda)
      ks = 1./(86400.*efolds)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            do i=1,ncol
               kt    = ka + (ks - ka)*cossqsq(i)*(etamid(k) - sigmab)/onemsig
               tmp   = kt/(1.+ ztodt*kt)
               trefc = 315. - 60.*sinsq(i)
               trefa = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**cappa
               trefa = max(t00,trefa)
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         else
            do i=1,ncol
               tmp     = ka/(1.+ ztodt*ka)
               pressmb = pmid(i,k)*0.01
               trefc   = 315. - 60.*sinsq(i)
               trefa   = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)** &
                                                                                     cappa
               trefa   = max(t00,trefa)
               if (pressmb <= 100. .and. pressmb > 1.) then
                  trefa = t00 + lapse*h0*coslat(i)*log(100./pressmb)
                  tmp   = (efold_strat-efold_meso)*log(pressmb)/log(100.)
                  tmp   = efold_meso + tmp
                  tmp   = 1./(86400.*tmp)
                  tmp   = tmp/(1.+ ztodt*tmp)
               endif
               if (pressmb <= 1. .and. pressmb > 0.01) then
                  trefa = t00 + lapse*h0*coslat(i)*log(100.*pressmb)
                  tmp   = 1./(86400.*efold_meso)
                  tmp   = tmp/(1.+ ztodt*tmp)
               endif
               if (pressmb <= 0.01) then
                  tmp   = 1./(86400.*efold_meso)
                  tmp   = tmp/(1.+ ztodt*tmp)
               endif
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0.
            ptend%v(i,k) = 0.
         end do
      end do
      do i=1,pcols
         taux(i) = 0.
         tauy(i) = 0.
      end do
!
      kf = 1./(86400.*efoldf)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            kv  = kf*(etamid(k) - sigmab)/onemsig
            tmp = -kv/(1.+ ztodt*kv)
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
               tend%dudt(i,k)  = tend%dudt(i,k) + ptend%u(i,k)
               tend%dvdt(i,k)  = tend%dvdt(i,k) + ptend%v(i,k)
            end do
         else
            do i=1,ncol
               pressmb  = pmid(i,k)*0.01
               if (pressmb <= 100.) then
                  kv       = 1./(86400.*efoldv)
                  tmp      = 1. + tanh(1.5*log10(p_infint/pressmb))
                  kv       = kv*tmp
                  tmp      = -kv/(1.+ ztodt*kv)
                  ptend%u(i,k) = tmp*state%u(i,k)
                  ptend%v(i,k) = tmp*state%v(i,k)
                  tend%dudt(i,k)  = tend%dudt(i,k) + ptend%u(i,k)
                  tend%dvdt(i,k)  = tend%dvdt(i,k) + ptend%v(i,k)
               endif
            end do
         endif
      end do
      elseif (idlflag == 5) then
!
!-----------------------------------------------------------------------
!
!  TITAN PARAMETERIZATION
!
!-----------------------------------------------------------------------
!
! 
! Gravitational Tides only 
!
!  Find  n_t :

      tim_fac=0.5
      if ( dycore_is('LR') ) tim_fac=1.0  ! ztodt=delta-t for LR dyn core !!      
      pi=3.14159
      n_titan=4.56e-6 !s^-1                     
!!!      nstp=get_nstep()
      call get_curr_date(yr, mon, day, ncsec)
      time_curr=3.1536e7*yr+2.592e6*mon+86400.0*day+ncsec !coeffs must be accurate to 3 decimal places
      call get_start_date(yr_s, mon_s, day_s, ncsec_s)
      time_s=3.1536e7*yr_s+2.592e6*mon_s+86400.0*day_s+ncsec_s 
      n_t=(time_curr-time_s)*n_titan
!!!      if (masterproc) then
!!!        print*,' yr= ',yr,' mon= ',mon,' day= ',day,' ncsec= ',ncsec, &
!!!               'time_curr= ',time_curr
!!!        print*,' yr_s= ',yr_s,' mon_s= ',mon_s,' day_s= ',day_s,' ncsec_s= ',ncsec_s, &
!!!               'time_s= ',time_s
!!!        print*,'n_t= ',n_t
!!!      endif
        
      s_n_t=sin(n_t)
      c_n_t=cos(n_t)

      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0.
            ptend%v(i,k) = 0.
         end do
      end do
!
!  Define tidal amplitude constant (in MKS units):
!             3 *     G        M_s  /        a^2      *  (R_T/a)        * ecc
      A_tide= 3.* 6.673e-11*5.685e26/(1.222e9*1.222e9)*(2.575e6/1.222e9)*0.029
      do k=1,pver
            do i=1,ncol
               utid=A_tide*coslat(i)*(-1.5*sn2lon(i)*c_n_t + 2.*cs2lon(i)*s_n_t )
               vtid=-A_tide*sn2lat(i)*(1.5*csqlon(i)*c_n_t +    sn2lon(i)*s_n_t )               
!  Add ad hoc eq torq here:
!!EJL 11-1-13 - run with only stratospheric forcing
              if (pmid(i,k) .lt. 500.0) then 
!              if (pmid(i,k) .lt. 500.0 .or. pmid(i,k) .ge. 1.e4) then !p<5mb or > 0.1 bar 
               utid=utid+g_ex(i)*u_eq(i,k)/t_eq  !small, steady acceleration
 
              endif
               
!!!               utid=utid+g_ex(i)*(u_eq(i,k)-state%u(i,k))/t_eq  !relaxation to obs u


               ptend%u(i,k) =  utid
               ptend%v(i,k) =  vtid
               rho_utid(i,k)=  state%pdel(i,k)/gravit*utid
               rho_vtid(i,k)=  state%pdel(i,k)/gravit*vtid
            end do
      end do

      do i=1,ncol  !Find KE boundary flux due to tides, ajf 9/17/09
       flx_bndry(i)=0.0
       do k=1,pver
       flx_bndry(i)=flx_bndry(i)+ &
         (state%u(i,k)*ptend%u(i,k)+state%v(i,k)*ptend%v(i,k)) &
          *state%pdel(i,k)/gravit
       enddo
      enddo

      ptend%name='tides_titan_eqtorq'
      ptend%lu=.true.
      ptend%lv=.true.

!!!      call outfld ('UTID', rho_utid, pcols, lchnk)
!!!      call outfld ('VTID', rho_vtid, pcols, lchnk)

   else
      write(6,*) 'TPHYSIDL: flag for choosing desired type of idealized ', &
                 'physics ("idlflag") is set incorrectly.'
      write(6,*) 'The valid options are 1, 2, or 3.'
      write(6,*) 'idlflag is currently set to: ',idlflag
      call endrun
   endif

   return
end subroutine tides_titan_eqtorq

