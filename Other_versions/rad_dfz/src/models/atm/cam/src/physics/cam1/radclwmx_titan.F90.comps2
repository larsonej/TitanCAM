#include <misc.h>
#include <params.h>

subroutine di_radflx(knu_ir,nzm,nzi,czar,emis,taui,wt, &
                      planckf,planckfg,fluxir,flns,flnt,flut,flntc,flutc,&
                      flwds,fcnl,fnl)
!=======================================================================
!
!  Compute the radiative flux at interfaces of the model (Direct Integration method)
!  Planck function calculated as a linear interpolarion from interface Planck functions
!
!=======================================================================
!  USES:

   use shr_kind_mod, only: i4 => SHR_KIND_I4, r8 => shr_kind_r8   
   
   implicit none
   
! Arguments
   
   integer(i4), intent(in) :: nzm,nzi, knu_ir
   real(r8), intent(in) :: taui(knu_ir,nzi)
   real(r8), intent(in) :: czar, emis(knu_ir)
   real(r8), intent(out) :: fluxir(knu_ir,nzi)  !IR Flux at interfaces
   real(r8), intent(in) :: planckf(knu_ir,nzi),planckfg(knu_ir)
   real(r8), parameter :: pi=3.141592654
   real(r8), intent(in) :: wt(knu_ir)  ! weights for spectral sums
   real(r8), intent(out) :: flns,flnt,flwds,fcnl(nzi),fnl(nzi),flut,flntc,flutc ! CAM variables
   real(r8)  temp

! qrl(loc,vl) --> longwave (lw) heating rate NOT COMPUTED HERE
! flns(loc)   --> surface cooling flux
! flnt(loc)   --> net outgoing flux at top
! flut(loc)   --> upward flux at model top
! flnsc(loc)  --> "clear-sky" surface cooling NOT COMPUTED HERE
! flntc(loc)  --> net "clear-sky" outgoing flux
! flutc(loc)  --> upward "clear-sky" flux at model top
! flwds(loc)  --> downward lw flux at surface
! fcnl(loc,vl)--> "clear-sky" net flux at interfaces
! fnl(loc,vl) --> net flux at interfaces

! Local variables:

   integer i,j,k,n,nz1,nz2, nz1i, nz2i
   real(r8) arad(knu_ir,nzi),brad(knu_ir,nzi),crad(knu_ir,nzi), rrad(nzi)
   real(r8) dtm,  vnu, cza2, tav
   real(r8) v   ! ~flux in two-stream approximation
   real(r8) aar(nzi),bbr(nzi),ccr(nzi) ! a,b,c used in tridiagonal matrix solver
   real(r8) reverse_planck,wave_number,angstrom,brightness_temp(knu_ir)
   logical cgs
   

! ---------------------------------------------------------------------

   real(r8)  b0(nzm),b1(nzm),mu,t1(nzm),expt1(nzm),df(nzi),uf(nzi),planckdiff(nzm)
   
   
    flns = 0.0
    flnt = 0.0
    flwds = 0.0
    fcnl(1:nzi) = 0.0
    fnl(1:nzi) = 0.0
   
   
   mu = czar
    do k=1,knu_ir
   
   t1 = taui(k,2:nzi)-taui(k,1:nzm)
   expt1 = exp(-t1/mu)
   planckdiff = (planckf(k,2:nzi)-planckf(k,1:nzm))
   where (t1 > 0.0)
   b1 = planckdiff/t1
   elsewhere
   b1 = 0.0
   endwhere
   b0 = planckf(k,1:nzm)
  
! df = downward flux, uf = upward flux
   df(1) = 0.0
   
   do j = 2,nzi
   i = j-1
   !print*,b0(i),b1(i)
   df(j) = df(i)*expt1(i) + 2*pi*mu*(b0(i)*(1-expt1(i)) + &
                            b1(i)*((expt1(i)-1)*mu + t1(i)))
   enddo
   uf(nzi) = 2.0*pi*mu*planckfg(k)*emis(k) + df(nzi)*(1.0-emis(k))
   !print*,'emis, df(nzi), un(nzi) = ',emis,df(nzi),uf(nzi)
   do j = 2,nzi
   
   n = nzi-j+1
   
    uf(n) = uf(n+1)*expt1(n) + &
        2*pi*mu*(b0(n)*(1-expt1(n)) + b1(n)*(mu -(mu+t1(n))*expt1(n)))         
            
   enddo
   fluxir(k,1:nzi) = uf-df

    fnl = fnl + fluxir(k,1:nzi)*wt(k)
    flwds = flwds + df(nzi)*wt(k)
    
 ! This section for diagnostic - comment out once code is confirmed
 ! calculate brightness temperature - check against IDL result
 ! 
 
 !cgs = .false.
 !wave_number = k*10.0
 !angstrom = 1.0E+08/wave_number
 !temp = fluxir(k,1)/(1.e+08/(wave_number-5.0) - 1.e+08/(wave_number + 5.0))
 !brightness_temp(k) = reverse_planck(angstrom,temp,cgs)
 

! end diagnostic section

 40 format(i5,1p2e10.2)  
 
    enddo

! print 50,brightness_temp ! more diagnostic
50 format(10f8.3)

    flns = fnl(nzi)
    flnt = fnl(1)
    flut = flnt
    flntc = flnt
    flutc = flnt
    fcnl = fnl


return
end subroutine di_radflx


function reverse_planck (wave,flux,cgs) Result (temp)

use shr_kind_mod, only: i4 => SHR_KIND_I4, r8 => shr_kind_r8   
   

logical, intent(IN), optional :: cgs 
! real, dimension(:), intent(IN) :: wave,flux
real(r8), intent(IN) :: wave,flux
!real, dimension(size(wave)):: w,f,temp
real(r8) w,f,temp
real(r8) :: c1,c2
logical :: units

!  Purpose:  to calculate the brightness temperature from the
!            wavelength (Angstroms) and flux watt/m2/A or Joule/m2/s/A
!  If argument cgs is present and .true. then flux is (ergs/cm2/s/A)

if (present(cgs)) then 
units = cgs 
else 
units = .false.
end if

w = wave / 1.E8                              ! Angstroms to cm
  c1 =  3.74185E-5             !constants appropriate to cgs units.
  C2 =  1.43883
  if(units) then 
  f = flux*1.e+08 ! input flux in erg/cm2/s/A
  else 
  !print*,'mks'
  f = flux*1.0e+11 ! input flux in Watt/m2/A = Joule/m2/s/A
  end if
!print*,f,wave,w,w**5,log(1.0 + c1/(f*w**5))
temp = (c2/w/log(1.0 + c1/(f*w**5)))

end function


subroutine radclwmx_titan(lchnk,ncol, mmr_gas, pint, tnm, planckf_in, & 
     planckfg_in, qrl,flns,flnt,flut,flnsc,flntc,flutc,flwds,fcnl,fnl)

   use shr_kind_mod, only: i4 => shr_kind_i4, r8 => shr_kind_r8
   !use shr_const_mod, only : shr_const_g
   use shr_const_mod
   use radcnst
   use ir_rad_mod
   use ppgrid
   use constituents, only: ppcnst
   use pmgrid, only: masterproc
  ! USE nrtype
! AJF 3-27-08: MODIFY FOR "INLINE", SPATIALLY VARIABLE LONGWAVE OD
   use taulw_cam, only: get_taulw_cam
  
   implicit none 

   integer,  intent(in) :: lchnk,ncol           ! number of chunk,atmospheric columns
   real(r8), intent(in) :: mmr_gas(pcols,pver,ppcnst) ! gas mass mixing ratios in chunk
   real(r8), intent(in) :: pint(pcols,pverp)    ! Model interface pressure
   real(r8), intent(in) :: planckf_in(Nf_lw,pverp,ncol) !atm planck function on interfcs
   real(r8), intent(in) :: planckfg_in(Nf_lw,ncol)     !surface planck function
   real(r8), intent(in) :: tnm(pcols,pver)       ! Layer temperatures

! FAO these are the inputs for the original radclwmx_titan:
!   integer, intent(in) :: lchnk                 ! chunk identifier
!   integer, intent(in) :: ncol                  ! number of atmospheric columns
!   integer, intent(in) :: nmxrgn(pcols)         ! Number of maximally overlapped regions
!
!   real(r8), intent(in) :: pmxrgn(pcols,pverp)  ! Maximum values of pmid for each
!   real(r8), intent(in) :: lwupcgs(pcols)       ! Longwave up flux in CGS units
!! Input arguments which are only passed to other routines
!   real(r8), intent(in) :: tnm(pcols,pver)      ! Level temperature
!   real(r8), intent(in) :: qnm(pcols,pver)      ! Level moisture field
!   real(r8), intent(in) :: o3vmr(pcols,pver)    ! ozone volume mixing ratio
!   real(r8), intent(in) :: pmid(pcols,pver)     ! Level pressure
!   real(r8), intent(in) :: pint(pcols,pverp)    ! Model interface pressure
!   real(r8), intent(in) :: pmln(pcols,pver)     ! Ln(pmid)
!   real(r8), intent(in) :: piln(pcols,pverp)    ! Ln(pint)
!   real(r8), intent(in) :: n2o(pcols,pver)      ! nitrous oxide mass mixing ratio
!   real(r8), intent(in) :: ch4(pcols,pver)      ! methane mass mixing ratio
!   real(r8), intent(in) :: cfc11(pcols,pver)    ! cfc11 mass mixing ratio
!   real(r8), intent(in) :: cfc12(pcols,pver)    ! cfc12 mass mixing ratio
!   real(r8), intent(in) :: cld(pcols,pver)      ! Cloud cover
!   real(r8), intent(in) :: emis(pcols,pver)     ! Cloud emissivity
!   real(r8), intent(in) :: aer_mass(pcols,pver) ! STRAER mass in layer


   real(r8), intent(out) :: qrl(pcols,pver)      ! Longwave heating rate
   real(r8), intent(out) :: flns(pcols)          ! Surface cooling flux
   real(r8), intent(out) :: flnt(pcols)          ! Net outgoing flux (at model top)
   real(r8), intent(out) :: flut(pcols)          ! Upward flux at top of model
   real(r8), intent(out) :: flnsc(pcols)         ! Clear sky surface cooling <NOT COMPUTED>
   real(r8), intent(out) :: flntc(pcols)         ! Net clear sky outgoing flux
   real(r8), intent(out) :: flutc(pcols)         ! Upward clear-sky flux at top of model
   real(r8), intent(out) :: flwds(pcols)         ! Down longwave flux at surface
   real(r8), intent(out) :: fcnl(pcols,pverp)    ! clear sky net flux at interfaces
   real(r8), intent(out) :: fnl(pcols,pverp)     ! net flux at interfaces

!  Local Variables:
   logical :: debug_rad
   real(r8), parameter :: czar=0.57
   real(r8), parameter :: pi8=3.141592654
   !real(r8), parameter :: emis_ir=1.0 ! ajf 7-2-07, fix dimensions of emis later
   real(r8), parameter :: emis_ir=TITAN_EMISSIVITY

   integer(i4) :: knu,nzm,nzi,loc
   real(r8) :: fnl_i(pverp), fcnl_i(pverp)
  
   integer(i4) :: i,j,k,b

   real(r8) :: planckf(Nf_lw,Ng_lw,pverp,ncol) !atm planck function on interfcs
   real(r8) :: planckfg(Nf_lw,Ng_lw,ncol)     !surface planck function
   real(r8) :: dnu(Nf_lw)

   integer(i4) :: Zf, Zg, i_nu, i_g
   
 
   ! set frequency bin size
   dnu(1) = wn_lw(2) - wn_lw(1)
   do Zf=2,Nf_lw-1
      dnu(Zf) = 0.5*(wn_lw(Zf+1) - wn_lw(Zf-1))
   enddo
   dnu(Nf_lw) = wn_lw(Nf_lw) - wn_lw(Nf_lw-1)

   do i=1,ncol
      do k=1,pverp
      do Zf=1,Nf_lw
         do Zg=1,Ng_lw
            i_nu = f_of_glob(Zg,Zf)
            i_g = g_of_glob(Zg,Zf)
            planckf(i_nu,i_g,k,i) = planckf_in(Zf,k,i) * weight_lw(i_g,i_nu) * dnu(i_nu)
         enddo
      enddo
      enddo
      do Zf=1,Nf_lw
         do Zg=1,Ng_lw
            i_nu = f_of_glob(Zg,Zf)
            i_g = g_of_glob(Zg,Zf)
            planckfg(i_nu,i_g,i) = planckfg_in(Zf,i)  * weight_lw(i_g,i_nu) * dnu(i_nu)
         enddo
      enddo
   enddo

   do i=1,ncol

     debug_rad=.false.
!!     if (( lchnk == 27 .and. i == 4) .or. ( lchnk == 28 .and. i == 4) ) then
!!      debug_rad=.true.
!!      print*,' DEBUG_RAD ON' 
!!      print*,' lchnk=, i= ',lchnk,i
!!!      k=1
!!!      do i_nu=1,10
!!!       print*,' i_nu, nu= ',i_nu,wn_lw(i_nu)
!!!       print*,' planckf_in= ',planckf_in(i_nu,k,i)
!!!      do i_g=1,num_g(i_nu)
!!!       print*,' i_g, weight_lw(i_g,i_nu), dnu(i_nu)= ',i_g, weight_lw(i_g,i_nu), dnu(i_nu)
!!!       print*,' i_g, planckfg(i_nu,i_g,i)= ',i_g, planckfg(i_nu,i_g,i)
!!!       print*,' planckf(i_nu,i_g,k,i)= ',planckf(i_nu,i_g,k,i)
!!!      enddo
!!!      enddo
!!     endif   
!!!     call t_startf("get_taulw_cam")
!!!!!       call get_taulw_cam(mmr_gas(i,:,:), pint(i,:), tnm(i,:), taui_lw, debug_rad)
!!!     call t_stopf("get_taulw_cam")

      call get_ir_totflxs(Nf_lw,Ng_lw,num_g,weight_lw,taui_lw(lchnk,i,:,:,:), &
           czar,emis_ir,planckf(:,:,:,i),planckfg(:,:,i), &
           fnl_i,flwds(i),flns(i),flnsc(i),flnt(i),flut(i),flntc(i),flutc(i),fcnl_i,debug_rad)
      
      do k=1,pverp
         fcnl(i,k) = fcnl_i(k)
         fnl(i,k) = fnl_i(k)
      enddo

!
! Compute thermal IR heating rate (W/kg)
!

      do k=1,pver
         qrl(i,k) = shr_const_g*(fnl(i,k)-fnl(i,k+1))/(pint(i,k) - pint(i,k+1))
      end do


!!!      if (debug_rad) then
!!!         print*,' IR heating rates: '
!!!         do k=1,pver
!!!          print*,' k,qrl(W/kg)= ',k,qrl(i,k)
!!!         enddo
!!!      endif

   enddo  ! end loop over i
   

   return

end subroutine radclwmx_titan


