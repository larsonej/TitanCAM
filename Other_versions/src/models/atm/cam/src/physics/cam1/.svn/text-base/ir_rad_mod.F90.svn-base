module ir_rad_mod

!=======================================================================
!
!  Contains routines to compute mean IR net flux and associated quantities
!  averaged over arbitrary spectral interval and a routine to sum band contributions
!  together to get total spectrally-averaged flux or irradiance quantity
! 
!
!=======================================================================

  use shr_kind_mod, only: i4 => SHR_KIND_I4, r8 => shr_kind_r8
  
  implicit none
  
  private
  save
!
! Public interfaces
!
  public get_ir_bndmean  ! computes mean irradiances over spectral bandpass
  public get_ir_totflxs  ! sum contributions of bandpasses to get total fluxes
  public planck

! Public data
 
  
CONTAINS

real(r8) function planck(vnu,temp)

!=======================================================================
!
! Computes planck function  at wavenumber vnu
! B = 2 h c^2 vnu^3 /( exp [hc vnu/kT] -1 )
! RETURN IN UNITS OF W m^-2/sr/cm^-1
!
!=======================================================================

!USES

  use shr_kind_mod, only: r8 => shr_kind_r8
  
   implicit none
!
! Input arguments
!
   real(r8), intent(in) :: vnu, temp  ! wavenumber, cm^-1 and T, K
!
! Local variables
!   
   real(r8), parameter :: cgs_to_mks=1.e-3
   real(r8), parameter :: c=3.e10 !speed of light
   real(r8), parameter :: hpl=6.626e-27 !planck's constant   
   real(r8), parameter :: c1=2.*hpl*c*c !coefficient in planck when using wavenumbers
   real(r8), parameter :: hck=hpl*c/1.38054e-16
   real(r8) xnu3,vot,xprtb

      xnu3=vnu*vnu*vnu
      vot=hck*vnu/temp
      xprtb=exp(vot)-1.0
      planck=cgs_to_mks*c1*xnu3/xprtb

return
end function planck


subroutine get_ir_bndmean(ns,nzm,nzi,mu,emis_ir,taui,gws, &
     planckf,planckfg,flwds,fnl, debug_rad)

!=======================================================================
!
!  Compute the radiative flux at interfaces of the model (Direct Integration method)
!  Source function linear-in-tau interpolation of interface Planck functions
!
!=======================================================================
!  USES:

   use shr_kind_mod, only: i4 => SHR_KIND_I4, r8 => shr_kind_r8 
   use ppgrid, only: pver, pverp

   implicit none
!!!!!!!!!!   integer(i4), parameter :: pver=26
!!!!!!!!!!   integer(i4), parameter :: pverp=27     
   

   
! Arguments
   
   logical, intent(in) :: debug_rad
   integer(i4), intent(in) :: nzm, nzi  !number of layers, interfaces
   integer(i4), intent(in) :: ns        !number of spectral points in band
   real(r8), intent(in) :: taui(ns,nzi) !opt. depths at interfaces 
   real(r8), intent(in) :: mu           ! cos-zenith-angle
   real(r8), intent(in) :: emis_ir     ! allow to vary with wn later 
!   real(r8), intent(in) :: emis(ns)     ! local surface emissivity
   real(r8), intent(in)  :: planckf(ns,nzi), planckfg(ns)
   real(r8), intent(in) ::  gws(ns)     ! band weights for spectral sums
! The following irradiances store averages over spectral interval [nu1,nu2]
   real(r8), intent(out) :: flwds    ! downward lw flux at surface
   real(r8), intent(out) :: fnl(nzi) ! net flux at interfaces


   real(r8), parameter :: pi=3.141592654

! Local variables:

   integer i,j,k,n,nz1,nz2, nz1i, nz2i

   real(r8)  b0(pver),b1(pver),t1(pver),expt1(pver),df(pverp),uf(pverp)
   real(r8)  planckdiff(pver)
   real(r8)  fluxir(pverp)  !IR Flux at interfaces      
   

! ---------------------------------------------------------------------
   

    flwds = 0.0
    fnl(1:nzi) = 0.0
   
   
  do k=1,ns
   
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

!!   uf(nzi) = 2.0*pi*mu*planckfg(k)*emis(k) + df(nzi)*(1.0-emis(k))
   ! FAO -- was:
!!   uf(nzi) = 2.0*pi*mu*planckfg(k)*emis + df(nzi)*(1.0-emis)
!!!! In general, we want emis to depend on wavenumber.  But we
!!!! currently have a dimensioning problem, so a scalar emis_ir=1 is
!!!! being passed thru from radclwmx_titan down to here:

   uf(nzi) = 2.0*pi*mu*planckfg(k)*emis_ir + df(nzi)*(1.0-emis_ir)

   !print*,'emis, df(nzi), un(nzi) = ',emis,df(nzi),uf(nzi)
   do j = 2,nzi
   
   n = nzi-j+1
   !uf(n) = uf(n+1)*exp(-t1(n)/mu) + &
   !     2*pi*mu*(b0(n)*(1-exp(-t1(n)/mu)) + b1(n)*(mu -(mu+t1(n))*exp(-t1(n)/mu)))
    uf(n) = uf(n+1)*expt1(n) + &
        2*pi*mu*(b0(n)*(1-expt1(n)) + b1(n)*(mu -(mu+t1(n))*expt1(n)))         

   enddo
   fluxir = uf-df

!!!  if (debug_rad) then
!!!   print*,' k, emis_ir, planckfg(k)= ',k,emis_ir,planckfg(k)
!!!   do j = 1,nzm
!!!    print*,' b0, b1, t1, expt1= ',b0(j),b1(j),t1(j),expt1(j)
!!!   enddo
!!!   do j = 1,nzi
!!!    print*,'j, df,uf, fluxir= ',j,df(j),uf(j),fluxir(j)
!!!   enddo
!!!   j=nzm
!!!   print*,'dfnzi1,dfnzi2,dfnzi3= ',df(j)*expt1(j), &
!!!      b0(j)*(1-expt1(j)),b1(j)*((expt1(j)-1)*mu + t1(j))
!!!  endif

    fnl = fnl + fluxir(1:nzi)
    flwds = flwds + df(nzi)
!!!   if (debug_rad) then
!!!    print*, 'k, flwds= ',k,flwds
!!!    print*,' k, fnl(1)= ',k,fnl(1)
!!!   endif

!!! fnl = fnl + fluxir(1:nzi)*gws(k) !don't need gws if planck already weighted
!!! flwds = flwds + df(nzi)*gws(k)
    
 ! This section for diagnostic - comment out once code is confirmed
 ! calculate brightness temperature - check against IDL result
 ! 
 
 !cgs = .false.
 !wave_number = k*10.0
 !angstrom = 1.0E+08/wave_number
 !temp = fluxir(k,1)/(1.e+08/(wave_number-5.0) - 1.e+08/(wave_number + 5.0))
 !brightness_temp(k) = reverse_planck(angstrom,temp,cgs)
 

! end diagnostic section
    

 enddo    !end loop over k (spectral points)


return
end subroutine get_ir_bndmean


!==================================================================================
!==================================================================================



subroutine get_ir_totflxs(num_ir_band,nsm,nsa,wts,taui,mu,emis_ir,  &
     planckf,planckfg,fnl,flwds,flns,flnsc,flnt,flut,flntc,flutc,fcnl,debug_rad)
     
!=======================================================================
!
!  Compute total spectral-mean irradiances required by CAM model
!   - calls get_ir_bndmean to get spectral mean in single bandpass 
!
!=======================================================================        
        
! USES

   use shr_kind_mod, only: i4 => SHR_KIND_I4, r8 => shr_kind_r8
   use ppgrid, only: pver, pverp

   implicit none
!!!!!!!!!!!   integer(i4), parameter :: pver=26
!!!!!!!!!!!   integer(i4), parameter :: pverp=27

!
! Input arguments
!
   logical, intent(in) ::  debug_rad
   integer(i4), intent(in) :: num_ir_band  !number of distinct spectral intervals in IR
   integer(i4), intent(in) :: nsm         !max dimension of spectral points in wts, etc.
   integer(i4), intent(in) :: nsa(num_ir_band) ! number of spectral points in each band
   real(r8), intent(in) :: wts(nsm,num_ir_band)  !gaussian wts of each spect interval
   real(r8), intent(in) :: taui(num_ir_band,nsm,pverp)  !optical depth array
   real(r8), intent(in) :: planckf(num_ir_band,nsm,pverp) !atm planck function on interfcs
   real(r8), intent(in) :: planckfg(num_ir_band,nsm)     !surface planck function
   real(r8), intent(in) :: mu                          !2-stream za
   real(r8), intent(in) :: emis_ir
!   real(r8), intent(in) :: emis(nsm)                   !surface emissivity    

!
! Output arguments
!
! The following irradiances store *total* spectral integration 
   real(r8), intent(out) :: flns     ! surface cooling flux
   real(r8), intent(out) :: flnsc    ! "clear-sky" surface cooling flux
   real(r8), intent(out) :: flnt     ! net outgoing flux at top
   real(r8), intent(out) :: flwds    ! downward lw flux at surface
   real(r8), intent(out) :: fcnl(pverp)! "clear-sky" net flux at interfaces
   real(r8), intent(out) :: fnl(pverp) ! net flux at interfaces
   real(r8), intent(out) :: flut     ! upward flux at model top
   real(r8), intent(out) :: flntc    ! net "clear-sky" outgoing flux
   real(r8), intent(out) :: flutc    ! upward "clear-sky" flux at model top
   
!  Local Variables

   integer(i4) ::  ibnd
    integer(i4) :: ns,k  !for debug purposes, ajf 6/29/07
   real(r8) :: fnl_bnd(pverp),flwds_bnd
   real(r8) :: fnl_sum(pverp),flwds_sum
! --------------------------------------------------------------------------------------

    
    fnl_sum(:)=0.0
    flwds_sum =0.0

!  Loop over each spectral interval

    do ibnd=1,num_ir_band

!!!    if (debug_rad) print*,'get_totflxs: ibnd= ',ibnd,' nsa(ibnd)= ',nsa(ibnd)

     call get_ir_bndmean(nsa(ibnd),pver,pverp,mu,emis_ir,  &
                         taui(ibnd,1:nsa(ibnd),:),wts(1:nsa(ibnd),ibnd), &
                         planckf(ibnd,1:nsa(ibnd),:),planckfg(ibnd,1:nsa(ibnd)), &
                         flwds_bnd,fnl_bnd,debug_rad)
                        
!!!     if (debug_rad) then 
!!!         print*,' ibnd,fnl_bnd(1)= ',ibnd,fnl_bnd(1)
!!!     endif

     fnl_sum(:) = fnl_sum(:) + fnl_bnd(:)
     flwds_sum  = flwds_sum  + flwds_bnd
     
    enddo !end of ibnd loop
    
    fnl = fnl_sum
    flwds = flwds_sum 
    flns = fnl(pverp)
!!!    if (debug_rad) then
!!!      print*,' get_totflxs: flwds= ',flwds              
!!!      print*,'get_totflxs: flns= ',flns
!!!      print*,'get_tot_flxs: nu-integred net flux:'
!!!      do k=1,pverp
!!!       print*,' k, fnl(k)= ',k,fnl(k)
!!!      enddo
!!!    endif

    flnsc= flns
    flnt = fnl(1)
    flut = flnt
    flntc = flnt
    flutc = flnt
    fcnl = fnl

    
return
end subroutine get_ir_totflxs
       

end module ir_rad_mod        
        
