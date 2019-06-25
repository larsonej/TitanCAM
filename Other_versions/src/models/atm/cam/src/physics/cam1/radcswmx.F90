#include <misc.h>
#include <params.h>



subroutine rayleigh_n2(lambda,density,sigma, sz)

  implicit none

  ! This procedure calculates the Rayleigh scattering cross section for molecular
  ! nitrogen = sigma (square cm) at wavelength lambda (micron), n2 density (molecules/cc)
  !   From Sneep and Ubachs, Journal of Quantitative Spectroscopy &
  !       Radiative Transfer 92 (2005) 293-310.
  
  integer sz
  double precision, intent(in) :: lambda(sz)  ! in microns
  double precision, intent(in) :: density     ! in molecules / cm^3
  double precision, intent(out) :: sigma(sz)  ! in cm^2

  ! local variables
  double precision :: pi
  
  double precision, dimension(sz) :: nu, nusq, Fk, n, top, n_nusq, temp
  
  pi = 4.0 * atan(1.0)

  !sz = size(lambda)

  !allocate(nu(sz), nusq(sz), Fk(sz), n(sz), top(sz), n_nusq(sz), temp(sz))

  ! try some default values for testing

  !density = 2.547e+19 ; Reference p. 304
  !lambda = [10000./18788.4, 1.0] ! lambda can be an array of wavelengths or just a single number
  !lambda = 10000./18788.4 ! see reference table 2

  nu = 10000.0 / lambda ! Converts lambda (micron) to frequency (cm^{-1})
  nusq = nu**2
  Fk = 1.034 + 3.17e-12*nusq ! Eqn. 12

  n = lambda 
  top = n
  n_nusq = n

  where (nu < 21360.0)
     temp = 1.0d-08*(6498.2 + 307.43305D12/(14.4d09-nusq))*density/2.547e+19 ! Eqn. 10
     n = 1.0 + temp
     top = 2.0*temp + temp*temp ! = n^2 - 1
     n_nusq = n*n
  elsewhere
     temp = 1.e-08*(5677.465 + 318.81874D12/(14.4e09-nusq))*density/2.547e+19 ! Eqn. 11
     n = 1.0 + temp
     top = temp*(2.0 + temp) ! = n^2 - 1
     n_nusq = n*n
  endwhere

  sigma = (24*pi**3*(nusq/density)**2) * ((top/(n_nusq+2.0))**2) * Fk ! Eqn. 2
  
  !deallocate(nu, nusq, Fk, n, top, n_nusq, temp)
  
  !if n_params() lt 3 then
  ! print, 'sigma = ',sigma
  ! print, ' the test case sigma should be 5.30E-027'
  !endif

end subroutine rayleigh_n2


subroutine tau_rayleigh_n2(tau, lambda, mmr, pmid, pflx, tmid, Nday, Nlambda)

  use physconst, only: gravit ! m/s^2
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only : SHR_CONST_BOLTZ
  use ppgrid
  use pmgrid      , only : masterproc

  implicit none
  
  real(r8), parameter :: m_N2 = 4.6495088e-23  ! mass of N2 (g)

  integer,  intent(in) :: Nday, Nlambda
  real(r8), intent(in) :: pflx(pcols,0:pverp) ! Interface press, including extra layer (dyne/cm^2)
  real(r8), intent(in) :: lambda(Nlambda)  ! wavelength
  real(r8), intent(in) :: mmr(pcols,pver) ! gas mass mixing ratio of CH4
  real(r8), intent(in) :: pmid(pcols,pver) ! Level pressure (dyne/cm^2)
  real(r8), intent(in) :: tmid(pcols,pver) ! Level temperature (K)

  real(r8), intent(out) :: tau(pcols,0:pver,Nlambda)

  integer i,k,f
  real(r8) :: pmid_0
  real(r8) :: u, density
  real(r8) :: sigma(Nlambda) ! scattering cross section

  tau(:,0,:) = 0.0d0

  do k=1,pver
  do i=1,Nday
     ! FAO:  kludge -- assume mmr(N2) = 1 - mmr(CH4)
     u = (1.0 - mmr(i,k)) * (pflx(i,k+1) - pflx(i,k)) / (1.0e2*gravit)   ! g / cm^2	 
     density = (0.1 * pmid(i,k)) / (SHR_CONST_BOLTZ * tmid(i,k) ) * 1.0e-6  ! molecules / cm^3
     call rayleigh_n2(lambda, density, sigma, Nlambda)

     do f=1,Nlambda
        tau(i,k,f) = (sigma(f) / m_N2) * u;

        if (i .eq. 1 .and. masterproc .and. lambda(f) .lt. 1.01 .and. lambda(f).gt.0.99) then
           ! grep FAO_N2R camrun.o | uniq | awk '{print $3,$4,$5,$6,$7,$8}'
           !write(6,'(a,i3,e12.4,e12.4,e12.4,e12.4,e12.4,e12.4)'), 'FAO_N2R', &
           !     k, pmid(i,k), pflx(i,k), density, lambda(f), sigma(f), tau(i,k,f)
        endif

     enddo
     
  enddo
  enddo
  
  ! FAO: treat k=0 as a special case
  do i=1,Nday
     pmid_0 = 0.5*(pflx(i,1) + pflx(i,0))
     ! FAO:  kludge -- assume mmr(N2) = 1 - mmr(CH4)
     u = (1.0 - mmr(i,1)) * (pflx(i,1) - pflx(i,0)) / (1.0e2*gravit)   ! g / cm^2
     density = (0.1 * pmid_0) / (SHR_CONST_BOLTZ * tmid(i,1) ) * 1.0e-6  ! molecules / cm^3
     call rayleigh_n2(lambda, density, sigma, Nlambda)

     do f=1,Nlambda
        tau(i,0,f) = (sigma(f) / m_N2) * u;
     enddo
     
  enddo
!open(unit=99, file='gravmks.txt', status='unknown')
!write(99,*) gravit, pflx(1,0),pflx(1,1), pflx(1,pverp)
!close(unit=99)


end subroutine tau_rayleigh_n2


! FAO:  tiny helper function for debugging
subroutine str_of_int_less_than_1000(str,n)
  character*(*) str
  integer n

  if (n<10) then
     write(str,'(i1)') n
  else if (n<100) then
     write(str,'(i2)') n
  else if (n<1000) then
     write(str,'(i3)') n
  else
     str=''
  endif

end subroutine str_of_int_less_than_1000

!=========================================================================


subroutine radcswmx(lchnk   ,ncol, E_pint  ,E_pmid    , E_t, &
                    E_gmmr  ,E_ammr  ,E_cld     ,E_cicewp  ,E_cliqwp  ,E_rel     , &
                    E_rei     ,eccf    ,E_coszrs  ,scon    ,solin   , &
                    E_asdir   ,E_asdif   ,E_aldir   ,E_aldif   ,nmxrgn  , &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsnirtoa,fsnrtoac,fsnrtoaq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   ,frc_day ,                   &
                    aertau  ,aerssa  ,aerasm  ,aerfwd  ,fns     , &
                    fcns)
					
					!,    carmatau, carmassa, carmaasm, carmafwd)
!-----------------------------------------------------------------------
! 
! Purpose: 
! Solar radiation code
! 
! Method: 
! Basic method is Delta-Eddington as described in:
! 
! Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
! 
! Five changes to the basic method described above are:
! (1) addition of sulfate aerosols (Kiehl and Briegleb, 1993)
! (2) the distinction between liquid and ice particle clouds 
! (Kiehl et al, 1996);
! (3) provision for calculating TOA fluxes with spectral response to
! match Nimbus-7 visible/near-IR radiometers (Collins, 1998);
! (4) max-random overlap (Collins, 2001)
! (5) The near-IR absorption by H2O was updated in 2003 by Collins, 
!     Lee-Taylor, and Edwards for consistency with the new line data in
!     Hitran 2000 and the H2O continuum version CKD 2.4.  Modifications
!     were optimized by reducing RMS errors in heating rates relative
!     to a series of benchmark calculations for the 5 standard AFGL 
!     atmospheres.  The benchmarks were performed using DISORT2 combined
!     with GENLN3.  The near-IR scattering optical depths for Rayleigh
!     scattering were also adjusted, as well as the correction for
!     stratospheric heating by H2O.
!
! The treatment of maximum-random overlap is described in the
! comment block "INDEX CALCULATIONS FOR MAX OVERLAP".
! 
! Divides solar spectrum into 19 intervals from 0.2-5.0 micro-meters.
! solar flux fractions specified for each interval. allows for
! seasonally and diurnally varying solar input.  Includes molecular,
! cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud, 
! and surface absorption. Computes delta-eddington reflections and
! transmissions assuming homogeneously mixed layers. Adds the layers 
! assuming scattering between layers to be isotropic, and distinguishes 
! direct solar beam from scattered radiation.
! 
! Longitude loops are broken into 1 or 2 sections, so that only daylight
! (i.e. coszrs > 0) computations are done.
! 
! Note that an extra layer above the model top layer is added.
! 
! cgs units are used.
! 
! Special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
! 
!-----------------------------------------------------------------------
!
! D. Parks (NEC) 09/11/03
! Restructuring of routine to support SX vector architecture.
!
! Possible improvements:
!
! 1. Look at vectorizing index calculations for maximum overlap.
!
! 2. Consider making innermost loop in flux computations the number
!    of spectral intervals.  Given that NS is fixed at 19, the trade-off
!    will be stride one memory accesses of length 19 versus indirect
!    addressing (list vector - gather/scatter) with potential vector
!    lenghts of the number of day light points.  Vectorizing on the number
!    of spectral intervals seems worthwhile for low resolution models (T42),
!    but might be inefficient with higher resolutions.
!
! 3. Move the linearization of daylight points (compression/expansion) out
!    of radcswmx and into d_p_coupling.  This would eliminate the cost of
!    routines CmpDayNite and ExpDayNite.
!
! 4. Look at expliciting computing all streams in upward propagation of
!    radiation. There would be additional floating point operations in
!    exchange for the elimination of indirect addressing.
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use ghg_surfvals, only: ghg_surfvals_get_co2mmr
   use prescribed_aerosols, only: idxBG, idxSUL, idxSSLT, idxOCPHO, idxBCPHO, idxOCPHI, idxBCPHI, &
     idxDUSTfirst, numDUST, idxVOLC, naer_all
   use aer_optics, only: nrh, ndstsz, ksul, wsul, gsul, &
     ksslt, wsslt, gsslt, kcphil, wcphil, gcphil, kcphob, wcphob, gcphob, &
     kcb, wcb, gcb, kdst, wdst, gdst, kbg, wbg, gbg, kvolc, wvolc, gvolc
   use abortutils, only: endrun
   use cmparray_mod, only: CmpDayNite, ExpDayNite
   use quicksort, only: quick_sort
! EJL   
   use carma, only: carma_get_mmr, carma_optics_init, ncarma_bins, ncarma_elems, kcarma,wcarma,gcarma
   use history, only: outfld
   ! FAO
   use radcnst
   use pmgrid      , only : iam
   ! FAO: for debuging
   !use phys_grid,       only: get_rlat_all_p, get_rlon_all_p, get_lon_all_p, get_lat_all_p
   use phys_grid
   use mpishorthand

   use commap, only: latdeg

   implicit none

!-------------Parameters for accelerating max-random solution-------------
! 
! The solution time scales like prod(j:1->N) (1 + n_j) where 
! N   = number of max-overlap regions (nmxrgn)
! n_j = number of unique cloud amounts in region j
! 
! Therefore the solution cost can be reduced by decreasing n_j.
! cldmin reduces n_j by treating cloud amounts < cldmin as clear sky.
! cldeps reduces n_j by treating cloud amounts identical to log(1/cldeps)
! decimal places as identical
! 
! areamin reduces the cost by dropping configurations that occupy
! a surface area < areamin of the model grid box.  The surface area
! for a configuration C(j,k_j), where j is the region number and k_j is the
! index for a unique cloud amount (in descending order from biggest to
! smallest clouds) in region j, is
! 
! A = prod(j:1->N) [C(j,k_j) - C(j,k_j+1)]
! 
! where C(j,0) = 1.0 and C(j,n_j+1) = 0.0.
! 
! nconfgmax reduces the cost and improves load balancing by setting an upper
! bound on the number of cloud configurations in the solution.  If the number
! of configurations exceeds nconfgmax, the nconfgmax configurations with the
! largest area are retained, and the fluxes are normalized by the total area
! of these nconfgmax configurations.  For the current max/random overlap 
! assumption (see subroutine cldovrlap), 30 levels, and cloud-amount 
! parameterization, the mean and RMS number of configurations are 
! both roughly 5.  nconfgmax has been set to the mean+2*RMS number, or 15.
! 
! Minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
! 
   real(r8) cldmin
   parameter (cldmin = 1.0e-80_r8)
! 
! Minimimum horizontal area (as a fraction of the grid-box area) to retain 
! for a unique cloud configuration in the max-random solution
! 
   real(r8) areamin
   parameter (areamin = 0.01_r8)
! 
! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
! 
   real(r8) cldeps
   parameter (cldeps = 0.0_r8)
! 
! Maximum number of configurations to include in solution
! 
   integer nconfgmax
   parameter (nconfgmax = 15)
!------------------------------Commons----------------------------------
#include <crdcon.h>
! 
! Input arguments
! 
   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: ncol              ! number of atmospheric columns

   real(r8), intent(in) :: E_pmid(pcols,pver) ! Level pressure
   real(r8), intent(in) :: E_pint(pcols,pverp) ! Interface pressure
   real(r8), intent(in) :: E_t(pcols,pver)        ! Model level temperatures  ! FAO
   real(r8), intent(in) :: E_gmmr(pcols,pver,Ng) ! gas mass mixing ratio    ! FAO
   real(r8), intent(in) :: E_ammr(pcols,pver,Na) ! dust mass mixing ratio    ! FAO
! 
   real(r8), intent(in) :: E_cld(pcols,pver)  ! Fractional cloud cover
   real(r8), intent(in) :: E_cicewp(pcols,pver) ! in-cloud cloud ice water path
   real(r8), intent(in) :: E_cliqwp(pcols,pver) ! in-cloud cloud liquid water path
   real(r8), intent(in) :: E_rel(pcols,pver)  ! Liquid effective drop size (microns)
   real(r8), intent(in) :: E_rei(pcols,pver)  ! Ice effective drop size (microns)
! 
   real(r8), intent(in) :: eccf             ! Eccentricity factor (1./earth-sun dist^2)
   real(r8), intent(in) :: E_coszrs(pcols)    ! Cosine solar zenith angle
   real(r8), intent(in) :: E_asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: E_aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: E_asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8), intent(in) :: E_aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad

! FAO: old, no longer supported arguments:
   ! real(r8), intent(in) :: E_o3mmr(pcols,pver) ! Ozone mass mixing ratio
   ! real(r8), intent(in) :: E_h2ommr(pcols,pver) ! Specific humidity (h2o mass mix ratio)


   real(r8), intent(in) :: scon             ! solar constant 
! 
! IN/OUT arguments
! 
   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!                                                 !    maximally overlapped region. 
!                                                 !    0->pmxrgn(i,1) is range of pressure for
!                                                 !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                                 !    2nd region, etc
   integer, intent(inout) ::  nmxrgn(pcols)    ! Number of maximally overlapped regions
! 
! Output arguments
! 
   real(r8) :: ejlalb(pcols,Nf) ! albedo - EJL
   real(r8) :: ejlflxup(pcols,24) !Fluxes up and used for calculating albedo
   real(r8) :: ejlflxdn(pcols,24) !
   real(r8) :: alb325(pcols) 
   real(r8) :: alb375(pcols)
   real(r8) :: alb425(pcols)
   real(r8) :: alb475(pcols)
   real(r8) :: alb525(pcols)
   real(r8) :: alb575(pcols)
   real(r8) :: alb642(pcols)
   real(r8) :: alb714(pcols)
   real(r8) :: alb784(pcols)
   real(r8) :: alb845(pcols)
   real(r8) :: alb891(pcols)
   real(r8) :: alb940(pcols)
   real(r8) :: alb985(pcols)
   real(r8) :: alb1070(pcols)
   real(r8) :: alb1140(pcols)
   real(r8) :: alb1220(pcols)
   real(r8) :: alb1290(pcols)
   real(r8) :: alb1380(pcols)
   real(r8) :: alb1490(pcols)
   real(r8) :: alb1610(pcols)
   real(r8) :: alb1750(pcols)
   real(r8) :: alb1910(pcols)
   real(r8) :: alb2110(pcols)
   real(r8) :: alb2350(pcols)
   
   real(r8), intent(out) :: solin(pcols)     ! Incident solar flux
   real(r8), intent(out) :: qrs(pcols,pver)  ! Solar heating rate
   real(r8), intent(out) :: fsns(pcols)      ! Surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)      ! Total column absorbed solar flux
   real(r8), intent(out) :: fsntoa(pcols)    ! Net solar flux at TOA
   real(r8), intent(out) :: fsds(pcols)      ! Flux shortwave downwelling surface
! 
   real(r8), intent(out) :: fsnsc(pcols)     ! Clear sky surface absorbed solar flux
   real(r8), intent(out) :: fsdsc(pcols)     ! Clear sky surface downwelling solar flux
   real(r8), intent(out) :: fsntc(pcols)     ! Clear sky total column absorbed solar flx
   real(r8), intent(out) :: fsntoac(pcols)   ! Clear sky net solar flx at TOA
   real(r8), intent(out) :: sols(pcols)      ! Direct solar rad on surface (< 0.7)
   real(r8), intent(out) :: soll(pcols)      ! Direct solar rad on surface (>= 0.7)
   real(r8), intent(out) :: solsd(pcols)     ! Diffuse solar rad on surface (< 0.7)
   real(r8), intent(out) :: solld(pcols)     ! Diffuse solar rad on surface (>= 0.7)
   real(r8), intent(out) :: fsnirtoa(pcols)  ! Near-IR flux absorbed at toa
   real(r8), intent(out) :: fsnrtoac(pcols)  ! Clear sky near-IR flux absorbed at toa
   real(r8), intent(out) :: fsnrtoaq(pcols)  ! Net near-IR flux at toa >= 0.7 microns

   real(r8) , intent(out) :: frc_day(pcols) ! = 1 for daylight, =0 for night columns
   real(r8) :: aertau(pcols,Nf,ncarma_bins) ! Aerosol column optical depth - EJL
   real(r8) :: aerssa(pcols,Nf,ncarma_bins) ! Aerosol column averaged single scattering albedo
   real(r8) :: aerasm(pcols,Nf,ncarma_bins) ! Aerosol column averaged asymmetry parameter
   real(r8) :: aerfwd(pcols,Nf,ncarma_bins) ! Aerosol column averaged forward scattering
   real(r8), intent(out) :: fns(pcols,pverp)   ! net flux at interfaces
   real(r8), intent(out) :: fcns(pcols,pverp)  ! net clear-sky flux at interfaces


! 
!---------------------------Local variables-----------------------------
!
! Local and reordered copies of the intent(in) variables
!
   real(r8) :: pmid(pcols,pver) ! Level pressure
   real(r8) :: pint(pcols,pverp) ! Interface pressure
   real(r8) :: t(pcols,pver) ! midpt temperature               ! FAO
   real(r8) :: gmmr(pcols,pver,Ng) ! gas mass mixing ratio        ! FAO
   real(r8) :: ammr(pcols,pver,ncarma_bins) ! dust mass mixing ratio        ! FAO
! 
   real(r8) :: cld(pcols,pver)  ! Fractional cloud cover
   real(r8) :: cicewp(pcols,pver) ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,pver) ! in-cloud cloud liquid water path
   real(r8) :: rel(pcols,pver)  ! Liquid effective drop size (microns)
   real(r8) :: rei(pcols,pver)  ! Ice effective drop size (microns)
! 
   real(r8) :: coszrs(pcols)    ! Cosine solar zenith angle
   real(r8) :: asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8) :: aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8) :: asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8) :: aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad


! 
! Max/random overlap variables
! 
   real(r8) asort(pverp)     ! 1 - cloud amounts to be sorted for max ovrlp.
   real(r8) atmp             ! Temporary storage for sort when nxs = 2
   real(r8) cld0             ! 1 - (cld amt) used to make wstr, cstr, nstr
   real(r8) totwgt(pcols)    ! Total of xwgts = total fractional area of 
!   grid-box covered by cloud configurations
!   included in solution to fluxes

   real(r8) wgtv(nconfgmax)  ! Weights for fluxes
!   1st index is configuration number
   real(r8) wstr(pverp,pverp) ! area weighting factors for streams
!   1st index is for stream #, 
!   2nd index is for region #

   real(r8) xexpt            ! solar direct beam trans. for layer above
   real(r8) xrdnd            ! diffuse reflectivity for layer above
   real(r8) xrupd            ! diffuse reflectivity for layer below
   real(r8) xrups            ! direct-beam reflectivity for layer below
   real(r8) xtdnt            ! total trans for layers above

   real(r8) xwgt             ! product of cloud amounts

   real(r8) yexpt            ! solar direct beam trans. for layer above
   real(r8) yrdnd            ! diffuse reflectivity for layer above
   real(r8) yrupd            ! diffuse reflectivity for layer below
   real(r8) ytdnd            ! dif-beam transmission for layers above
   real(r8) ytupd            ! dif-beam transmission for layers below

   real(r8) zexpt            ! solar direct beam trans. for layer above
   real(r8) zrdnd            ! diffuse reflectivity for layer above
   real(r8) zrupd            ! diffuse reflectivity for layer below
   real(r8) zrups            ! direct-beam reflectivity for layer below
   real(r8) ztdnt            ! total trans for layers above

   logical new_term          ! Flag for configurations to include in fluxes
   logical region_found      ! flag for identifying regions

   integer ccon(nconfgmax,0:pverp,pcols)                                
! flags for presence of clouds
!   1st index is for level # (including 
!    layer above top of model and at surface)
!   2nd index is for configuration #
   integer cstr(0:pverp,pverp)                                
! flags for presence of clouds
!   1st index is for level # (including 
!    layer above top of model and at surface)
!   2nd index is for stream #
   integer icond(nconfgmax,0:pverp,pcols)
! Indices for copying rad. properties from
!     one identical downward cld config.
!     to another in adding method (step 2)
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer iconu(nconfgmax,0:pverp,pcols)
! Indices for copying rad. properties from
!     one identical upward configuration
!     to another in adding method (step 2)
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer iconfig           ! Counter for random-ovrlap configurations
   integer irgn              ! Index for max-overlap regions
   integer is0               ! Lower end of stream index range
   integer is1               ! Upper end of stream index range
   integer isn               ! Stream index
   integer istr(pverp+1)     ! index for stream #s during flux calculation
   integer istrtd(0:nconfgmax+1,0:pverp,pcols)
! indices into icond 
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer istrtu(0:nconfgmax+1,0:pverp,pcols)
! indices into iconu 
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer j                 ! Configuration index
   integer jj                ! Configuration index
   integer k1                ! Level index
   integer k2                ! Level index
   integer ksort(pverp)      ! Level indices of cloud amounts to be sorted
   integer ktmp              ! Temporary storage for sort when nxs = 2
   integer kx1(0:pverp)      ! Level index for top of max-overlap region
   integer kx2(0:pverp)      ! Level index for bottom of max-overlap region
   integer l                 ! Index 
   integer l0                ! Index
   integer mrgn              ! Counter for nrgn
   integer mstr              ! Counter for nstr
   integer n0                ! Number of configurations with ccon(:,k,:)==0
   integer n1                ! Number of configurations with ccon(:,k,:)==1
   integer nconfig(pcols)    ! Number of random-ovrlap configurations
   integer nconfigm          ! Value of config before testing for areamin,
!    nconfgmax
   integer npasses           ! number of passes over the indexing loop
   integer nrgn              ! Number of max overlap regions at current 
!    longitude
   integer nstr(pverp)       ! Number of unique cloud configurations
!   ("streams") in a max-overlapped region
!   1st index is for region #
   integer nuniq             ! # of unique cloud configurations
   integer nuniqd(0:pverp,pcols)   ! # of unique cloud configurations: TOA 
!   to level k
   integer nuniqu(0:pverp,pcols)   ! # of unique cloud configurations: surface
!   to level k 
   integer nxs               ! Number of cloudy layers between k1 and k2 
   integer ptr0(nconfgmax)   ! Indices of configurations with ccon(:,k,:)==0
   integer ptr1(nconfgmax)   ! Indices of configurations with ccon(:,k,:)==1
   integer ptrc(nconfgmax)   ! Pointer for configurations sorted by wgtv
   integer, dimension(1) :: min_idx  ! required for return val of func minloc

! 
! Other
! 
   integer ns                ! Spectral loop index
   integer i                 ! Longitude loop index
   integer k                 ! Level loop index
   integer km1               ! k - 1
   integer kp1               ! k + 1
   integer n                 ! Loop index for daylight
   integer indxsl            ! Index for cloud particle properties
   integer ksz               ! dust size bin index
   integer kaer              ! aerosol group index
   integer aaa

#ifdef TURN_ON_CLOUDS
! 
! A. Slingo's data for cloud particle radiative properties (from 'A GCM
! Parameterization for the Shortwave Properties of Water Clouds' JAS
! vol. 46 may 1989 pp 1419-1427)
! 
   real(r8) abarl(4)         ! A coefficient for extinction optical depth
   real(r8) bbarl(4)         ! B coefficient for extinction optical depth
   real(r8) cbarl(4)         ! C coefficient for single scat albedo
   real(r8) dbarl(4)         ! D coefficient for single  scat albedo
   real(r8) ebarl(4)         ! E coefficient for asymmetry parameter
   real(r8) fbarl(4)         ! F coefficient for asymmetry parameter

   save abarl, bbarl, cbarl, dbarl, ebarl, fbarl

   data abarl/ 2.817e-02, 2.682e-02,2.264e-02,1.281e-02/
   data bbarl/ 1.305    , 1.346    ,1.454    ,1.641    /
   data cbarl/-5.62e-08 ,-6.94e-06 ,4.64e-04 ,0.201    /
   data dbarl/ 1.63e-07 , 2.35e-05 ,1.24e-03 ,7.56e-03 /
   data ebarl/ 0.829    , 0.794    ,0.754    ,0.826    /
   data fbarl/ 2.482e-03, 4.226e-03,6.560e-03,4.353e-03/

   real(r8) abarli           ! A coefficient for current spectral band
   real(r8) bbarli           ! B coefficient for current spectral band
   real(r8) cbarli           ! C coefficient for current spectral band
   real(r8) dbarli           ! D coefficient for current spectral band
   real(r8) ebarli           ! E coefficient for current spectral band
   real(r8) fbarli           ! F coefficient for current spectral band
! 
! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
! greater than 20 micro-meters
! 
! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
! 
   real(r8) abari(4)         ! a coefficient for extinction optical depth
   real(r8) bbari(4)         ! b coefficient for extinction optical depth
   real(r8) cbari(4)         ! c coefficient for single scat albedo
   real(r8) dbari(4)         ! d coefficient for single scat albedo
   real(r8) ebari(4)         ! e coefficient for asymmetry parameter
   real(r8) fbari(4)         ! f coefficient for asymmetry parameter

   save abari, bbari, cbari, dbari, ebari, fbari

   data abari/ 3.448e-03, 3.448e-03,3.448e-03,3.448e-03/
   data bbari/ 2.431    , 2.431    ,2.431    ,2.431    /
   data cbari/ 1.00e-05 , 1.10e-04 ,1.861e-02,.46658   /
   data dbari/ 0.0      , 1.405e-05,8.328e-04,2.05e-05 /
   data ebari/ 0.7661   , 0.7730   ,0.794    ,0.9595   /
   data fbari/ 5.851e-04, 5.665e-04,7.267e-04,1.076e-04/

   real(r8) abarii           ! A coefficient for current spectral band
   real(r8) bbarii           ! B coefficient for current spectral band
   real(r8) cbarii           ! C coefficient for current spectral band
   real(r8) dbarii           ! D coefficient for current spectral band
   real(r8) ebarii           ! E coefficient for current spectral band
   real(r8) fbarii           ! F coefficient for current spectral band
! 
!   real(r8) delta            ! Pressure (in atm) for stratos. h2o limit
!   real(r8) o2mmr            ! O2 mass mixing ratio:

!   save delta, o2mmr

!
! UPDATE TO H2O NEAR-IR: Delta optimized for Hitran 2K and CKD 2.4
!
!   data delta / 0.0014257179260883 /
!
! END UPDATE
!
!   data o2mmr / .23143 /
#endif

   real(r8) albdir(pcols,Nf) ! Current spc intrvl srf alb to direct rad
   real(r8) albdif(pcols,Nf) ! Current spc intrvl srf alb to diffuse rad

! 
! Next series depends on spectral interval
! 
   real(r8) wgtint           ! Weight for specific spectral interval


! H2O ~ 1, O3 ~ 2, CO2 ~ 3, O2 ~ 4

! 
! Diagnostic and accumulation arrays; note that sfltot, fswup, and
! fswdn are not used in the computation,but are retained for future use.
! 
   real(r8) solflx(pcols)    ! Solar flux in current interval
   real(r8) sfltot(pcols)    ! Spectrally summed total solar flux
   real(r8) totfld(pcols,0:pver)   ! Spectrally summed flux divergence
   real(r8) fswup(pcols,0:pverp)   ! Spectrally summed up flux
   real(r8) fswdn(pcols,0:pverp)   ! Spectrally summed down flux

! 
! Cloud radiative property arrays
! 

   real(r8) tauxcl(pcols,0:pver) ! water cloud extinction optical depth
   real(r8) tauxci(pcols,0:pver) ! ice cloud extinction optical depth
   real(r8) wcl(pcols,0:pver) ! liquid cloud single scattering albedo
   real(r8) gcl(pcols,0:pver) ! liquid cloud asymmetry parameter
   real(r8) fcl(pcols,0:pver) ! liquid cloud forward scattered fraction
   real(r8) wci(pcols,0:pver) ! ice cloud single scattering albedo
   real(r8) gci(pcols,0:pver) ! ice cloud asymmetry parameter
   real(r8) fci(pcols,0:pver) ! ice cloud forward scattered fraction

!
! Aerosol mass paths by species
!
!  real(r8) usul(pcols,pver)   ! sulfate (SO4)
!  real(r8) ubg(pcols,pver)    ! background aerosol
!  real(r8) usslt(pcols,pver)  ! sea-salt (SSLT)
!  real(r8) ucphil(pcols,pver) ! hydrophilic organic carbon (OCPHI)
!  real(r8) ucphob(pcols,pver) ! hydrophobic organic carbon (OCPHO)
!  real(r8) ucb(pcols,pver)    ! black carbon (BCPHI + BCPHO)
!  real(r8) uvolc(pcols,pver) ! volcanic mass
!  real(r8) udst(pcols,ndstsz,pver) ! dust


!
! local variables used for the external mixing of aerosol species
!
  real(r8) tau_dst_tot         ! optical depth, total dust
  real(r8) tau_w_dst_tot       ! optical depth * single scattering albedo, total dust
  real(r8) tau_w_g_dst_tot     ! optical depth * single scattering albedo * asymmetry parameter, total dust

  real(r8) f_sul               ! forward scattering fraction, sulfate
  real(r8) f_bg                ! forward scattering fraction, background aerosol
  real(r8) f_sslt              ! forward scattering fraction, sea-salt
  real(r8) f_cphil             ! forward scattering fraction, hydrophilic carbon
  real(r8) f_cphob             ! forward scattering fraction, hydrophobic carbon
  real(r8) f_cb                ! forward scattering fraction, black carbon
  real(r8) f_volc              ! forward scattering fraction, volcanic
  real(r8) f_dst(ndstsz)       ! forward scattering fraction, dust, by size
  real(r8) f_dst_tot           ! forward scattering fraction, total dust

  real(r8) tau_w_f_sul         ! optical depth * forward scattering fraction * single scattering albedo, sulfate
  real(r8) tau_w_f_bg          ! optical depth * forward scattering fraction * single scattering albedo, background
  real(r8) tau_w_f_sslt        ! optical depth * forward scattering fraction * single scattering albedo, sea-salt
  real(r8) tau_w_f_cphil       ! optical depth * forward scattering fraction * single scattering albedo, hydrophilic C
  real(r8) tau_w_f_cphob       ! optical depth * forward scattering fraction * single scattering albedo, hydrophobic C
  real(r8) tau_w_f_cb          ! optical depth * forward scattering fraction * single scattering albedo, black C
  real(r8) tau_w_f_volc        ! optical depth * forward scattering fraction * single scattering albedo, volcanic
  real(r8) tau_w_f_dst(ndstsz) ! optical depth * forward scattering fraction * single scattering albedo, dust, by size
  real(r8) tau_w_f_dst_tot     ! optical depth * forward scattering fraction * single scattering albedo, total dust

  real(r8) tau_tot             ! optical depth, total aerosol
  real(r8) tau_w_tot           ! optical depth * single scattering albedo, total aerosol
  real(r8) tau_w_g_tot         ! optical depth * single scattering albedo * asymmetry parameter, total aerosol
  real(r8) f_tot               ! forward scattering fraction, total aerosol
  real(r8) tau_w_f_tot         ! optical depth * forward scattering fraction * single scattering albedo, total aerosol

  real(r8) w_dst_tot           ! single scattering albedo, total dust
  real(r8) w_tot               ! single scattering albedo, total aerosol
  real(r8) g_dst_tot           ! asymmetry parameter, total dust
  real(r8) g_tot               ! asymmetry parameter, total aerosol
  real(r8) ksuli               ! specific extinction interpolated between rh look-up-table points, sulfate
  real(r8) ksslti              ! specific extinction interpolated between rh look-up-table points, sea-salt
  real(r8) kcphili             ! specific extinction interpolated between rh look-up-table points, hydrophilic carbon
  real(r8) wsuli               ! single scattering albedo interpolated between rh look-up-table points, sulfate
  real(r8) wsslti              ! single scattering albedo interpolated between rh look-up-table points, sea-salt
  real(r8) wcphili             ! single scattering albedo interpolated between rh look-up-table points, hydrophilic carbon
  real(r8) gsuli               ! asymmetry parameter interpolated between rh look-up-table points, sulfate
  real(r8) gsslti              ! asymmetry parameter interpolated between rh look-up-table points, sea-salt
  real(r8) gcphili             ! asymmetry parameter interpolated between rh look-up-table points, hydrophilic carbon

! 
! Aerosol radiative property arrays
! 
   real(r8) tauxar(pcols,0:pver) ! aerosol extinction optical depth
   real(r8) wa(pcols,0:pver) ! aerosol single scattering albedo
   real(r8) ga(pcols,0:pver) ! aerosol assymetry parameter
   real(r8) fa(pcols,0:pver) ! aerosol forward scattered fraction

! 
! Various arrays and other constants:
! 
   real(r8) pflx(pcols,0:pverp) ! Interface press, including extra layer
   real(r8) pmid_copy(pcols,0:pver) ! midlevel pressure, including extra layer
   real(r8) zenfac(pcols)    ! Square root of cos solar zenith angle
   real(r8) sqrco2           ! Square root of the co2 mass mixg ratio
   real(r8) tmp1             ! Temporary constant array
   real(r8) tmp2             ! Temporary constant array
   real(r8) pdel             ! Pressure difference across layer
   real(r8) path             ! Mass path of layer
   real(r8) ptop             ! Lower interface pressure of extra layer
   real(r8) ptho2            ! Used to compute mass path of o2
   real(r8) ptho3            ! Used to compute mass path of o3
   real(r8) pthco2           ! Used to compute mass path of co2
   real(r8) pthh2o           ! Used to compute mass path of h2o
   real(r8) h2ostr           ! Inverse sq. root h2o mass mixing ratio

#ifdef JPE_VMATH
   real(r8) v_h2ostr(pcols,pver)           ! Inverse sq. root h2o mass mixing ratio

   real(r8) v_rtotwgt(pcols)           ! recipricle totwgt
#endif

   real(r8) wavmid(Nf)   ! Spectral interval middle wavelength
   real(r8) tmp1l            ! Temporary constant array
   real(r8) tmp2l            ! Temporary constant array
   real(r8) tmp3l            ! Temporary constant array
   real(r8) tmp1i            ! Temporary constant array
   real(r8) tmp2i            ! Temporary constant array
   real(r8) tmp3i            ! Temporary constant array
   real(r8) rdenom           ! Multiple scattering term
   real(r8) rdirexp          ! layer direct ref times exp transmission
   real(r8) tdnmexp          ! total transmission - exp transmission
   real(r8) psf(Nf)      ! Frac of solar flux in spect interval
! 
! Layer absorber amounts; note that 0 refers to the extra layer added
! above the top model layer
! 
   real(r8) uh2o(pcols,0:pver) ! Layer absorber amount of h2o
   real(r8) uo3(pcols,0:pver) ! Layer absorber amount of  o3
   real(r8) uco2(pcols,0:pver) ! Layer absorber amount of co2
   real(r8) uo2(pcols,0:pver) ! Layer absorber amount of  o2
   real(r8) uaer(pcols,0:pver) ! Layer aerosol amount 

!  FAO
   integer :: Za, Zg
   real(r8) abund_g(pcols,0:pver,Ng)         ! Layer absorber amount
   real(r8) abund_a(pcols,0:pver,ncarma_bins)         ! Layer absorber amount
   real(r8) tau_sum
   real(r8) tau_w_sum
   real(r8) tau_w_g_sum
   real(r8) tau_w_f_sum
   real(r8) tau_a, tau_w, tau_w_g, tau_w_f
   real(r8) taugab(pcols,0:pver)        ! Layer total gas absorption optical depth


! 
! Total column absorber amounts:
! 
   real(r8) uth2o(pcols)     ! Total column  absorber amount of  h2o
   real(r8) uto3(pcols)      ! Total column  absorber amount of  o3
   real(r8) utco2(pcols)     ! Total column  absorber amount of  co2
   real(r8) uto2(pcols)      ! Total column  absorber amount of  o2

! 
! These arrays are defined for pver model layers; 0 refers to the extra
! layer on top:
! 
   real(r8) rdir(Nf,pcols,0:pver) ! Layer reflectivity to direct rad
   real(r8) rdif(Nf,pcols,0:pver) ! Layer reflectivity to diffuse rad
   real(r8) tdir(Nf,pcols,0:pver) ! Layer transmission to direct rad
   real(r8) tdif(Nf,pcols,0:pver) ! Layer transmission to diffuse rad
   real(r8) explay(Nf,pcols,0:pver) ! Solar beam exp trans. for layer

   real(r8) rdirc(Nf,pcols,0:pver) ! Clear Layer reflec. to direct rad
   real(r8) rdifc(Nf,pcols,0:pver) ! Clear Layer reflec. to diffuse rad
   real(r8) tdirc(Nf,pcols,0:pver) ! Clear Layer trans. to direct rad
   real(r8) tdifc(Nf,pcols,0:pver) ! Clear Layer trans. to diffuse rad
   real(r8) explayc(Nf,pcols,0:pver) ! Solar beam exp trans. clear layer

!  FAO
   real(r8) rdir_0(Nf,pcols,0:pver) ! Layer reflectivity to direct rad
   real(r8) rdif_0(Nf,pcols,0:pver) ! Layer reflectivity to diffuse rad
   real(r8) tdir_0(Nf,pcols,0:pver) ! Layer transmission to direct rad
   real(r8) tdif_0(Nf,pcols,0:pver) ! Layer transmission to diffuse rad
   real(r8) explay_0(Nf,pcols,0:pver) ! Solar beam exp trans. for layer

   real(r8) rdirc_0(Nf,pcols,0:pver) ! Clear Layer reflec. to direct rad
   real(r8) rdifc_0(Nf,pcols,0:pver) ! Clear Layer reflec. to diffuse rad
   real(r8) tdirc_0(Nf,pcols,0:pver) ! Clear Layer trans. to direct rad
   real(r8) tdifc_0(Nf,pcols,0:pver) ! Clear Layer trans. to diffuse rad
   real(r8) explayc_0(Nf,pcols,0:pver) ! Solar beam exp trans. clear layer

   real(r8) flxdiv           ! Flux divergence for layer

!
! Temporary arrays for either clear or cloudy values.
!
   real(r8), dimension(Nf) :: Trdir
   real(r8), dimension(Nf) :: Trdif
   real(r8), dimension(Nf) :: Ttdir
   real(r8), dimension(Nf) :: Ttdif
   real(r8), dimension(Nf) :: Texplay
!cdir vreg(Trdir)
!cdir vreg(Trdif)
!cdir vreg(Ttdir)
!cdir vreg(Ttdif)
!cdir vreg(Texplay)
! 
! 
! Radiative Properties:
! 
! There are 1 classes of properties:
! (1. All-sky bulk properties
! (2. Clear-sky properties
! 
! The first set of properties are generated during step 2 of the solution.
! 
! These arrays are defined at model interfaces; in 1st index (for level #),
! 0 is the top of the extra layer above the model top, and
! pverp is the earth surface.  2nd index is for cloud configuration
! defined over a whole column.
! 
   real(r8) exptdn(Nf,0:pverp,nconfgmax,pcols) ! Sol. beam trans from layers above
   real(r8) rdndif(Nf,0:pverp,nconfgmax,pcols) ! Ref to dif rad for layers above
   real(r8) rupdif(Nf,0:pverp,nconfgmax,pcols) ! Ref to dif rad for layers below
   real(r8) rupdir(Nf,0:pverp,nconfgmax,pcols) ! Ref to dir rad for layers below
   real(r8) tdntot(Nf,0:pverp,nconfgmax,pcols) ! Total trans for layers above
! 
! Bulk properties used during the clear-sky calculation.
! 
   real(r8) exptdnc(pcols,0:pverp) ! clr: Sol. beam trans from layers above
   real(r8) rdndifc(pcols,0:pverp) ! clr: Ref to dif rad for layers above
   real(r8) rupdifc(pcols,0:pverp) ! clr: Ref to dif rad for layers below
   real(r8) rupdirc(pcols,0:pverp) ! clr: Ref to dir rad for layers below
   real(r8) tdntotc(pcols,0:pverp) ! clr: Total trans for layers above

   real(r8) fluxup(Nf,0:pverp,pcols)  ! Up   flux at model interface
   real(r8) fluxdn(Nf,0:pverp,pcols)  ! Down flux at model interface
   real(r8) wexptdn(Nf,pcols)   ! Direct solar beam trans. to surface
!
! Scalars used in vectorization
!
  integer :: kk
  integer :: Nday                      ! Number of daylight columns
  integer :: Nnite                     ! Number of night columns
  integer, dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
  integer, dimension(pcols) :: IdxNite ! Indicies of night coumns
!
! Arrays used in vectorization
!
  real(r8), dimension(pcols,ndstsz) :: v_tau_dst             ! optical depth, dust, by size category
  real(r8), dimension(pcols,ndstsz) :: v_tau_w_dst           ! optical depth * single scattering albedo, dust, by size
  real(r8), dimension(pcols,ndstsz) :: v_tau_w_g_dst         ! optical depth * single scattering albedo * asymmetry parameter, dust, by size
  real(r8), dimension(pcols,ndstsz) :: v_tau_w_f_dst         ! optical depth * forward scattering fraction * single scattering albedo, dust, by size
  real(r8), dimension(pcols)        :: v_tau_dst_tot         ! optical depth, total dust
  real(r8), dimension(pcols)        :: v_tau_w_dst_tot       ! optical depth * single scattering albedo, total dust
  real(r8), dimension(pcols)        :: v_tau_w_g_dst_tot     ! optical depth * single scattering albedo * asymmetry parameter, total dust
  real(r8), dimension(pcols)        :: v_tau_w_f_dst_tot     ! optical depth * forward scattering fraction * single scattering albedo, total dust

   real(r8) :: Tv_tau_dst
   real(r8) :: Tv_tau_w_dst
   real(r8) :: Tv_tau_w_g_dst
   real(r8) :: Tv_tau_w_f_dst

   real(r8) v_wgtv(nconfgmax,pcols)  ! Weights for fluxes

   logical :: lg_tot_gt1(pver)        ! Logical g_tot > 1.0_r8
   logical :: lg_tot_ltm1(pver)       ! Logical g_tot < -1.0_r8
   logical :: lf_tot_gt1(pver)        ! Logical f_tot > 1.0_r8
   logical :: lf_tot_lt0(pver)        ! Logical f_tot < 0.0_r8

!   real(r8) :: rdiff, ro, rn
!   rdiff(ro,rn) = abs((ro-rn)/merge(ro,1.0_r8,ro /= 0.0_r8))

   real(r8) :: tau_rayleigh(pcols,0:pver,Nf)

   ! FAO: for debugging
   integer ierr, Zp
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)
   real(r8) :: clat_day(pcols)                   ! current latitudes(radians)
   real(r8) :: clon_day(pcols)                   ! current longitudes(radians)
   character(len=4) :: iam_str
   real(r8) :: ssa_scaled

!#define EXAMPLE_MPI
#ifdef EXAMPLE_MPI
   integer :: Ilat(pcols)                   ! current latitudes(index)
   integer :: Ilon(pcols)                   ! current longitudes(index)
   real(r8) :: snd(plat, pverp+1), rcv(plat, pverp+1), x_sum(plat,pverp+1)

   save x_sum
#endif

   ! FAO: for debugging
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)
   call str_of_int_less_than_1000(iam_str,iam)


#ifdef EXAMPLE_MPI
   call get_lat_all_p(lchnk, ncol, Ilat)
   call get_lon_all_p(lchnk, ncol, Ilon)
#endif



! 
!-----------------------------------------------------------------------
! START OF CALCULATION
!-----------------------------------------------------------------------
! 
!  write (6, '(a, x, i3)') 'radcswmx : chunk identifier', lchnk

! 
! Initialize output fields:
! 
   fsds(1:ncol)     = 0.0_r8

   fsnirtoa(1:ncol) = 0.0_r8
   fsnrtoac(1:ncol) = 0.0_r8
   fsnrtoaq(1:ncol) = 0.0_r8

   fsns(1:ncol)     = 0.0_r8
   fsnsc(1:ncol)    = 0.0_r8
   fsdsc(1:ncol)    = 0.0_r8

   fsnt(1:ncol)     = 0.0_r8
   fsntc(1:ncol)    = 0.0_r8
   fsntoa(1:ncol)   = 0.0_r8
   fsntoac(1:ncol)  = 0.0_r8

   solin(1:ncol)    = 0.0_r8

   sols(1:ncol)     = 0.0_r8
   soll(1:ncol)     = 0.0_r8
   solsd(1:ncol)    = 0.0_r8
   solld(1:ncol)    = 0.0_r8

   frc_day(1:ncol) = 0.0_r8

   ejlalb(1:ncol,1:Nf) = 0.0_r8
   qrs(1:ncol,1:pver) = 0.0_r8
   fns(1:ncol,1:pverp) = 0.0_r8
   fcns(1:ncol,1:pverp) = 0.0_r8

   ! initialize aerosol diagnostic fields to 0.0 
   ! Average can be obtained by dividing <aerod>/<frc_day>

   aertau(1:ncol,1:Nf,1:ncarma_bins) = 0.0_r8
   aerssa(1:ncol,1:Nf,1:ncarma_bins) = 0.0_r8
   aerasm(1:ncol,1:Nf,1:ncarma_bins) = 0.0_r8
   aerfwd(1:ncol,1:Nf,1:ncarma_bins) = 0.0_r8

   ! FAO
   rdir = 0
   rdif = 0
   tdir = 0
   tdif = 0
   explay = 0
   rdirc = 0
   rdifc = 0
   tdirc = 0
   tdifc = 0
   explayc = 0

! 
! Compute starting, ending daytime loop indices:
!  *** Note this logic assumes day and night points are contiguous so
!  *** will not work in general with chunked data structure.
! 

   Nday = 0
   Nnite = 0
   do i = 1, ncol
      if ( E_coszrs(i) > 0.0_r8 ) then
         Nday = Nday + 1
         IdxDay(Nday) = i
         !print*, 'FAO_DAY', clat(i), clon(i), E_coszrs(i)
      else
         Nnite = Nnite + 1
         IdxNite(Nnite) = i
         !print*, 'FAO_NIT', clat(i), clon(i), E_coszrs(i)
      end if
   end do

! 
! If night everywhere, return:
! 
   if ( Nday == 0 ) return
!
! Rearrange input arrays
!
   call CmpDayNite(E_pmid, pmid,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_pint, pint,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call CmpDayNite(E_t,    t, 	        Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_gmmr, gmmr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver, 1, Ng)    ! FAO
   call CmpDayNite(E_ammr, ammr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver, 1, ncarma_bins)    ! FAO
   call CmpDayNite(E_cld, cld,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_cicewp, cicewp,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_cliqwp, cliqwp,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_rel, rel, 		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_rei, rei,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_coszrs, coszrs,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_asdir, asdir,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aldir, aldir,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_asdif, asdif,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aldif, aldif,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)

   call CmpDayNite(pmxrgn,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call CmpDayNite(nmxrgn,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)

   !FAO: for debugging
   call CmpDayNite(clat, clat_day,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(clon, clon_day,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)

#ifdef EXAMPLE_MPI
   call CmpDayNite(Ilat, Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(Ilon, Nday, IdxDay, Nnite, IdxNite, 1, pcols)
#endif

! 
! Perform other initializations
! 
   do k=1,pverp
      do i=1,Nday
         pflx(i,k) = pint(i,k)
      enddo
   enddo
   ! NB: pressure is set to zero, at interface above the top one
   do i=1,Nday
      pflx(i,0) = 0._r8  
   enddo

   ! FAO: this bit of code has been added, b/c we need pmid(k=0) later on
   do k=1,pver
      do i=1,Nday
         pmid_copy(i,k) = pmid(i,k)
      enddo
   enddo
   do i=1,Nday
      pmid_copy(i,0) = 0.5*(pflx(i,0) + pflx(i,1))
   enddo

!!!!!
!!!!! FAO: old code that has been replaced
!!!!!

!   tmp1   = 0.5_r8/(gravit*sslp)
!   tmp2   = delta/gravit

!   sqrco2 = sqrt(ghg_surfvals_get_co2mmr())

!#ifdef JPE_VMATH
!   call vrsqrt(v_h2ostr,h2ommr,pver*pcols)
!#endif


!   do i=1,Nday
!! 
!! Define solar incident radiation and interface pressures:
!! 
!         solin(i)  = scon*eccf*coszrs(i)
!         pflx(i,0) = 0._r8
!! 
!! Compute optical paths:
!! 
!         ptop      = pflx(i,1)
!         ptho2     = o2mmr * ptop / gravit
!         ptho3     = o3mmr(i,1) * ptop / gravit
!         pthco2    = sqrco2 * (ptop / gravit)
!#ifdef JPE_VMATH
!         zenfac(i) = sqrt(coszrs(i))
!         pthh2o    = ptop**2*tmp1 + (ptop*rga)* &
!                    (v_h2ostr(i,1)*zenfac(i)*delta)
!#else
!         h2ostr    = sqrt( 1._r8 / h2ommr(i,1) )
!         zenfac(i) = sqrt(coszrs(i))
!         pthh2o    = ptop**2*tmp1 + (ptop*rga)* &
!                    (h2ostr*zenfac(i)*delta)
!#endif
!         uh2o(i,0) = h2ommr(i,1)*pthh2o
!         uco2(i,0) = zenfac(i)*pthco2
!         uo2 (i,0) = zenfac(i)*ptho2
!         uo3 (i,0) = ptho3
!         uaer(i,0) = 0.0_r8
!! 
!! End  do i=1,Nday
!! 
!   end do
!
!   do k=1,pver
!!cdir nodep
!      do i=1,Nday
!         pdel      = pflx(i,k+1) - pflx(i,k)
!         path      = pdel / gravit
!         ptho2     = o2mmr * path
!         ptho3     = o3mmr(i,k) * path
!         pthco2    = sqrco2 * path
!         h2ostr    = sqrt(1.0_r8/h2ommr(i,k))
!         pthh2o    = (pflx(i,k+1)**2 - pflx(i,k)**2)*tmp1 + pdel*h2ostr*zenfac(i)*tmp2
!         uh2o(i,k) = h2ommr(i,k)*pthh2o
!         uco2(i,k) = zenfac(i)*pthco2
!         uo2 (i,k) = zenfac(i)*ptho2
!         uo3 (i,k) = ptho3
!         usul(i,k) = aermmr(i,k,idxSUL) * path 
!         ubg(i,k) = aermmr(i,k,idxBG) * path 
!         usslt(i,k) = aermmr(i,k,idxSSLT) * path
!         if (usslt(i,k) .lt. 0.0) then  ! usslt is sometimes small and negative, will be fixed
!           usslt(i,k) = 0.0
!         end if
!         ucphil(i,k) = aermmr(i,k,idxOCPHI) * path
!         ucphob(i,k) = aermmr(i,k,idxOCPHO) * path
!         ucb(i,k) = ( aermmr(i,k,idxBCPHO) + aermmr(i,k,idxBCPHI) ) * path
!         uvolc(i,k) =  aermmr(i,k,idxVOLC)
!!cdir expand=ndstsz
!         do ksz = 1, ndstsz
!           udst(i,ksz,k) = aermmr(i,k,idxDUSTfirst-1+ksz) * path
!         end do
!      end do
!   end do
!
!! FAO:  total column absorber amounts do not seem to be used, so delete
!! 
!! Compute column absorber amounts for the clear sky computation:
!! 
!   do i=1,Nday
!
!      uth2o(i) = 0.0_r8
!      uto3(i)  = 0.0_r8
!      utco2(i) = 0.0_r8
!      uto2(i)  = 0.0_r8
!
!!cdir expand=pver
!      do k=1,pver
!         uth2o(i) = uth2o(i) + uh2o(i,k)
!         uto3(i)  = uto3(i)  + uo3(i,k)
!         utco2(i) = utco2(i) + uco2(i,k)
!         uto2(i)  = uto2(i)  + uo2(i,k)
!      end do
!   end do


   ! FAO:  replace individual u's with abundances (units of g/cm^2)
   do i=1,Nday
      ! Define solar incident radiation and interface pressures:
      solin(i)  = scon*eccf*coszrs(i)


      ! NB: gravit used here is in cgs units 
      path = (pflx(i,1) - pflx(i,0)) / gravit
      do Zg=1,Ng
         abund_g(i,0,Zg) = gmmr(i,1,Zg) * path  ! layers 0,1 assumed to have same mmr
		 
!		 open(unit=99,file='gmmr.txt',status='unknown')
!         write(99,*) i,Ng,gmmr(i,:,Ng)
!         close(unit=99)
	  enddo
      do Za=1,ncarma_bins
         abund_a(i,0,Za) = 0.0_r8
      enddo
   enddo

   do k=1,pver
   do i=1,Nday
      path = (pflx(i,k+1) - pflx(i,k)) / gravit
      do Zg=1,Ng
         abund_g(i,k,Zg) = gmmr(i,k,Zg) * path
      enddo
      do Za=1,ncarma_bins
         abund_a(i,k,Za) = ammr(i,k,Za) * path ! units g/cm2
      enddo
   enddo
   enddo

!open(unit=99, file='gravcgs.txt', status='unknown')
!write(99,*) gravit, pflx(1,0), pflx(1,1), pflx(1,pver)
!close(unit=99)

   do ns=1,Nf
      wavmid(ns) = 0.5_r8*(lambda_min(ns) + lambda_max(ns)) ! midpoint wavelength
   enddo

   ! FAO
   call tau_rayleigh_n2(tau_rayleigh, wavmid, gmmr, pmid_copy, pflx, t, Nday, Nf)

!#define FAO_DO_NOT_INCLUDE_RAYLEIGH
#ifdef FAO_DO_NOT_INCLUDE_RAYLEIGH
   tau_rayleigh = 1d-10
#endif

! 
! Set cloud properties for top (0) layer; so long as tauxcl is zero,
! there is no cloud above top of model; the other cloud properties
! are arbitrary:
! 
      do i=1,Nday

         tauxcl(i,0)  = 0._r8
         wcl(i,0)     = 0.999999_r8
         gcl(i,0)     = 0.85_r8
         fcl(i,0)     = 0.725_r8
         tauxci(i,0)  = 0._r8
         wci(i,0)     = 0.999999_r8
         gci(i,0)     = 0.85_r8
         fci(i,0)     = 0.725_r8
! 
! Aerosol 
! 
         tauxar(i,0)  = 0._r8
         wa(i,0)      = 0.925_r8
         ga(i,0)      = 0.850_r8
         fa(i,0)      = 0.7225_r8

      end do ! End  do i=1,Nday


   ! FAO:  not sure this is necessary, but initialize all cloud properties to zero
   tauxcl  = 0._r8
   wcl     = 0.999999_r8
   gcl     = 0._r8
   fcl     = 0._r8
   tauxci  = 0._r8
   wci    = 0.999999_r8
   gci    = 0._r8
   fci    = 0._r8

!EJL
!			   open(unit=99, file='optical_properties.txt', status='unknown')
!               write(99,*) aerssa(1,1,1), aertau(1,1,1)
!			   write(99,*) kcarma
!			   write(99,*) ''
!			   write(99,*) wcarma, gcarma
!			   close(unit=99)


! 
! Begin spectral loop
! 
   do ns=1,Nf

#ifdef TURN_ON_CLOUDS
! 
! Set index for cloud particle properties based on the wavelength,
! according to A. Slingo (1989) equations 1-3:
! Use index 1 (0.25 to 0.69 micrometers) for visible
! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
! 
! Note that the minimum wavelength is encoded (with .001, .002, .003)
! in order to specify the index appropriate for the near-infrared
! cloud absorption properties
! 
      if(lambda_max(ns) <= 0.7_r8) then
         indxsl = 1
      else if(lambda_min(ns) == 0.700_r8) then
         indxsl = 2
      else if(lambda_min(ns) == 0.701_r8) then
         indxsl = 3
      else if(lambda_min(ns) == 0.702_r8 .or. lambda_min(ns) > 2.38_r8) then
         indxsl = 4
      end if
! 
! Set cloud extinction optical depth, single scatter albedo,
! asymmetry parameter, and forward scattered fraction:
! 
      abarli = abarl(indxsl)
      bbarli = bbarl(indxsl)
      cbarli = cbarl(indxsl)
      dbarli = dbarl(indxsl)
      ebarli = ebarl(indxsl)
      fbarli = fbarl(indxsl)
! 
      abarii = abari(indxsl)
      bbarii = bbari(indxsl)
      cbarii = cbari(indxsl)
      dbarii = dbari(indxsl)
      ebarii = ebari(indxsl)
      fbarii = fbari(indxsl)
#endif ! TURN_ON_CLOUDS

! 
! adjust fraction within spectral interval to allow for the possibility of
! sub-divisions within a particular interval:
! 
      psf(ns) = 1.0_r8
      do Zg=1,Ng
         if (spectralweight(ns,Zg)/=0._r8) then
            psf(ns) = psf(ns)*spectralweight(ns,Zg)
         endif
      enddo

      do i=1,Nday
         frc_day(i) = 1.0_r8
      end do ! End do i=1,Nday

!      f_cphob  = gcphob(ns) * gcphob(ns)
!      f_cb     = gcb(ns) * gcb(ns)
!      f_volc   = gvolc(ns) * gvolc(ns)
!      f_bg     = gbg(ns) * gbg(ns)
!      f_dst(:) = gdst(:,ns) * gdst(:,ns)

!CSD$ PARALLEL DO PRIVATE( v_tau_dst_tot ) &
!CSD$ PRIVATE( v_tau_w_dst_tot, v_tau_w_g_dst_tot, v_tau_w_f_dst_tot, Tv_tau_dst, Tv_tau_w_dst, Tv_tau_w_g_dst ) &
!CSD$ PRIVATE( Tv_tau_w_f_dst, tmp1l, tmp2l, tmp3l, tmp1i, tmp2i, tmp3i, rhtrunc, krh, wrh, ksuli, ksslti ) &
!CSD$ PRIVATE( kcphili, wsuli, wsslti, wcphili, gsuli, gsslti, gcphili, tau_sul, tau_sslt, tau_cphil, tau_cphob ) &
!CSD$ PRIVATE( tau_cb, tau_volc, tau_bg, tau_w_sul, tau_w_sslt, tau_w_cphil, tau_w_cphob, tau_w_cb, tau_w_volc, tau_w_bg ) &
!CSD$ PRIVATE( tau_w_g_sul, tau_w_g_sslt, tau_w_g_cphil, tau_w_g_cphob, tau_w_g_cb, tau_w_g_volc, tau_w_g_bg, f_sul,f_sslt ) &
!CSD$ PRIVATE(  f_cphil, tau_w_f_sul, tau_w_f_bg, tau_w_f_sslt, tau_w_f_cphil, tau_w_f_cphob, tau_w_f_cb, tau_w_f_volc ) &
!CSD$ PRIVATE( tau_dst_tot, tau_w_dst_tot, tau_w_g_dst_tot, tau_w_f_dst_tot, w_dst_tot, g_dst_tot, f_dst_tot, tau_tot ) &
!CSD$ PRIVATE( tau_w_tot, tau_w_g_tot, tau_w_f_tot, w_tot, g_tot, f_tot, k, i, kk )

      do k=1,pver

!
! The following logicals are used to check to see whether we have invalid
! values for rhtrunc, g_tot, and f_tot.
! The above values are checked for every level of K and column of I.
! However, since the out of bounds conditions are *extremely* rare (in fact
! we abort the run if we find any one), let's reduce the number of 'if' tests
! by 1/K, and only check for aborting the run after scanning all Ks and Is.
!
         lg_tot_gt1(k)   = .false.
         lg_tot_ltm1(k)  = .false.
         lf_tot_gt1(k)   = .false.
         lf_tot_lt0(k)   = .false.
         
#ifdef UNUSED_CODE
         ! Note(FAO):  This code sets v_tau_*_tot, which sets tau_*_dst_tot, 
         !             which sets *_dst_tot, which is not used
         v_tau_dst_tot(1:Nday)     = 0.0_r8
         v_tau_w_dst_tot(1:Nday)   = 0.0_r8
         v_tau_w_g_dst_tot(1:Nday) = 0.0_r8
         v_tau_w_f_dst_tot(1:Nday) = 0.0_r8
         do i=1,Nday
!cdir expand
            do kk = 1, ndstsz
               Tv_tau_dst           = 1.e4 * kdst(kk,ns) * udst(i,kk,k)
               Tv_tau_w_dst         = Tv_tau_dst         * wdst(kk,ns)
               Tv_tau_w_g_dst       = Tv_tau_w_dst       * gdst(kk,ns)
               Tv_tau_w_f_dst       = Tv_tau_w_dst       * f_dst(kk)
               v_tau_dst_tot(i)     = v_tau_dst_tot(i)     + Tv_tau_dst
               v_tau_w_dst_tot(i)   = v_tau_w_dst_tot(i)   + Tv_tau_w_dst
               v_tau_w_g_dst_tot(i) = v_tau_w_g_dst_tot(i) + Tv_tau_w_g_dst
               v_tau_w_f_dst_tot(i) = v_tau_w_f_dst_tot(i) + Tv_tau_w_f_dst
            end do
         end do ! End do i=1,Nday
#endif ! UNUSED_CODE

         do i=1,Nday

#ifdef TURN_ON_CLOUDS
! 
! liquid
! 

               tmp2l = 1._r8 - cbarli - dbarli*rel(i,k)
               tmp3l = fbarli*rel(i,k)
! 
! ice
! 

               tmp2i = 1._r8 - cbarii - dbarii*rei(i,k)
               tmp3i = fbarii*rei(i,k)

               if (cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
! liquid
                  tmp1l = abarli + bbarli/rel(i,k)
! ice
                  tmp1i = abarii + bbarii/rei(i,k)
                  tauxcl(i,k) = cliqwp(i,k)*tmp1l
                  tauxci(i,k) = cicewp(i,k)*tmp1i
               else
                  tauxcl(i,k) = 0.0
                  tauxci(i,k) = 0.0
               endif
! 
! Do not let single scatter albedo be 1.  Delta-eddington solution
! for non-conservative case has different analytic form from solution
! for conservative case, and raddedmx is written for non-conservative case.
! 
               wcl(i,k) = min(tmp2l,.999999_r8)
               gcl(i,k) = ebarli + tmp3l
               fcl(i,k) = gcl(i,k)*gcl(i,k)
! 
               wci(i,k) = min(tmp2i,.999999_r8)
               gci(i,k) = ebarii + tmp3i
               fci(i,k) = gci(i,k)*gci(i,k)
open(unit=99,file='cloudson.txt',status='unknown')
write(99,*) 'clouds on'
close(unit=99)			   
			   
#endif ! TURN_ON_CLOUDS

! 
! Set aerosol properties
! Conversion factor to adjust aerosol extinction (m2/g)
! 

!
! mix dust aerosol size bins
!

               tau_sum = 0
               tau_w_sum = 0
               tau_w_g_sum = 0
               tau_w_f_sum = 0

!EJL - aerosol optical properties have dimensions (spectral bins=96, size bins=40, n_elements=1)
!
#define TURN_ON_AEROSOLS
#ifdef TURN_ON_AEROSOLS
               !EJL - change Na to 40 for the aerosol bins and replace, k, w, g, f, etc with carma values.
               !    - Are the units right on tau_a? I have extinctions of m^2/g?			   
               do Za=1,ncarma_bins
                  ! FAO: prefactor needed to convert m^2/kg --> cm^2/g
                   !tau_a = 10.0 * kaero(ns,Za) * abund_a(i,k,Za)
                   tau_a = 10000.0 * kcarma(1,Za,ns) * abund_a(i,k,Za) ! The prefactor is due to units of 
                              ! g/cm2 for the abundance and m^2/g for k.
                 
                     !if (i .eq. 1 .and. masterproc) then
                     ! grep FAO_N2R camrun.o | uniq | awk '{print $3,$4,$5,$6,$7,$8}'
                     !write(6,'(a,i3,i3,e12.4,e12.4,e12.4)'), 'FAO_AERO', &
                     !     k, Za, wavmid(ns), abund_a(i,k,Za), tau_a
                  !endif

!#define FAO_RESCALE_SSA
!#ifdef FAO_RESCALE_SSA
                  ! FAO: play with single scattering albedo -- experiment 1
                  !if (.true.) then
                     !tau_w = tau_a * (0.90d0*waero(ns,Za))
                     !tau_w = tau_a * (0.90d0*wcarma(ns,Za,1)) !EJL replacing waero with wcarma
                  !else
                  !   ssa_scaled = 0.5d0
                     !tau_w = tau_a * (1.0d0 - ssa_scaled*(1.0d0-waero(ns,Za)))
                     !tau_w = tau_a * (1.0d0 - ssa_scaled*(1.0d0 - wcarma(ns,Za,1)))
                  !endif
!#else ! FAO_RESCALE_SSA
                  !tau_w = tau_a * waero(ns,Za)
                  tau_w = tau_a * wcarma(1,Za,ns)
! EJL - I think this part of the code is not being used. 
!			  open(unit=99,file='wtf.txt',status='unknown')
!              write(99,*) tau_w, tau_a, wcarma(ns,za,1)
!              close(unit=99)
!#endif ! FAO_RESCALE_SSA

                  !tau_w_g = tau_w * gaero(ns,Za)
                  !tau_w_f = tau_w * gaero(ns,Za) * gaero(ns,Za)
                  tau_w_g = tau_w * gcarma(1,Za,ns)
                  tau_w_f = tau_w * gcarma(1,Za,ns) * gcarma(1,Za,ns)

                  tau_sum = tau_sum + tau_a
                  tau_w_sum = tau_w_sum + tau_w
                  tau_w_g_sum = tau_w_g_sum + tau_w_g
                  tau_w_f_sum = tau_w_f_sum + tau_w_f

                  aertau(i,ns,Za) = aertau(i,ns,Za) + tau_a
                  aerssa(i,ns,Za) = aerssa(i,ns,Za) + tau_w
                  aerasm(i,ns,Za) = aerasm(i,ns,Za) + tau_w_g
                  aerfwd(i,ns,Za) = aerfwd(i,ns,Za) + tau_w_f

               enddo !ncarma_bins
#endif  ! TURN_ON_AEROSOLS
!
! mix aerosols
               tau_tot = tau_sum
               tau_w_tot = tau_w_sum
               tau_w_g_tot = tau_w_g_sum
               tau_w_f_tot = tau_w_f_sum

               if (tau_tot .gt. 0.0) then
                 w_tot = tau_w_tot / tau_tot
               else
                 w_tot = 0.0
               endif

               if (tau_w_tot .gt. 0.0) then
                 g_tot = tau_w_g_tot / tau_w_tot
                 f_tot = tau_w_f_tot / tau_w_tot
               else
                 g_tot = 0.0
                 f_tot = 0.0
               endif

               if ( g_tot > 1.0_r8 )  lg_tot_gt1(k)  = .true.
               if ( g_tot < -1.0_r8 ) lg_tot_ltm1(k) = .true.
               if ( f_tot > 1.0_r8 )  lf_tot_gt1(k)  = .true.
               if ( f_tot < 0.0_r8 )  lf_tot_lt0(k)  = .true.

               tauxar(i,k) = tau_tot
               wa(i,k)     = min(w_tot, 0.999999_r8)
               ga(i,k)     = g_tot
               fa(i,k)     = f_tot

         end do  ! i=1,Nday

      end do ! k=1,pver
!EJL
!              open(unit=99,file='tau.txt',status='unknown')
!              write(99,*) tau_tot
!              close(unit=99)

!CSD$ END PARALLEL 

      if (any( lg_tot_gt1(:) )) then
         write(6,*) "g_tot > 1.0"
         call endrun('RADCSWMX')
      end if
      if (any( lg_tot_ltm1(:) )) then
         write(6,*) "g_tot < -1.0"
         call endrun('RADCSWMX')
      end if
      if (any( lf_tot_gt1(:) )) then
         write(6,*) "f_tot > 1.0"
         call endrun('RADCSWMX')
      end if
      if (any( lf_tot_lt0(:) )) then
         write(6,*) "f_tot < 0.0"
         call endrun('RADCSWMX')
      end if

      ! normalize aerosol optical diagnostic fields
      do Za = 1, ncarma_bins
         do i=1,Nday

            if (aerssa(i,ns,Za) .gt. 0.0) then   ! aerssa currently holds product of tau and ssa
               aerasm(i,ns,Za) = aerasm(i,ns,Za) / aerssa(i,ns,Za)
               aerfwd(i,ns,Za) = aerfwd(i,ns,Za) / aerssa(i,ns,Za)
            else
               aerasm(i,ns,Za) = 0.0_r8
               aerfwd(i,ns,Za) = 0.0_r8
            end if
            
            if (aertau(i,ns,Za) .gt. 0.0) then
               aerssa(i,ns,Za) = aerssa(i,ns,Za) / aertau(i,ns,Za)
            else
               aerssa(i,ns,Za) = 0.0_r8
            end if

         end do ! End do i=1,Nday

      end do ! End do Za = 1, Na

! 
! Set reflectivities for long and short wavelength cases
! 
      if (wavmid(ns) < 0.7_r8 ) then
         do i=1,Nday
               albdir(i,ns) = asdir(i)
               albdif(i,ns) = asdif(i)
         end do
      else  ! Wavelength greater than 0.7 micro-meter
         do i=1,Nday
               albdir(i,ns) = aldir(i)
               albdif(i,ns) = aldif(i)
         end do
      end if

! 
! Layer input properties now completely specified; compute the
! delta-Eddington solution reflectivities and transmissivities
! for each layer
! 

      !do i=1,Nday
      !write(*,'(a15,i6,i6,f9.5,f9.5,f9.5,f9.5)') 'FAO_ASDI{R,F}', &
      !       i, ns, asdir(i), aldir(i), asdif(i), aldif(i)
      !enddo

      ! FAO
      do k=0,pver
      do i=1,Nday
         taugab(i,k) = 0.0_r8
         do Zg=1,Ng
            taugab(i,k) = taugab(i,k) + &
                 kgas(ns,Zg) * (pmid_copy(i,k)/1.01325e6)**p_exp(ns,Zg) * abund_g(i,k,Zg)
         enddo

!         if (i .eq. 1 .and. masterproc) then
            ! grep FAO_N2R camrun.o | uniq | awk '{print $3,$4,$5,$6,$7,$8}'
            !write(6,'(a,i3,e12.4,e12.4,e12.4,e12.4,e12.4)') 'FAO_GAS', &
            !     k, wavmid(ns), taugab(i,k),kgas(ns,1),abund_g(i,k,1), gmmr(i,k,1)
!         endif
      enddo
      enddo

!#define DEBUG_TAU
!#ifdef DEBUG_TAU !EJL - uncomment this and endif to stop the writing of these outputs
!      open(unit=7325, file='tau_' // iam_str, access='append')
!      do i=1,1 !Nday
!         do k=0,pver
!            write(7325,'(e16.7e3, e16.7e3, i5, i5, e16.7e3, e16.7e3, e16.7e3, e16.7e3, e16.7e3, e16.7e3)') &
!                 clat_day(i), clon_day(i), k, ns, tau_rayleigh(i,k,ns), taugab(i,k), tauxar(i,k), wa(i,k), ga(i,k), fa(i,k)
!         end do
!      enddo
!
!      close(7325)
!#endif

!open(unit=99, file='tauxar.txt', status='unknown')
!write(99,*) tauxar, 
!close(unit=99)
!          if (ns .eq. 1) then
!          open(unit=99, file='tau1.txt', status='unknown')
!          write(99,'(e16.7e3,e16.7e3,e16.7e3,e16.7e3)') tauxcl
!	  write(99,*) ''
!	  write(99,'(e16.7e3,e16.7e3,e16.7e3,e16.7e3)') tauxci
!	  write(99,*) ''
!	  write(99,'(e16.7e3,e16.7e3,e16.7e3,e16.7e3)') tau_rayleigh
!	  write(99,*) ''
!	  write(99,'(e16.7e3,e16.7e3,e16.7e3,e16.7e3)') tauxar
!	  write(99,*) ''
!	  write(99,'(e16.7e3,e16.7e3,e16.7e3,e16.7e3)') taugab
!	  close(unit=99)
!          endif
!	  
!	  open(unit=99,file='wg.txt',status='unknown')
!	  write(99,'(e16.7e3,e16.7e3,e16.7e3)') wcl, wci, wa
!	  write(99,*) ''!
!	  write(99,'(e16.7e3,e16.7e3,e16.7e3)') gcl, gci, ga
!	  close(unit=99)
!	  endif
!      tauxar(:,:)=0.0

      ! FAO:  to turn off clouds, should have set taux{cl,ci} <-- 0
      call raddedmx(coszrs   ,Nday    , &
              tau_rayleigh(:,:,ns) ,pflx     ,ns       , &
              tauxcl   ,wcl      ,gcl      ,fcl      , &
              tauxci   ,wci      ,gci      ,fci      , &
              tauxar   ,wa       ,ga       ,fa       , &
              rdir_0   ,rdif_0   ,tdir_0   ,tdif_0   ,explay_0  , &
              rdirc_0  ,rdifc_0  ,tdirc_0  ,tdifc_0  ,explayc_0, &
              taugab )

      ! FAO
      do k=0,pver
         do i=1,Nday
            rdir(ns,i,k) = rdir(ns,i,k) + rdir_0(ns,i,k)
            rdif(ns,i,k) = rdif(ns,i,k) + rdif_0(ns,i,k)
            tdir(ns,i,k) = tdir(ns,i,k) + tdir_0(ns,i,k)
            tdif(ns,i,k) = tdif(ns,i,k) + tdif_0(ns,i,k)
            explay(ns,i,k) = explay(ns,i,k) + explay_0(ns,i,k)
            rdirc(ns,i,k) = rdirc(ns,i,k) + rdirc_0(ns,i,k)
            rdifc(ns,i,k) = rdifc(ns,i,k) + rdifc_0(ns,i,k)
            tdirc(ns,i,k) = tdirc(ns,i,k) + tdirc_0(ns,i,k)
            tdifc(ns,i,k) = tdifc(ns,i,k) + tdifc_0(ns,i,k)
            explayc(ns,i,k) = explayc(ns,i,k) + explayc_0(ns,i,k)
         enddo
      enddo

   end do   ! End spectral loop

! 
!----------------------------------------------------------------------
! 
! Solution for max/random cloud overlap.  
! 
! Steps:
! (1. delta-Eddington solution for each layer (called above)
! 
! (2. The adding method is used to
! compute the reflectivity and transmissivity to direct and diffuse
! radiation from the top and bottom of the atmosphere for each
! cloud configuration.  This calculation is based upon the
! max-random overlap assumption.
! 
! (3. to solve for the fluxes, combine the
! bulk properties of the atmosphere above/below the region.
! 
! Index calculations for steps 2-3 are performed outside spectral
! loop to avoid redundant calculations.  Index calculations (with
! application of areamin & nconfgmax conditions) are performed 
! first to identify the minimum subset of terms for the configurations 
! satisfying the areamin & nconfgmax conditions. This minimum set is 
! used to identify the corresponding minimum subset of terms in 
! steps 2 and 3.
! 
   do iconfig = 1, nconfgmax
      ccon(iconfig,0,1:Nday)      = 0
      ccon(iconfig,pverp,1:Nday)  = 0

      icond(iconfig,0,1:Nday)     = iconfig
      iconu(iconfig,pverp,1:Nday) = iconfig
   end do
! 
! Construction of nuniqu/d, istrtu/d, iconu/d using binary tree 
! 
         nuniqd(0,1:Nday) = 1
         nuniqu(pverp,1:Nday) = 1

         istrtd(1,0,1:Nday) = 1
         istrtu(1,pverp,1:Nday) = 1


!CSD$ PARALLEL DO PRIVATE( npasses, kx2, mrgn, region_found, k1, k2, kx1, nxs, ksort, asort ) &
!CSD$ PRIVATE ( ktmp, atmp, cstr, mstr, nstr, cld0, wstr, nrgn, nconfigm, istr, new_term, xwgt ) &
!CSD$ PRIVATE ( j, ptrc, wgtv, km1, nuniq, is0, is1, n0, n1, ptr0, ptr1, kp1, i, irgn ) &
!CSD$ PRIVATE ( k, l, iconfig, l0, isn )
   do i=1,Nday

!----------------------------------------------------------------------
! INDEX CALCULATIONS FOR MAX OVERLAP
! 
! The column is divided into sets of adjacent layers, called regions, 
! in which the clouds are maximally overlapped.  The clouds are
! randomly overlapped between different regions.  The number of
! regions in a column is set by nmxrgn, and the range of pressures
! included in each region is set by pmxrgn.  
! 
! The following calculations determine the number of unique cloud 
! configurations (assuming maximum overlap), called "streams",
! within each region. Each stream consists of a vector of binary
! clouds (either 0 or 100% cloud cover).  Over the depth of the region, 
! each stream requires a separate calculation of radiative properties. These
! properties are generated using the adding method from
! the radiative properties for each layer calculated by raddedmx.
! 
! The upward and downward-propagating streams are treated
! separately.
! 
! We will refer to a particular configuration of binary clouds
! within a single max-overlapped region as a "stream".  We will 
! refer to a particular arrangement of binary clouds over the entire column
! as a "configuration".
! 
! This section of the code generates the following information:
! (1. nrgn    : the true number of max-overlap regions (need not = nmxrgn)
! (2. nstr    : the number of streams in a region (>=1)
! (3. cstr    : flags for presence of clouds at each layer in each stream
! (4. wstr    : the fractional horizontal area of a grid box covered
! by each stream
! (5. kx1,2   : level indices for top/bottom of each region
! 
! The max-overlap calculation proceeds in 3 stages:
! (1. compute layer radiative properties in raddedmx.
! (2. combine these properties between layers 
! (3. combine properties to compute fluxes at each interface.  
! 
! Most of the indexing information calculated here is used in steps 2-3
! after the call to raddedmx.
! 
! Initialize indices for layers to be max-overlapped
! 
! Loop to handle fix in totwgt=0. For original overlap config 
! from npasses = 0.
! 
         npasses = 0
         do
!cdir novector
            do irgn = 0, nmxrgn(i)
               kx2(irgn) = 0
            end do
            mrgn = 0
! 
! Outermost loop over regions (sets of adjacent layers) to be max overlapped
! 
            do irgn = 1, nmxrgn(i)
! 
! Calculate min/max layer indices inside region.  
! 
               region_found = .false.
               if (kx2(irgn-1) < pver) then
                  k1 = kx2(irgn-1)+1
                  kx1(irgn) = k1
                  kx2(irgn) = k1-1
!cdir novector
                  do k2 = pver, k1, -1
                     if (pmid(i,k2) <= pmxrgn(i,irgn)) then
                        kx2(irgn) = k2
                        mrgn = mrgn+1
                        region_found = .true.
                        exit
                     end if
                  end do
               else
                  exit
               endif

               if (region_found) then
! 
! Sort cloud areas and corresponding level indices.  
! 
                  nxs = 0
                  if (cldeps > 0) then 
                     do k = k1,k2
                        if (cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
                           nxs = nxs+1
                           ksort(nxs) = k
! 
! We need indices for clouds in order of largest to smallest, so
! sort 1-cld in ascending order
! 
                           asort(nxs) = 1.0_r8-(floor(cld(i,k)/cldeps)*cldeps)
                        end if
                     end do
                  else
!cdir novector
                     do k = k1,k2
                        if (cld(i,k) >= cldmin) then
                           nxs = nxs+1
                           ksort(nxs) = k
! 
! We need indices for clouds in order of largest to smallest, so
! sort 1-cld in ascending order
! 
                           asort(nxs) = 1.0_r8-cld(i,k)
                        end if
                     end do
                  endif
! 
! If nxs eq 1, no need to sort. 
! If nxs eq 2, sort by swapping if necessary
! If nxs ge 3, sort using local sort routine
! 
                  if (nxs == 2) then
                     if (asort(2) < asort(1)) then
                        ktmp = ksort(1)
                        ksort(1) = ksort(2)
                        ksort(2) = ktmp

                        atmp = asort(1)
                        asort(1) = asort(2)
                        asort(2) = atmp
                     endif
                  else if (nxs >= 3) then
                     call quick_sort(asort(1:nxs),ksort(1:nxs))
                  endif
! 
! Construct wstr, cstr, nstr for this region
! 
!cdir novector
                  cstr(k1:k2,1:nxs+1) = 0
                  mstr = 1
                  cld0 = 0.0_r8
                  do l = 1, nxs
                     if (asort(l) /= cld0) then
                        wstr(mstr,mrgn) = asort(l) - cld0
                        cld0 = asort(l)
                        mstr = mstr + 1
                     endif
!cdir novector
                     cstr(ksort(l),mstr:nxs+1) = 1
                  end do
                  nstr(mrgn) = mstr
                  wstr(mstr,mrgn) = 1.0_r8 - cld0
! 
! End test of region_found = true
! 
               endif
! 
! End loop over regions irgn for max-overlap
! 
            end do
            nrgn = mrgn
! 
! Finish construction of cstr for additional top layer
! 
!cdir novector
            cstr(0,1:nstr(1)) = 0
! 
! INDEX COMPUTATIONS FOR STEP 2-3
! This section of the code generates the following information:
! (1. totwgt     step 3     total frac. area of configurations satisfying
! areamin & nconfgmax criteria
! (2. wgtv       step 3     frac. area of configurations 
! (3. ccon       step 2     binary flag for clouds in each configuration
! (4. nconfig    steps 2-3  number of configurations
! (5. nuniqu/d   step 2     Number of unique cloud configurations for
! up/downwelling rad. between surface/TOA
! and level k
! (6. istrtu/d   step 2     Indices into iconu/d
! (7. iconu/d    step 2     Cloud configurations which are identical
! for up/downwelling rad. between surface/TOA
! and level k
! 
! Number of configurations (all permutations of streams in each region)
! 
            nconfigm = product(nstr(1: nrgn))
! 
! Construction of totwgt, wgtv, ccon, nconfig
! 
!cdir novector
            istr(1: nrgn) = 1
            nconfig(i) = 0
            totwgt(i) = 0.0_r8
            new_term = .true.
            do iconfig = 1, nconfigm
               xwgt = 1.0_r8
!cdir novector
               do mrgn = 1,  nrgn
                  xwgt = xwgt * wstr(istr(mrgn),mrgn)
               end do
               if (xwgt >= areamin) then
                  nconfig(i) = nconfig(i) + 1
                  if (nconfig(i) <= nconfgmax) then
                     j = nconfig(i)
                     ptrc(nconfig(i)) = nconfig(i)
                  else
                     nconfig(i) = nconfgmax
                     if (new_term) then
                        min_idx = minloc(wgtv)
                        j = min_idx(1)
                     endif
                     if (wgtv(j) < xwgt) then
                        totwgt(i) = totwgt(i) - wgtv(j)
                        new_term = .true.
                     else
                        new_term = .false.
                     endif
                  endif
                  if (new_term) then
                     wgtv(j) = xwgt
                     totwgt(i) = totwgt(i) + xwgt
!cdir novector
                     do mrgn = 1, nrgn
!cdir novector
                        ccon(j,kx1(mrgn):kx2(mrgn),i) = cstr(kx1(mrgn):kx2(mrgn),istr(mrgn))
                     end do
                  endif
               endif

               mrgn =  nrgn
               istr(mrgn) = istr(mrgn) + 1
               do while (istr(mrgn) > nstr(mrgn) .and. mrgn > 1)
                  istr(mrgn) = 1
                  mrgn = mrgn - 1
                  istr(mrgn) = istr(mrgn) + 1
               end do
! 
! End do iconfig = 1, nconfigm
! 
            end do
! 
! If totwgt(i) = 0 implement maximum overlap and make another pass
! if totwgt(i) = 0 on this second pass then terminate.
! 
            if (totwgt(i) > 0.) then
               exit
            else
               npasses = npasses + 1
               if (npasses >= 2 ) then
                  write(6,*)'RADCSWMX: Maximum overlap of column ','failed'
                  call endrun('RADCSWMX')
               endif
               nmxrgn(i)=1
               pmxrgn(i,1)=1.0e30
            end if
!
! End npasses = 0, do
!
         end do
! 
! Finish construction of ccon
! 

         istrtd(2,0,i) = nconfig(i)+1
         istrtu(2,pverp,i) = nconfig(i)+1

         do k = 1, pverp
            km1 = k-1
            nuniq = 0
            istrtd(1,k,i) = 1
!cdir novector
            do l0 = 1, nuniqd(km1,i)
               is0 = istrtd(l0,km1,i)
               is1 = istrtd(l0+1,km1,i)-1
               n0 = 0
               n1 = 0
!cdir novector
               do isn = is0, is1
                  j = icond(isn,km1,i)
                  if (ccon(j,k,i) == 0) then
                     n0 = n0 + 1
                     ptr0(n0) = j
                  else       ! if (ccon(j,k,i) == 1) then
                     n1 = n1 + 1
                     ptr1(n1) = j
                  endif
               end do
               if (n0 > 0) then
                  nuniq = nuniq + 1
                  istrtd(nuniq+1,k,i) = istrtd(nuniq,k,i)+n0
!cdir novector
                  icond(istrtd(nuniq,k,i):istrtd(nuniq+1,k,i)-1,k,i) =  ptr0(1:n0)
               endif
               if (n1 > 0) then
                  nuniq = nuniq + 1
                  istrtd(nuniq+1,k,i) = istrtd(nuniq,k,i)+n1
!cdir novector
                  icond(istrtd(nuniq,k,i):istrtd(nuniq+1,k,i)-1,k,i) =  ptr1(1:n1)
               endif
            end do
            nuniqd(k,i) = nuniq
         end do
!
!  Find 'transition point' in downward configurations where the number
!  of 'configurations' changes from 1.  This is used to optimize the
!  construction of the upward configurations.
!  Note: k1 == transition point
!

         do k = pverp,0,-1
           if ( nuniqd(k,i) == 1) then
              k1 = k
              exit
           end if
         end do

         do k = pver, k1+1,-1
            kp1 = k+1
            nuniq = 0
            istrtu(1,k,i) = 1
!cdir novector
            do l0 = 1, nuniqu(kp1,i)
               is0 = istrtu(l0,kp1,i)
               is1 = istrtu(l0+1,kp1,i)-1
               n0 = 0
               n1 = 0
!cdir novector
               do isn = is0, is1
                  j = iconu(isn,kp1,i)
                  if (ccon(j,k,i) == 0) then
                     n0 = n0 + 1
                     ptr0(n0) = j
                  else       ! if (ccon(j,k,i) == 1) then
                     n1 = n1 + 1
                     ptr1(n1) = j
                  endif
               end do
               if (n0 > 0) then
                  nuniq = nuniq + 1
                  istrtu(nuniq+1,k,i) = istrtu(nuniq,k,i)+n0
!cdir novector
                  iconu(istrtu(nuniq,k,i):istrtu(nuniq+1,k,i)-1,k,i) =  ptr0(1:n0)
               endif
               if (n1 > 0) then
                  nuniq = nuniq + 1
                  istrtu(nuniq+1,k,i) = istrtu(nuniq,k,i)+n1
!cdir novector
                  iconu(istrtu(nuniq,k,i):istrtu(nuniq+1,k,i)-1,k,i) = ptr1(1:n1)
               endif
            end do
            nuniqu(k,i) = nuniq
         end do
!
!  Copy identical configurations from 'transition point' to surface.
!
         k1 = min(pverp-1,k1)
         nuniq = nuniqu(k1+1,i)
         do k = k1,0,-1
            nuniqu(k,i) = nuniq
!cdir novector
            iconu(1:nuniq,k,i) = iconu(1:nuniq,k1+1,i)
!cdir novector
            istrtu(1:nuniq+1,k,i) = istrtu(1:nuniq+1,k1+1,i)
         end do

!cdir novector
         v_wgtv(1:nconfig(i),i) = wgtv(1:nconfig(i))

! 
!----------------------------------------------------------------------
! End of index calculations
!----------------------------------------------------------------------
! 
! End do i=1,Nday
! 
   end do
!CSD$ END PARALLEL 

!----------------------------------------------------------------------
! Start of flux calculations
!----------------------------------------------------------------------
!
! Initialize spectrally integrated totals:
! 
         totfld(1:Nday,0:pver) = 0.0_r8
         fswup (1:Nday,0:pver) = 0.0_r8
         fswdn (1:Nday,0:pver) = 0.0_r8

         sfltot(1:Nday)        = 0.0_r8
         fswup (1:Nday,pverp)  = 0.0_r8
         fswdn (1:Nday,pverp)  = 0.0_r8
! 
! Start spectral interval
! 
!old   do ns = 1,Nf
!old     wgtint = nirwgt(ns)

     do i=1,Nday

!----------------------------------------------------------------------
! STEP 2
! 
! 
! Apply adding method to solve for radiative properties
! 
! first initialize the bulk properties at toa
! 

! Nf, 0:pverp, nconfgmax, pcols

            rdndif(:,0,1:nconfig(i),i) = 0.0_r8
            exptdn(:,0,1:nconfig(i),i) = 1.0_r8
            tdntot(:,0,1:nconfig(i),i) = 1.0_r8
! 
! End do i=1,Nday
! 
     end do
! 
! solve for properties involving downward propagation of radiation.
! the bulk properties are:
! 
! (1. exptdn   sol. beam dwn. trans from layers above
! (2. rdndif   ref to dif rad for layers above
! (3. tdntot   total trans for layers above
! 

!CSD$ PARALLEL DO PRIVATE( km1, is0, is1, j, jj, Ttdif, Trdif, Trdir, Ttdir, Texplay ) &
!CSD$ PRIVATE( xexpt, xrdnd, tdnmexp,  ytdnd, yrdnd, rdenom, rdirexp, zexpt, zrdnd, ztdnt ) &
!CSD$ PRIVATE( i, k, l0, ns, isn )
         do i = 1, Nday
            do k = 1, pverp
               km1 = k - 1
!cdir nodep
               do l0 = 1, nuniqd(km1,i)
                  is0 = istrtd(l0,km1,i)
                  is1 = istrtd(l0+1,km1,i)-1

                  j = icond(is0,km1,i)

! 
! If cloud in layer, use cloudy layer radiative properties (ccon == 1)
! If clear layer, use clear-sky layer radiative properties (ccon /= 1)
! 
                  if ( ccon(j,km1,i) == 1 ) then
                     Ttdif(:) = tdif(:,i,km1)
                     Trdif(:) = rdif(:,i,km1)
                     Trdir(:) = rdir(:,i,km1)
                     Ttdir(:) = tdir(:,i,km1)
                     Texplay(:) = explay(:,i,km1)
                  else
                     Ttdif(:) = tdifc(:,i,km1)
                     Trdif(:) = rdifc(:,i,km1)
                     Trdir(:) = rdirc(:,i,km1)
                     Ttdir(:) = tdirc(:,i,km1)
                     Texplay(:) = explayc(:,i,km1)
                  end if

                  do ns = 1, Nf
                  xexpt   = exptdn(ns,km1,j,i)
                  xrdnd   = rdndif(ns,km1,j,i)
                  tdnmexp = tdntot(ns,km1,j,i) - xexpt

                  ytdnd = Ttdif(ns)
                  yrdnd = Trdif(ns)

                  rdenom  = 1._r8/(1._r8-yrdnd*xrdnd)
                  rdirexp = Trdir(ns)*xexpt

                  zexpt = xexpt * Texplay(ns)
                  zrdnd = yrdnd + xrdnd*(ytdnd**2)*rdenom
                  ztdnt = xexpt*Ttdir(ns) + ytdnd* &
                          (tdnmexp + xrdnd*rdirexp)*rdenom

                  exptdn(ns,k,j,i) = zexpt
                  rdndif(ns,k,j,i) = zrdnd
                  tdntot(ns,k,j,i) = ztdnt
                  end do ! ns = 1, Nf
!
! If 2 or more configurations share identical properties at a given level k,
! the properties (at level k) are computed once and copied to
! all the configurations for efficiency.
!
                  do isn = is0+1, is1
                     jj = icond(isn,km1,i)
                     exptdn(:,k,jj,i) = exptdn(:,k,j,i)
                     rdndif(:,k,jj,i) = rdndif(:,k,j,i)
                     tdntot(:,k,jj,i) = tdntot(:,k,j,i)
                  end do

! 
! end do l0 = 1, nuniqd(k,i)
! 
               end do
! 
! end do k = 1, pverp
! 
            end do
! 
! end do i = 1, Nday
! 
         end do
!CSD$ END PARALLEL 
! 
! Solve for properties involving upward propagation of radiation.
! The bulk properties are:
! 
! (1. rupdif   Ref to dif rad for layers below
! (2. rupdir   Ref to dir rad for layers below
! 
! Specify surface boundary conditions (surface albedos)
! 


! Nf, 0:pverp, nconfgmax, pcols
   rupdir = 0._r8
   rupdif = 0._r8
   do i = 1, Nday
      do ns = 1, Nf
         rupdir(ns,pverp,1:nconfig(i),i) = albdir(i,ns)
         rupdif(ns,pverp,1:nconfig(i),i) = albdif(i,ns)
      end do
   end do

         do i = 1, Nday
            do k = pver, 0, -1
               do l0 = 1, nuniqu(k,i)
                  is0 = istrtu(l0,k,i)
                  is1 = istrtu(l0+1,k,i)-1

                  j = iconu(is0,k,i)

! 
! If cloud in layer, use cloudy layer radiative properties (ccon == 1)
! If clear layer, use clear-sky layer radiative properties (ccon /= 1)
! 
                  if ( ccon(j,k,i) == 1 ) then
                     Ttdif(:) = tdif(:,i,k)
                     Trdif(:) = rdif(:,i,k)
                     Trdir(:) = rdir(:,i,k)
                     Ttdir(:) = tdir(:,i,k)
                     Texplay(:) = explay(:,i,k)
                  else
                     Ttdif(:) = tdifc(:,i,k)
                     Trdif(:) = rdifc(:,i,k)
                     Trdir(:) = rdirc(:,i,k)
                     Ttdir(:) = tdirc(:,i,k)
                     Texplay(:) = explayc(:,i,k)
                  end if

                  do ns = 1, Nf
                  xrupd = rupdif(ns,k+1,j,i)
                  xrups = rupdir(ns,k+1,j,i)

! 
! If cloud in layer, use cloudy layer radiative properties (ccon == 1)
! If clear layer, use clear-sky layer radiative properties (ccon /= 1)
! 
                  yexpt = Texplay(ns)
                  yrupd = Trdif(ns)
                  ytupd = Ttdif(ns)

                  rdenom  = 1._r8/( 1._r8 - yrupd*xrupd)
                  tdnmexp = (Ttdir(ns)-yexpt)
                  rdirexp = xrups*yexpt

                  zrupd = yrupd + xrupd*(ytupd**2)*rdenom
                  zrups = Trdir(ns) + ytupd*(rdirexp + xrupd*tdnmexp)*rdenom

                  rupdif(ns,k,j,i) = zrupd
                  rupdir(ns,k,j,i) = zrups
                  end do ! ns = 1, Nf
!
! If 2 or more configurations share identical properties at a given level k,
! the properties (at level k) are computed once and copied to
! all the configurations for efficiency.
!
                  do isn = is0+1, is1
                     jj = iconu(isn,k,i)
                     rupdif(:,k,jj,i) = rupdif(:,k,j,i)
                     rupdir(:,k,jj,i) = rupdir(:,k,j,i)
                  end do

! 
! end do l0 = 1, nuniqu(k,i)
! 
               end do
! 
! end do k = pver,0,-1
! 
            end do
! 
! end do i = 1, Nday
! 
         end do

! 
!----------------------------------------------------------------------
! 
! STEP 3
! 
! Compute up and down fluxes for each interface k.  This requires
! adding up the contributions from all possible permutations
! of streams in all max-overlap regions, weighted by the
! product of the fractional areas of the streams in each region
! (the random overlap assumption).  The adding principle has been
! used in step 2 to combine the bulk radiative properties 
! above and below the interface.
! 

! 
! Initialize the fluxes
! 
            fluxup = 0.0_r8
            fluxdn = 0.0_r8

            do i = 1, Nday
!cdir novector
            do iconfig = 1, nconfig(i)
               xwgt = v_wgtv(iconfig,i)

!cdir collapse
               do k = 0, pverp
                  do ns = 1, Nf
                  xexpt = exptdn(ns,k,iconfig,i)
                  xtdnt = tdntot(ns,k,iconfig,i)
                  xrdnd = rdndif(ns,k,iconfig,i)
                  xrupd = rupdif(ns,k,iconfig,i)
                  xrups = rupdir(ns,k,iconfig,i)
! 
! Flux computation
! 
                  rdenom = 1._r8/(1._r8 - xrdnd * xrupd)

                  fluxup(ns,k,i) = fluxup(ns,k,i) + xwgt *  &
                              ((xexpt * xrups + (xtdnt - xexpt) * xrupd) * rdenom)
                  fluxdn(ns,k,i) = fluxdn(ns,k,i) + xwgt *  &
                              (xexpt + (xtdnt - xexpt + xexpt * xrups * xrdnd) * rdenom)

                  !if (k > 60 .and. ns > 92 .and. clat_day(i) < -.52 .and. clat_day(i) > -.53 .and. &
                  !    clon_day(i) > 3.66 .and. clon_day(i) < 3.67)  &
                  !write(*,'(A5,I4,I4,I4,F9.6,F9.6)') 'CCC', i, k, ns, rdenom, fluxup(ns,k,i)
                  end do ! do ns = 1, Nf
               end do
! 
! End do iconfig = 1, nconfig(i)
! 
            end do
! 
! End do iconfig = 1, Nday
! 
            end do


! 
! Normalize by total area covered by cloud configurations included
! in solution
! 
#ifdef JPE_VMATH
            call vrec(v_rtotwgt,totwgt,Nday)
#endif
            do i = 1, Nday
            do k = 0, pverp
            do ns = 1, Nf
#ifdef JPE_VMATH
               fluxup(ns,k,i)=fluxup(ns,k,i) * v_rtotwgt(i)
               fluxdn(ns,k,i)=fluxdn(ns,k,i) * v_rtotwgt(i)
#else
               fluxup(ns,k,i)=fluxup(ns,k,i) / totwgt(i)
               fluxdn(ns,k,i)=fluxdn(ns,k,i) / totwgt(i)
#endif
            end do ! do i = 1, nday
            end do ! do k = 0, pverp
            end do ! do i = 1, nday

!zero albedo arrays
alb325(:)=0.0_r8
alb375(:)=0.0_r8
alb425(:)=0.0_r8
alb475(:)=0.0_r8
alb525(:)=0.0_r8
alb575(:)=0.0_r8
alb642(:)=0.0_r8
alb714(:)=0.0_r8
alb784(:)=0.0_r8
alb845(:)=0.0_r8
alb891(:)=0.0_r8
alb940(:)=0.0_r8
alb985(:)=0.0_r8
alb1070(:)=0.0_r8
alb1140(:)=0.0_r8
alb1220(:)=0.0_r8
alb1290(:)=0.0_r8
alb1380(:)=0.0_r8
alb1490(:)=0.0_r8
alb1610(:)=0.0_r8
alb1750(:)=0.0_r8
alb1910(:)=0.0_r8
alb2110(:)=0.0_r8
alb2350(:)=0.0_r8

!
! 
! Initialize the direct-beam flux at surface
! 
            wexptdn(:,1:Nday) = 0.0_r8

   do ns = 1,Nf
      wgtint = nirwgt(ns)


      do i=1,Nday
      do iconfig = 1, nconfig(i)
!
! Note: exptdn can be directly indexed by iconfig at k=pverp.
!
         wexptdn(ns,i) =  wexptdn(ns,i) + v_wgtv(iconfig,i) * exptdn(ns,pverp,iconfig,i)
      end do
      end do

      do i=1,Nday
#ifdef JPE_VMATH
         wexptdn(ns,i) = wexptdn(ns,i) * v_rtotwgt(i)
#else
         wexptdn(ns,i) = wexptdn(ns,i) / totwgt(i)
#endif
! 
! Monochromatic computation completed; accumulate in totals
! 
         ! FAO:  kludge -- one gas only (CH4) at this time!!!
         !solflx(i)   = solin(i)*frcsol(ns)*psf(ns)
         solflx(i)   = solin(i)*frcsol(ns)*spectralweight(ns,1)
		 
         
		 fsnt(i)  = fsnt(i) + solflx(i)*(fluxdn(ns,1,i) - fluxup(ns,1,i))
         fsntoa(i)= fsntoa(i) + solflx(i)*(fluxdn(ns,0,i) - fluxup(ns,0,i))
         fsns(i)  = fsns(i) + solflx(i)*(fluxdn(ns,pverp,i)-fluxup(ns,pverp,i))
         sfltot(i)   = sfltot(i) + solflx(i)
         fswup(i,0) = fswup(i,0) + solflx(i)*fluxup(ns,0,i)
         fswdn(i,0) = fswdn(i,0) + solflx(i)*fluxdn(ns,0,i)
! 
! Down spectral fluxes need to be in mks; thus the .001 conversion factors
! 
         if (wavmid(ns) < 0.7_r8) then
            sols(i)  = sols(i) + wexptdn(ns,i)*solflx(i)*0.001_r8
            solsd(i) = solsd(i)+(fluxdn(ns,pverp,i)-wexptdn(ns,i))*solflx(i)*0.001_r8 
         else
            soll(i)  = soll(i) + wexptdn(ns,i)*solflx(i)*0.001_r8
            solld(i) = solld(i)+(fluxdn(ns,pverp,i)-wexptdn(ns,i))*solflx(i)*0.001_r8 
            fsnrtoaq(i) = fsnrtoaq(i) + solflx(i)*(fluxdn(ns,0,i) - fluxup(ns,0,i))
         end if
         fsnirtoa(i) = fsnirtoa(i) + wgtint*solflx(i)*(fluxdn(ns,0,i) - fluxup(ns,0,i))

! 
! End do i=1,Nday
! 
      end do

      do k=0,pver
      do i=1,Nday
! 
! Compute flux divergence in each layer using the interface up and down
! fluxes:
! 
         kp1 = k+1

         flxdiv = (fluxdn(ns,k,i) - fluxdn(ns,kp1,i)) + (fluxup(ns,kp1,i) - fluxup(ns,k,i))
         totfld(i,k)  = totfld(i,k)  + solflx(i)*flxdiv
         fswdn(i,kp1) = fswdn(i,kp1) + solflx(i)*fluxdn(ns,kp1,i)
         fswup(i,kp1) = fswup(i,kp1) + solflx(i)*fluxup(ns,kp1,i)
         fns(i,kp1)   = fswdn(i,kp1) - fswup(i,kp1)
      end do
      end do
! 
! Perform clear-sky calculation
! 

      exptdnc(1:Nday,0) =   1.0_r8
      rdndifc(1:Nday,0) =   0.0_r8
      tdntotc(1:Nday,0) =   1.0_r8
      rupdirc(1:Nday,pverp) = albdir(1:Nday,ns)
      rupdifc(1:Nday,pverp) = albdif(1:Nday,ns)

!cdir expand=pverp
      do k = 1, pverp
      do i=1,Nday
         km1 = k - 1
         xexpt = exptdnc(i,km1)
         xrdnd = rdndifc(i,km1)
         yrdnd = rdifc(ns,i,km1)
         ytdnd = tdifc(ns,i,km1)

         exptdnc(i,k) = xexpt*explayc(ns,i,km1)

         rdenom  = 1._r8/(1._r8 - yrdnd*xrdnd)
         rdirexp = rdirc(ns,i,km1)*xexpt
         tdnmexp = tdntotc(i,km1) - xexpt

         tdntotc(i,k) = xexpt*tdirc(ns,i,km1) + ytdnd*(tdnmexp + xrdnd*rdirexp)* &
              rdenom
         rdndifc(i,k) = yrdnd + xrdnd*(ytdnd**2)*rdenom
      end do
      end do

      do k=pver,0,-1
      do i=1,Nday
         xrupd = rupdifc(i,k+1)
         yexpt = explayc(ns,i,k)
         yrupd = rdifc(ns,i,k)
         ytupd = tdifc(ns,i,k)

         rdenom = 1._r8/( 1._r8 - yrupd*xrupd)

         rupdirc(i,k) = rdirc(ns,i,k) + ytupd*(rupdirc(i,k+1)*yexpt + &
              xrupd*(tdirc(ns,i,k)-yexpt))*rdenom
         rupdifc(i,k) = yrupd + xrupd*ytupd**2*rdenom
      end do
      end do

      do k=0,pverp
      do i=1,Nday
         rdenom    = 1._r8/(1._r8 - rdndifc(i,k)*rupdifc(i,k))
         fluxup(ns,k,i) = (exptdnc(i,k)*rupdirc(i,k) + (tdntotc(i,k)-exptdnc(i,k))*rupdifc(i,k))* &
              rdenom
         fluxdn(ns,k,i) = exptdnc(i,k) + &
              (tdntotc(i,k) - exptdnc(i,k) + exptdnc(i,k)*rupdirc(i,k)*rdndifc(i,k))* &
              rdenom
			  
		 !EJL - albedo = fluxup/(fluxdn) at top of atmosphere
		 ejlalb(i,ns) = fluxup(ns,0,i) / fluxdn(ns,0,i)	 ! Did I check to see if 0 was toa? 
      end do
      end do

      do i=1,Nday
         fsntc(i)    = fsntc(i)+solflx(i)*(fluxdn(ns,1,i)-fluxup(ns,1,i))
         fsntoac(i)  = fsntoac(i)+solflx(i)*(fluxdn(ns,0,i)-fluxup(ns,0,i))
         fsnsc(i)    = fsnsc(i)+solflx(i)*(fluxdn(ns,pverp,i)-fluxup(ns,pverp,i))
         fsdsc(i)    = fsdsc(i)+solflx(i)*(fluxdn(ns,pverp,i))
         fsnrtoac(i) = fsnrtoac(i)+wgtint*solflx(i)*(fluxdn(ns,0,i)-fluxup(ns,0,i))
      end do


      do k = 1,pverp
      do i=1,Nday
         fcns(i,k)=fcns(i,k) + solflx(i)*(fluxdn(ns,k,i)-fluxup(ns,k,i))
      enddo
      enddo
! 
! End of clear sky calculation
! 
! 
! End of spectral interval loop
!    
   end do
   
   
!open(unit=99, file='flux2.txt', status='unknown')
!write(99,*) Nf, Nday,fluxup(:,0,1)
!write(99,*) ''
!write(99,*) fluxup(1,0,:)
!close(unit=99)   
   

   
!EJL - calculating the albedo. sum of (fluxup*spectralweight*frcsol) for each wavelength 
 !     divided by the sum of the flux down...
   do i=1,Nday
     do aaa=1,24
     ejlflxup(i,aaa) = frcsol(4*(aaa-1)+1)*spectralweight(4*(aaa-1)+1,1)*fluxup(4*(aaa-1)+1,0,i) &
	 +frcsol(4*(aaa-1)+2)*spectralweight(4*(aaa-1)+2,1)*fluxup(4*(aaa-1)+2,0,i) &
	 +frcsol(4*(aaa-1)+3)*spectralweight(4*(aaa-1)+3,1)*fluxup(4*(aaa-1)+3,0,i) &
	 +frcsol(4*(aaa-1)+4)*spectralweight(4*(aaa-1)+4,1)*fluxup(4*(aaa-1)+4,0,i)
	 
     ejlflxdn(i,aaa) = frcsol(4*(aaa-1)+1)*spectralweight(4*(aaa-1)+1,1)*fluxdn(4*(aaa-1)+1,0,i) &
	 +frcsol(4*(aaa-1)+2)*spectralweight(4*(aaa-1)+2,1)*fluxdn(4*(aaa-1)+2,0,i) &
	 +frcsol(4*(aaa-1)+3)*spectralweight(4*(aaa-1)+3,1)*fluxdn(4*(aaa-1)+3,0,i) &
	 +frcsol(4*(aaa-1)+4)*spectralweight(4*(aaa-1)+4,1)*fluxdn(4*(aaa-1)+4,0,i)
     enddo

       alb325(i)=ejlflxup(i,1)/ejlflxdn(i,1)
	   alb375(i)=ejlflxup(i,2)/ejlflxdn(i,2)
	   alb425(i)=ejlflxup(i,3)/ejlflxdn(i,3)
	   alb475(i)=ejlflxup(i,4)/ejlflxdn(i,4)
	   alb525(i)=ejlflxup(i,5)/ejlflxdn(i,5)
	   alb575(i)=ejlflxup(i,6)/ejlflxdn(i,6)
	   alb642(i)=ejlflxup(i,7)/ejlflxdn(i,7)
	   alb714(i)=ejlflxup(i,8)/ejlflxdn(i,8)
	   alb784(i)=ejlflxup(i,9)/ejlflxdn(i,9)
	   alb845(i)=ejlflxup(i,10)/ejlflxdn(i,10)
	   alb891(i)=ejlflxup(i,11)/ejlflxdn(i,11)
	   alb940(i)=ejlflxup(i,12)/ejlflxdn(i,12)
	   alb985(i)=ejlflxup(i,13)/ejlflxdn(i,13)
	   alb1070(i)=ejlflxup(i,14)/ejlflxdn(i,14)
	   alb1140(i)=ejlflxup(i,15)/ejlflxdn(i,15)
	   alb1220(i)=ejlflxup(i,16)/ejlflxdn(i,16)
	   alb1290(i)=ejlflxup(i,17)/ejlflxdn(i,17)
	   alb1380(i)=ejlflxup(i,18)/ejlflxdn(i,18)
	   alb1490(i)=ejlflxup(i,19)/ejlflxdn(i,19)
	   alb1610(i)=ejlflxup(i,20)/ejlflxdn(i,20)
	   alb1750(i)=ejlflxup(i,21)/ejlflxdn(i,21)
	   alb1910(i)=ejlflxup(i,22)/ejlflxdn(i,22)
	   alb2110(i)=ejlflxup(i,23)/ejlflxdn(i,23)
	   alb2350(i)=ejlflxup(i,24)/ejlflxdn(i,24)

   end do 
!open(unit=99,file='alb.txt',status='unknown')
!write(99,*) alb375
!close(unit=99)

   do i=1,Nday

! 
! Compute solar heating rate (J/kg/s)
! 
!cdir expand=pver
      do k=1,pver
         qrs(i,k) = -1.E-4*gravit*totfld(i,k)/(pint(i,k) - pint(i,k+1))
      end do

! 
! Set the downwelling flux at the surface 
! 
      fsds(i) = fswdn(i,pverp)

   end do  ! End do i=1,Nday
   
!open(unit=99,file='qrs.txt',status='unknown') !EJL
!write(99,*) qrs
!close(unit=99)

!   FAO: debug short wave output
!   integer, intent(inout) ::  nmxrgn(pcols)    ! Number of maximally overlapped regions
!   real(r8), intent(out) :: solin(pcols)     ! Incident solar flux
!   real(r8), intent(out) :: fsns(pcols)      ! Surface absorbed solar flux
!   real(r8), intent(out) :: fsnt(pcols)      ! Total column absorbed solar flux
!   real(r8), intent(out) :: fsntoa(pcols)    ! Net solar flux at TOA
!   real(r8), intent(out) :: fsds(pcols)      ! Flux shortwave downwelling surface
!   real(r8), intent(out) :: fsnsc(pcols)     ! Clear sky surface absorbed solar flux
!   real(r8), intent(out) :: fsdsc(pcols)     ! Clear sky surface downwelling solar flux
!   real(r8), intent(out) :: fsntc(pcols)     ! Clear sky total column absorbed solar flx
!   real(r8), intent(out) :: fsntoac(pcols)   ! Clear sky net solar flx at TOA
!   real(r8), intent(out) :: sols(pcols)      ! Direct solar rad on surface (< 0.7)
!   real(r8), intent(out) :: soll(pcols)      ! Direct solar rad on surface (>= 0.7)
!   real(r8), intent(out) :: solsd(pcols)     ! Diffuse solar rad on surface (< 0.7)
!   real(r8), intent(out) :: solld(pcols)     ! Diffuse solar rad on surface (>= 0.7)
!   real(r8), intent(out) :: fsnirtoa(pcols)  ! Near-IR flux absorbed at toa
!   real(r8), intent(out) :: fsnrtoac(pcols)  ! Clear sky near-IR flux absorbed at toa
!   real(r8), intent(out) :: fsnrtoaq(pcols)  ! Net near-IR flux at toa >= 0.7 microns
!   real(r8), intent(out) :: frc_day(pcols) ! = 1 for daylight, =0 for night columns

!   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!                                                 !    maximally overlapped region. 
!                                                 !    0->pmxrgn(i,1) is range of pressure for
!                                                 !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                                 !    2nd region, etc
!   real(r8), intent(out) :: qrs(pcols,pver)  ! Solar heating rate
!   real(r8), intent(out) :: fns(pcols,pverp)   ! net flux at interfaces
!   real(r8), intent(out) :: fcns(pcols,pverp)  ! net clear-sky flux at interfaces

!   real(r8) :: aertau(pcols,Nf,Na) ! Aerosol column optical depth
!   real(r8) :: aerssa(pcols,Nf,Na) ! Aerosol column averaged single scattering albedo
!   real(r8) :: aerasm(pcols,Nf,Na) ! Aerosol column averaged asymmetry parameter
!   real(r8) :: aerfwd(pcols,Nf,Na) ! Aerosol column averaged forward scattering


#ifdef EXAMPLE_MPI
   snd = 0
   rcv = 0
   do i=1,Nday
      ! if dlats are all the same ...
      ! area = dphi/(2*pi) * (cos(clat(i) - dlat/2)  - cos(clat(i) + dlat/2))
      ! area = dphi/(2*pi) * (2*sin(lat)*sin(dlat/2))
      do k=1,pverp
         !snd(Ilat(i),k) = snd(Ilat(i),k) + fswdn(i,k) * area
      enddo
      ! tack on area to snd buffer
      !snd(Ilat(i),pverp+1) = snd(Ilat(i),pverp+1) + area
   enddo
   call mpi_allreduce(snd, rcv, plat*(pverp+1), MPI_DOUBLE_PRECISION, MPI_SUM, mpicom, ierr)
   x_sum = x_sum + rcv
#endif

!#define VERBOSE_DEBUG
#ifdef VERBOSE_DEBUG
   open(unit=7321, file='solt_' // iam_str, access='append')
   open(unit=7322, file='solz_' // iam_str, access='append')
   open(unit=7323, file='solznu_' // iam_str, access='append')
   open(unit=7324, file='sol_' // iam_str, access='append')

   do i=1,Nday
      write(7321,'(f8.2, f8.2, e16.7e3, e16.7e3, e16.7e3, e16.7e3, e16.7e3, e16.7e3, e16.7e3)') &
           clat_day(i)*180/3.14159, clon_day(i)*180/3.14159, coszrs(i), solin(i), solflx(i), &
           fsns(i), fsnt(i), fsntoa(i), fsds(i)
           !sols(i),soll(i),solsd(i),solld(i),fsnirtoa(i),fsnrtoaq(i),frc_day(i)
   enddo
   do i=1,Nday
      do k=0,pverp
         write(7322,'(f8.2,f8.2, i5, e16.7e3, e16.7e3)') &
              clat_day(i)*180/3.14159, clon_day(i)*180/3.14159, k, fswdn(i,k), fswup(i,k)
      end do
   enddo

   do i=1,Nday
   do ns=1,Nf
   do k=0,pverp
      if (k==pverp) then
         kp1 = k
      else
         kp1 = k+1
      endif

      write(7323,'(f8.2,f8.2, i5, i5, e16.7e3, e16.7e3, e16.7e3)') &
           clat_day(i)*180/3.14159, clon_day(i)*180/3.14159, ns, k, fluxup(ns,k,i), fluxdn(ns,k,i), &
           (fluxdn(ns,k,i) - fluxdn(ns,kp1,i)) + (fluxup(ns,kp1,i) - fluxup(ns,k,i))
   enddo
   enddo
   enddo

   do i=1,Nday
   do ns=1,Nf
      write(7324,'(f8.2,f8.2, i5, e16.7e3, e16.7e3, e16.7e3)') &
           clat_day(i)*180/3.14159, clon_day(i)*180/3.14159, ns, solin(i), frcsol(ns), spectralweight(ns,1)
   enddo
   enddo

   close(7324)
   close(7323)
   close(7322)
   close(7321)
#endif

!
! Rearrange output arrays.
!
! intent(inout)
!
   call ExpDayNite(pmxrgn,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call ExpDayNite(nmxrgn,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
!
! intent(out)
!
   call ExpDayNite(solin,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(qrs,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)

   call ExpDayNite(ejlalb,  Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, Nf)
   call ExpDayNite(alb325,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb375,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb425,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb475,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb525,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb575,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb642,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb714,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb784,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb845,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb891,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb940,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb985,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1070,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1140,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1220,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1290,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1380,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1490,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1610,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1750,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb1910,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb2110,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(alb2350,  Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   
   call ExpDayNite(fns,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call ExpDayNite(fcns,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call ExpDayNite(fsns,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnt,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsds,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnsc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsdsc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntoac,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(sols,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(soll,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(solsd,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(solld,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnirtoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnrtoac,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnrtoaq,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(frc_day,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)

   call ExpDayNite(aertau,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, Nf, 1, ncarma_bins)
   call ExpDayNite(aerssa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, Nf, 1, ncarma_bins)
   call ExpDayNite(aerasm,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, Nf, 1, ncarma_bins)
   call ExpDayNite(aerfwd,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, Nf, 1, ncarma_bins)


! EJL - ouputting ejlalb (albedo) for each wavelength
    call outfld('EJLALB  ',ejlalb, pcols,lchnk)
    call outfld('ALB325  ',alb325 ,pcols,lchnk)
    call outfld('ALB375  ',alb375 ,pcols,lchnk)
    call outfld('ALB425  ',alb425 ,pcols,lchnk)
    call outfld('ALB475  ',alb475 ,pcols,lchnk)
    call outfld('ALB525  ',alb525 ,pcols,lchnk)
    call outfld('ALB575  ',alb575 ,pcols,lchnk)
	call outfld('ALB642  ',alb642 ,pcols,lchnk)
	call outfld('ALB714  ',alb714 ,pcols,lchnk)
	call outfld('ALB784  ',alb784 ,pcols,lchnk)
	call outfld('ALB845  ',alb845 ,pcols,lchnk)
	call outfld('ALB891  ',alb891 ,pcols,lchnk)
	call outfld('ALB940  ',alb940 ,pcols,lchnk)
	call outfld('ALB985  ',alb985 ,pcols,lchnk)
	call outfld('ALB1070 ',alb1070 ,pcols,lchnk)
	call outfld('ALB1140 ',alb1140 ,pcols,lchnk)
	call outfld('ALB1220 ',alb1220 ,pcols,lchnk)
	call outfld('ALB1290 ',alb1290 ,pcols,lchnk)
    call outfld('ALB1380 ',alb1380 ,pcols,lchnk)
    call outfld('ALB1490 ',alb1490 ,pcols,lchnk)
	call outfld('ALB1610 ',alb1610 ,pcols,lchnk)
    call outfld('ALB1750 ',alb1750 ,pcols,lchnk)
	call outfld('ALB1910 ',alb1910 ,pcols,lchnk)
	call outfld('ALB2110 ',alb2110 ,pcols,lchnk)
	call outfld('ALB2350 ',alb2350 ,pcols,lchnk)


   return
end subroutine radcswmx
