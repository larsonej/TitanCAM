#include <misc.h>
#include <params.h>

module tracers_suite

!---------------------------------------------------------------
!
! Implements suite of advected tracers
!
!  D Bundy June 2003
!  modified Feb 2004 to include TT_UN and smooting
!  modified Feb 2008 by E. Wilson for application in Titan model
!
! This is a module that contains a suite of tracers that is interfaced
! by tracers.F90. The details of the suite are contained entirely
! in this file, so the public routines are all very generic. 
!
!  ------------  calling tree --------------
!  Initialize the tracer mixing ratio field
!  	-> tracers.F90: tracers_init_cnst 
!  		-> tracers_suite.F90:init_cnst_tr
!
!  Initialize data set, things that need to be done at the beginning of a 
!  run (whether restart or initial)
!  	-> tracers.F90: tracers_init
!  		-> tracers_suite.F90:init_tr
!
!  Timestepping:
!  	-> tracers_timestep_init
!  		-> tracers_suite.F90:timestep_init_tr
!
!  	-> tracers_timestep_tend
!  		-> tracers_suite.F90:flux_tr
!
!  		-> tracers_suite.F90:tend_tr
!
!---------------------------------------------------------------


  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid !,       only: plon, plev, plat
  use ppgrid !,       only: pcols, pver
  use abortutils, only: endrun


  implicit none
  private
!  integer setpindxtr

  save

! Public interfaces
  public get_tracer_name  ! store names of tracers
  public get_tracer_char  ! stores characteristics of tracers
!  public init_cnst_tr      ! initialize tracer fields
  public init_tr      ! initialize data sets need for tracer
  public tend_tr      ! tracer tendency
  public flux_tr      ! surface flux of tracer
  public timestep_init_tr  ! interpolate tracer emissions data set

! Public module data
  integer, public, parameter :: trac_ncnst=18

! Private module data
  logical, parameter :: smooth = .false.
  
contains

!======================================================================
function get_tracer_name(n)

!------------------------------------------------------------------------
! Purpose:
!
! The tracer names are only defined in this module. This function is for
! outside programs to grab the name for each tracer number. 
!    If n > trac_ncst calls endrun.
!
!------------------------------------------------------------------------

! -----------------------------Arguments---------------------------------

  integer, intent(in) :: n
  character(len=8) :: get_tracer_name  ! return value
  character(len=8) :: name   ! constituent name
  integer :: type_adv
  real(r8)         :: mw, cp
  character(len=18) :: longname


!write(6,*) n
!write(6,*) trac_ncnst
  if ( n > trac_ncnst ) then
     write(6,*) 'tracers_suite:get_tracer_name()','requested tracer',n
     write(6,*) 'only ',trac_ncnst,' tracers available'
     call endrun
  else
!write(6,*) 'calling get_tracer_char'
     call get_tracer_char(n,name,type_adv,mw,cp,longname)
!write(6,*) 'return from get_tracer_char'
     get_tracer_name = name
  endif

  return

end function get_tracer_name


!======================================================================
!======================================================================
subroutine get_tracer_char(n,name,type_adv,mw,cp,longname)

!------------------------------------------------------------------------
! Purpose:
!
! To set the characteristics of tracers to be used in tracers_register
!    If n > trac_ncst calls endrun.
!
!------------------------------------------------------------------------
  use constituents, only: advected, nonadvec
  
  Implicit none
  
! -----------------------------Arguments---------------------------------

  integer, intent(in) :: n
  character(len=8) :: name
  character(len=*) :: longname
  integer, intent(out) :: type_adv      !advected or nonadvec
  real(r8)         :: mw, cp

! ----------------------------- Local ---------------------------------
  character(len=8), dimension(trac_ncnst), parameter :: & ! constituent names
	trac_name  =  (/ 'CH4TR', 'H2TR', 'C2H2TR', 'C2H4TR' , &
	'C2H6TR' , 'CH3C2HTR' , 'CH2CCH2' , 'C3H6TR' , 'C3H8TR' ,&
	'C4H2TR' , 'C6H2TR' , 'HCNTR' , 'HC3NTR' , 'C2H3CNTR' , &
	'HAZE_0', 'HAZE_A', 'HAZE_B', 'HAZE_C' /)
  integer, dimension(trac_ncnst), parameter :: & ! advection type
	trac_type  =  (/ advected, advected, advected, advected , &
	advected , advected , advected , advected , advected ,&
	advected , advected , advected , advected , advected , &
	nonadvec, nonadvec, nonadvec, nonadvec /)
	
  character(len=*), dimension(trac_ncnst), parameter :: &!constituent longnames
	trac_longname = (/ 'Methane' , 'Molecular_Hydrogen' , &
	'Acetylene' , 'Ethylene' , 'Ethane' , 'Methylacetylene' , &
	'Allene' , 'Propylene' , 'Propane' , 'Diacetylene' , &
	'Triacetylene' , 'Hydrogen_Cyanide' , 'Cyanoacetylene' , &
	'Acrylonitrile' , &
	'HAZE_0', 'HAZE_A', 'HAZE_B', 'HAZE_C' /)
  real(r8), dimension(trac_ncnst), parameter :: & !constituent mol. weights
	trac_mw = (/ 16.,2.,26.,28.,30.,40.,40.,42.,44.,50.,74., &
	29.,51.,53. , &
	2.5e12, 2.5e12, 2.5e12, 2.5e12 /) 
  real(r8), dimension(trac_ncnst), parameter :: & !constituent specific heats
	trac_cp = (/ 2.087e3_r8,1.531e4_r8,1.47e3_r8,1.241e3_r8, &
	1.535e3_r8,1.061e3_r8,1.021e3_r8,1.056e3_r8,1.109e3_r8, &
	0.947e3_r8,0.895e3_r8,1.445e3_r8,1.445e3_r8,1.445e3_r8, &
	1.445e3_r8, 1.445e3_r8, 1.445e3_r8, 1.445e3_r8 /)
!-----------------------------------------------------------------------

         if ( n > trac_ncnst ) then
            write(6,*) 'tracers_suite:get_tracer_name()','requested tracer',n
            write(6,*) 'only ',trac_ncnst,' tracers available'
            call endrun
         else
!write(6,*) n
            name = trac_name(n)
!write(6,*) name
!write(6,*) trac_longname
            longname = trac_longname(n)
!write(6,*) longname          
            mw = trac_mw(n)
!write(6,*) mw
            cp = trac_cp(n)
!write(6,*) cp
            type_adv = trac_type(n)
!write(6,*) type_adv
         endif	

  return

end subroutine get_tracer_char
  
!-----------------------------------------------------------------------
!======================================================================
subroutine init_cnst_tr(m,q)

!----------------------------------------------------------------------- 
!
! Purpose:
! calls initialization routine for tracer m, returns mixing ratio in q
!
! This routine must be consistent with trac_ncnst.
!
!----------------------------------------------------------------------- 

   implicit none

   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air
   integer ,intent(in) :: m                      ! index of tracer

   if ( m > trac_ncnst ) then
      write(6,*) 'tracers_suite:init_cnst_tr()'
      write(6,*) ' asked to initialize tracer number ',m
      write(6,*) ' but there are only trac_ncnst = ',trac_ncnst,' tracers'
      call endrun
   endif

      write(6,*) 'tracers_suite:init_cnst_tr()'
      write(6,*) 'no initialization routine specified for tracer',m
      call endrun
      
end subroutine init_cnst_tr



!======================================================================

subroutine init_tr

  !----------------------------------------------------------------------- 
  ! Purpose:
  !-----------------------------------------------------------------------


!---------------------------Local workspace-----------------------------
    integer :: ix                                ! tracer index
!-----------------------------------------------------------------------


end subroutine init_tr

!======================================================================

subroutine timestep_init_tr
!----------------------------------------------------------------------- 
! 
! Purpose: call tracer timestep init processes
!----------------------------------------------------------------------- 


end subroutine timestep_init_tr

!======================================================================

subroutine setsmoothtr(indx,q,width,rev_in)

  implicit none

#include <comhyb.h>

  !Arguments
  integer, intent(in)     :: indx               ! k index of pressure level
  real(r8), intent(inout) :: q(plon,plev,plat)  ! kg tracer/kg dry air
  real(r8), intent(in)    :: width              ! eta difference from unit level where q = 0.1
  integer,  intent(in), optional :: rev_in      ! reverse the mixing ratio


  !Local variables
  integer k
  real(r8) alpha ! guassian width, determined by width, T
  real(r8) pdist ! pressure distance (eta.e4) from k=indx
  real(r8) T     ! desired m.r. in level specified by pdiff from k=indx
  integer rev  ! = 1 then reverse (q = 1, q(k=indx) = 0 )




  rev = 0
  if (present(rev_in)) then
     if (rev_in == 1) then
        rev = 1
     endif
  endif


  write(6,*)'DRB TR SMOOTH indx hypm(indx)',indx,hypm(indx)
  write(6,67)'DRB TR SMOOTH ','k','hypm(k)','pdist','-a*(pd^2)','e^-a*(pd^2)'

  T = 0.1
  alpha = -log(T)/(width*1.e4)**2  ! s.t. in level width from indx, mr = T

!  alpha = 3.e-8  ! m.r. ~ 0.1 in adjacent levels, where change eta ~ 0.08

  do k=1,pver
     pdist = hypm(k) - hypm(indx)

     if ( rev == 1 ) then
        q(:,k,:) = 1.0 - exp(-alpha*(pdist**2))
     else
        q(:,k,:) =       exp(-alpha*(pdist**2))
  endif
     write(6,66)'DRB TR SMOOTH ',k,hypm(k),pdist,-alpha*pdist**2,q(1,k,1)
  end do
  
66 format (a15,i4,3f15.2,g15.6)
67 format (a15,a4,4a15)
  
!  call endrun('tracers_suite.F90 L278 DRBDBG')

end subroutine setsmoothtr


!======================================================================

integer function setpindxtr(pmb)

  ! find the index of layer nearest pmb

  use pmgrid, only: plev, plevp

  implicit none
  
#include <comhyb.h>
  
  real(r8) pmb
  real(r8) pmin, pdist
  integer indx, k
  
  indx = 0
  pmin = 1.e36
  pdist = 1.e36
  do k=1,pver
     pdist = abs(hypm(k) - pmb*100.)
     if (pdist < pmin) then
        indx = k
        pmin = pdist
     end if
  end do
  write (6,*) ' index near', pmb, ' is ', indx
  setpindxtr = indx

end function setpindxtr

!======================================================================
!======================================================================

subroutine tend_tr(m,ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of radon tracer
! 
! Method:
!-------------------------Code History----------------------------------
!
! tracers.F90:rndecay()
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
!
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: m                 ! tracer number
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! radon mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                             ! (kg tracer /(s kg moist air))
!--------------------------Local Variables------------------------------

  if ( m > trac_ncnst ) then
      write(6,*) 'tracers_suite:tend_tr()'
      write(6,*) ' asked to calculate tendency for tracer number ',m
      write(6,*) ' but there are only trac_ncnst = ',trac_ncnst,' tracers'
      call endrun('tracers_suite.F90:tend_tr() L457')
   endif



end subroutine tend_tr

!======================================================================

subroutine flux_tr(m,ncol, lchnk, landfrac, flux )
!----------------------------------------------------------------------- 
!

!-----------------------------------------------------------------------
   implicit none
!--------------------------Arguments-------------------------------------

   integer,  intent(in)  :: m               ! tracer number
   integer,  intent(in)  :: ncol            ! number of atmospheric columns in chunk
   integer , intent(in)  :: lchnk           ! current identifier
   real(r8), intent(in)  :: landfrac(pcols) ! landfraction
   real(r8), intent(out) :: flux(pcols)     ! specified radon flux in kg/m^2/s

!--------------------------Local Variables------------------------------


   if ( m > trac_ncnst ) then
      write(6,*) 'tracers_suite:flux_tr()'
      write(6,*) ' asked to calculate flux for tracer number ',m
      write(6,*) ' but there are only trac_ncnst = ',trac_ncnst,' tracers'
      call endrun
   endif
   
   ! flux is null for all tracers
   flux = 0.


end subroutine flux_tr


end module tracers_suite
