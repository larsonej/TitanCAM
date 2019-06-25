#include <misc.h>
subroutine initindx
!----------------------------------------------------------------------- 
! 
! Purpose: Register constituents and physics buffer fields.
! 
! Author:    CSM Contact: M. Vertenstein, Aug. 1997
!            B.A. Boville, Oct 2001
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst, ppcnst, cnst_add, advected, nonadvec, cnst_chk_dim, cnst_name
  use phys_buffer,  only: pbuf_init
  use carma,	    only: carma_register
  use chemistry,    only: trace_gas, chem_register
  use cldcond,      only: cldcond_register
  use physconst,    only: mwdry, cpair, mwh2o, cph2o, mwch4, mwc2h2, mwc2h4, mwc2h6, mwhcn
  use tracers, only: tracers_register
!  use constituents, only: dcconnam, sflxnam, hadvnam, vadvnam, fixcnam, 
  use constituents, only: dcconnam, sflxnam, tendnam, tottnam
  use check_energy, only: check_energy_register
  use aerosol_intr, only: aerosol_register_cnst
! EJL 7-7-13 commented out aerosol register
  use abortutils,   only : endrun

#if ( defined BFB_CAM_SCAM_IOP )
  use iop
#endif
  implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!---------------------------Local variables-----------------------------
!
  integer m            ! loop index
  integer mm           ! constituent index 
!-----------------------------------------------------------------------

! Initialize physics buffer
  call pbuf_init()

! Register water vapor.
! ***** N.B. ***** This must be the first call to cnst_add so that
!                  water vapor is constituent 1.
  call cnst_add('Q', advected, mwh2o, cph2o, 1.E-12_r8, mm, &
                longname='Specific humidity', readiv=.true.)

  ! FAO: Add new constituents =====================================
  ! Cp is @ T=200K
  call cnst_add('CH4', nonadvec, mwch4, 2.087e3_r8, 0.0, mm, &
                longname='Methane', readiv=.true., mixtype='wet')
  call cnst_add('C2H2', nonadvec, mwc2h2, 1.47e3_r8, 0.0, mm, &
                longname='Acetylene', readiv=.true., mixtype='wet')
  ! Cp is @ T=175K
  call cnst_add('C2H4', nonadvec, mwc2h4, 1.241e3_r8, 0.0, mm, &
                longname='Ethylene', readiv=.true., mixtype='wet')
  ! Cp is @ T=250K
  call cnst_add('C2H6', nonadvec, mwc2h6, 1.535e3_r8, 0.0, mm, &
                longname='Ethane', readiv=.true., mixtype='wet')
  ! look up Cp later
  call cnst_add('HCN', nonadvec, mwhcn, 2.087e3_r8, 0.0, mm, &
                longname='Hydrogen_Cyanide', readiv=.true., mixtype='wet')

!
! Register cloud water
  call cldcond_register()
!
! Register chemical constituents
  if (trace_gas) then
     call chem_register()
  endif
!
! register aerosols
  call aerosol_register_cnst()

! aerosol microphysics
  call carma_register()

! Register advected test tracers and determine starting index
  call tracers_register()

!
! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim()
!
! Set default names for non-water advected and non-advected tracers
! Set names of advected and non-advected tracer diagnostics
!
  do m=1,ppcnst
     dcconnam(m) = 'DC'//cnst_name(m)
     sflxnam(m)  = 'SF'//cnst_name(m)
  end do
  do m=1,pcnst
!     hadvnam(m)  = 'HA'//cnst_name(m)
!     vadvnam(m)  = 'VA'//cnst_name(m)
!     fixcnam(m)  = 'DF'//cnst_name(m)
     tendnam(m)  = 'TE'//cnst_name(m)
     tottnam(m)  = 'TA'//cnst_name(m)
  end do

#if ( defined BFB_CAM_SCAM_IOP )
  do m=1,pcnst
     alphanam(m) = 'AFIX'//cnst_name(m)
     alphanam(m)=to_lower(alphanam(m))
     dqfxnam(m) = 'DQFX'//cnst_name(m)
     dqfxnam(m) = to_lower(dqfxnam(m))
  end do
#endif
  call check_energy_register()

end subroutine initindx
