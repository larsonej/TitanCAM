!------------------------------------------------------------------------
! File: initialize_model.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id: init_model.F90 62 2008-04-23 22:59:18Z cam_titan $
!------------------------------------------------------------------------
!  -- note that runtype.h is the same include file used in the ui code 
#include <runtype.h>
#include <params.h>
#include <max.h>
subroutine init_model( &
   gui_switches, &
   latitude, &
   longitude, &
   base_date, &
   base_secs, &
   runtype, &
   steplength, &
   restart, &
   guiseedval, &
   out_error_code, &
   AnalFile_, &
   ModelFile_, &
   IOPFile_, &
   LsmIniFile_, &
   LsmSurfFile_, &
   OzoneFile_, &
   PressFile_, &
   SicFile_, &
   SstFile_, &
   AbsEmsFile_, &
   AerOpticsFile_, &
   AerMassFile_, &
   LsmPftFile_, &
   UserFile_ &
   )

!=====================================================================
!-----------------------------------------------------------------------
!   
!     Interfaces with the GUI to set model run parameters
!
!---------------------------Code history--------------------------------
!
!     Written by John Truesdale    August, 1996
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use phys_grid
   use pmgrid,    only: masterproc, beglat, endlat, plat, plon, plev, plevp, dyngrid_set
   use comsrf
   use prognostics
   use buffer
   use comozp
   use shr_orb_mod
   use filenames
   use tracers
   use constituents, only: readtrace
   use chemistry,          only: trace_gas, chem_init
   use time_manager, only: calendar, dtime, nestep, nelapse, &
      start_ymd, start_tod, stop_ymd, stop_tod, ref_ymd, ref_tod, &
      perpetual_run, perpetual_ymd, tm_aqua_planet
   use history, only: scm_intht,outfld
   use iop
   use runtime_opts, only: runtime_options
   use scamMod
   implicit none
#if ( defined RS6000 )
   implicit automatic (a-z)
#endif

!------------------------------Includes---------------------------------
#include <comadj.h>
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comfrc.h>
!-----------------------------------------------------------------------
#include <comtfc.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!------------------------------Externals--------------------------------

!------------------------------Inputs-----------------------------------

   character*(MAX_PATH_LEN) IOPFile_
   character*(MAX_PATH_LEN) AnalFile_
   character*(MAX_PATH_LEN) UserFile_
   character*(MAX_PATH_LEN) SicFile_
   character*(MAX_PATH_LEN) ModelFile_
   character*(MAX_PATH_LEN) OzoneFile_
   character*(MAX_PATH_LEN) PressFile_
   character*(MAX_PATH_LEN) SstFile_
   character*(MAX_PATH_LEN) LsmSurfFile_
   character*(MAX_PATH_LEN) LsmIniFile_
   character*(MAX_PATH_LEN) AbsEmsFile_
   character*(MAX_PATH_LEN) AerOpticsFile_
   character*(MAX_PATH_LEN) AerMassFile_
   character*(MAX_PATH_LEN) LsmPftFile_

   integer base_date
   integer base_secs
   integer out_error_code
   integer steplength
   integer guiseedval
   logical restart
   integer runtype

   logical gui_switches( NUM_SWITCHES )  ! user and option switches
!  set in the GUI

!------------------------------Locals-----------------------------------
   real(r8) day  ! place holder
   real(r8) latitude
   real(r8) longitude

   integer i,m,k

   real(r8) hbuf
   integer lat            ! placeholders in outfld calls
   character*3 trnum         ! Advected species number

!-----------------------------------------------------------------------

!   dp1 = 0.
!   zp1 = 0.
!   psm2 = 0.
!   fsns = 0.

!============================================================
!     Fill in common blocks with the input parameters
!============================================================
!
!     Set comgui variables
!

   columnLat = latitude
   columnLon = longitude
!   write(6,*)'latitude=',latitude
!-----------------------------------------------------------------------
!

   call runtime_options()
   call t_initializef()
!
!
!-----------------------------------------------------------------------
!!$!
!!$!     Set comctl variables
!!$!
!!$   itsst = 1
!!$   iradlw = -1
!!$   iradsw = -1
!!$   iradae = -12
!!$   anncyc = .TRUE.
!!$
!!$   ldebug = .FALSE.
!!$   ozncyc = .TRUE.
!!$   sstcyc = .TRUE.
!!$   nlend =  .FALSE.

!
! NOTE: the following is only used in versions ccm3.3.20 and later
! 
! Initialize tracer related info 
! 
   trace_gas   = .false. ! greenhouse gas code not implemented
   readtrace   = .true. ! do not initialize from initial conditions file
!!$   scon        = 1.367e6 ! solar constant
   iyear_AD    = 1950    ! year AD to calculate the orbital parameters for.  
   anncyc = .TRUE.       ! annual cycling for sst
!!$!
!!$! Initialize volume mixing ratios - comvmr variables
!!$!
!!$   ch4vmr = 1.714e-6
!!$   n2ovmr = 0.311e-6
!!$   f11vmr = 0.280e-9
!!$   f12vmr = 0.503e-9
!!$   co2vmr = 3.550e-4
!!$
!
! Initialize Orbital parameters.
!
   obliq    = SHR_ORB_UNDEF_REAL ! Use default in orb_params
   eccen    = SHR_ORB_UNDEF_REAL ! Use default in orb_params
   mvelp    = SHR_ORB_UNDEF_REAL ! Use default in orb_params
!!$!
!!$! Initialize visible optical depth 
!!$!
!!$   tauvis = .14              
!
!     Set comtim variables
!
   seedval = guiseedval
   dtime = steplength

   start_ymd = base_date        ! base date of run 
   start_tod  = base_secs        ! base seconds of run

   iopTimeIdx = 1

   dyngrid_set = .true.
!
!-----------------------------------------------------------------------
   modelfile = ModelFile_
   iopfile = IOPFile_
   analysisfile = AnalFile_
   userfile = UserFile_
   sicfile = SicFile_
   ozonefile = OzoneFile_
   bndtvo=ozonefile
   pressfile = PressFile_
   sstfile = SstFile_
   bndtvs=sstfile
   lsmpftfile = LsmPftFile_
   lsmsurffile = LsmSurfFile_
   lsminifile = LsmIniFile_
   absems_data = AbsEmsFile_
   aeropticsfile = AerOpticsFile_
   aeroptics = AerOpticsFile_
   aermassfile = AerMassFile_
   bndtvaer = AerMassFile_
   absemsfile = AbsEmsFile_
!    
!  determine which type of run we're doing
!  for a saved initial conditions run, the value of the runtype is 
!  added to SIC, e.g., a saved initial conditions analysis run 
!  would be = ANAL + SIC
!
   use_userdata = .FALSE.
   use_saveinit = .FALSE.
   use_analysis = .FALSE.
   use_iop =      .FALSE.
   isrestart =    .FALSE.

   if ( runtype .GE. SIC ) then
      use_saveinit = .TRUE.
      runtype = runtype - SIC  
   endif

   if ( runtype .EQ. ANAL ) then
      use_analysis = .TRUE.
   endif

   if ( runtype .EQ. IOP ) then
      use_iop = .TRUE.
      use_analysis = .TRUE.
   endif

!
!  user data will be treated as if it were iop data, with the difference
!  that missing variables will be ignored
!
   if ( runtype .EQ. USER ) then
      use_analysis = .TRUE.
      use_userdata = .TRUE.
      use_iop = .TRUE.
      iopfile = userfile
   endif

!
!    set the latitude and longitude indexes into the datasets
!      
   call setLatLonIdx()


!
!     set the iop and option switches 
!      (the iop switches come after the option switches)
!
   isrestart=restart
   if (.not.isrestart) then
      do i=1, NUM_USER_SWITCHES
         switch(i) = gui_switches(i)
      end do

      use_srfprop = gui_switches(SRFPROP_SW+1)
      use_pert_init = gui_switches(PERT_INIT_SW+1)
      use_pert_frc = gui_switches(PERT_FRC_SW+1)
      use_relax   = gui_switches(RELAX_SW+1)
      use_3dfrc   = gui_switches(FRC3D_SW+1)
      fix_div3dfrc   = gui_switches(FIX_DIV3D_SW+1)
      use_diurnal_avg   = gui_switches(DIURNAL_AVG_SW+1)
   end if	
!
!-----------------------------------------------------------------------
!
! On restart reset the time indices
!
   n3   = 3
   n3m1 = 2
   n3m2 = 1

!
! Set up indexing for tracers - trcindx variables
!
   if (.not.isrestart) call initindx

! Initialize run control variables, time manager, timestep

   do m=1,pcnst
     alphanam(m) = 'AFIX'//cnst_name(m)
     dqfxnam(m) = 'DQFX'//cnst_name(m)
   end do
!
!     Read in boundary datasets to define physical initial conditions
!
!   call deallocate_for_restart()

     if(isrestart) then
      deallocate (ozlon)
      deallocate (ozlat)
      deallocate (date_oz)
      deallocate (sec_oz)
      deallocate (pin)
      deallocate (ozmixm(1)%val)
      deallocate (ozmixm(2)%val)
      deallocate (ozmix)
!!$      deallocate (ps)
!!$      deallocate (u3)
!!$      deallocate (v3)
!!$      deallocate (t3)
!!$      deallocate (q3)
!!$      deallocate (qminus)
!!$      deallocate (vort)
!!$      deallocate (div )
!!$      deallocate (dpsl)
!!$      deallocate (dpsm)
!!$      deallocate (dps )
!!$      deallocate (phis)
!!$      deallocate (omga)
!!$
      endif

!!$   write(6,*)'before inital t3=',t3
   call scam_inital(out_error_code )
   if ( out_error_code .ne. 0 ) return
!   write(6,*)'after inital t3=',t3,n3,n3m1,n3m2
!
! Set parameterization-specific constants
!
   if (.not.isrestart) call inti

!
! Initialize external models or datasets
!
   call initext
!   write(6,*)'after initext t3=',t3,n3,n3m1,n3m2
   if ( out_error_code .ne. 0 ) return


!
!     initialize history tape handler
!
   if (.not.isrestart)  then 
	call scm_intht
   else
      call outfld('PHIS',  phis, plond, lat)
      call outfld('PS',    ps(1,1,n3m2),   plond, lat)
      call outfld('Q',     q3(1,1,1,1,n3m2),   plond, lat)
      call outfld('T',     t3(1,1,1,n3m2),   plond, lat)
      call outfld('U',     u3(1,1,1,n3m2),   plond, lat)
      call outfld('V',     v3(1,1,1,n3m2),   plond, lat)
      call outfld('OMEGA', wfld, plond, lat)
      call outfld('DIVQ',  divq, plond, lat)
      call outfld('DIVT',  divt, plond, lat)
      call outfld('DIVQ3D',divq3d, plond, lat)
      call outfld('DIVT3D',divt3d, plond, lat)
      call outfld('DIVU',  divu, plond, lat)
      call outfld('DIVV',  divv, plond, lat)
   endif
!   write(6,*)'after intht t3=',t3,n3,n3m1,n3m2

!
!========================================================================
!
! Numerical scheme default values
!
   eps    = 0.06
   nlvdry = 3 
! 
!     Convert iradsw and iradlw to iterations if necessary
!     
   if (iradsw.lt.0) iradsw = nint((-iradsw*3600.)/dtime)
   if (iradlw.lt.0) iradlw = nint((-iradlw*3600.)/dtime)
!     
!     Convert iradae to iterations if necessary
!     
   if (iradae.lt.0) iradae = nint((-iradae*3600.)/dtime)
!     
!     iradae must be an even multiple of iradlw
!     
   if (mod(iradae,iradlw).ne.0) then
      print *, 'ERROR: initialize_model.F:iradae must be ', &
               'an even multiple of iradlw.'
      write(6,*)'     iradae = ',iradae,', iradlw = ',iradlw
      stop
   end if



!$$$JP      print * 
!$$$JP      print *,  'Column latitude = ',columnLat
!$$$JP      print *,  'Column longitude = ',columnLon
!$$$JP      print *,  ' use_analysis = ',use_analysis
!$$$JP      print *,  ' use_iop = ',use_iop
!$$$JP      print *,  ' use_saveinit = ',use_saveinit
!$$$JP      print *,  ' use_userdata = ',use_userdata
!$$$JP      print *,  ' base_date = ',base_date
!$$$JP      print *,  ' base_secs = ',base_secs
!$$$JP      print *,  ' steplen = ', steplength
!$$$JP      print *,  ' modelfile = ',modelfile
!$$$JP      print *,  ' analysisfile = ',analysisfile
!$$$JP      print *,  ' iopfile = ',iopfile
!$$$JP      print *,  ' userfile = ',userfile
!$$$JP      print *,  ' sicfile = ',sicfile
!$$$JP      print *,  ' ozonefile = ',ozonefile
!$$$JP      print *,  ' absemsfile = ',absems_data
!$$$JP      print *,  ' aeropticsfile = ',aeropticsfile
!$$$JP      print *,  ' aermassfile = ',aermassfile
!$$$JP      print *,  ' pressfile = ',pressfile
!$$$JP      print *,  ' sstfile = ',sstfile
!$$$JP      print *,  ' timeinvfile = ',timeinvfile
!$$$JP      print *,  ' lsmsurffile = ',lsmsurffile
!$$$JP      print *, ' use_srfprop = ', use_srfprop
!$$$JP      print *, ' use_pert_init = ', use_pert_init
!$$$JP      print *, ' use_pert_frc = ', use_pert_frc
!$$$JP      print *, ' use_relax   = ', use_relax
!$$$JP      print *, ' use_3dfrc   = ', use_3dfrc
!$$$JP      print *, ' use_diurnal_avg   = ', use_diurnal_avg
!$$$JP      print *  
   print *  
   print *  
   print *,'*******************************************************'
   print *  
   print *  
   print *,'********  SCAM MODEL INITIALIZATION COMPLETED  ********'
   print *  
   print *  
   print *,'*******************************************************'

   if ( use_ccmiop ) then
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,*) ' USING CCM GENERATED IOP DATA: ', iopfile
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * *'
   endif

   return
end subroutine init_model


