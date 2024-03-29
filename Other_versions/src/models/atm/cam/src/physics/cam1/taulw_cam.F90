
! taulw_cam:
!  1. Combines the ck coefficients of multiple gases into single one
!     appropriate for given mixture in atmospheric column
!  2. Computes ck, cia optical thicknesses, combines them 
!  3. Calculates longwave optical depth for each atmospheric column
!
! authors:
!   David Ye/Hamik Mukelyan/AJF

! AJF, 3-20-08, modified for CAM -- memory allocations moved to radcnst
! --------------------------------------------------------------------------
!    **NOTE:  This version assumes get_taulw_cam will be called within
!             loop over columns from inside radclwmx_titan
! --------------------------------------------------------------------------

Module taulw_cam

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: masterproc
!  use radcnst, only: radcnst_initialize, indexx

Contains

Subroutine get_taulw_cam(mmr_tr, pint, tlay, tau_tot, debug_rad)

! AJF, 3-20-08:  Combines Ck coefficients of different gases together,
!  combines the effective Ck absorption with Cia absorption in each
!  wavenumber interval, and computes optical depth array for each
!  atmospheric column separately.  CkComb contains the combined-gas Ck coefs.
!  CiaComb contains the CIA coefs interpolated to the CAM wavenumber-temperature 
!  grid.

  use radcnst,  only: N_gas_ck, Ng_ck, Nwn_ck, Ck, Ck_m, CkComb, CiaComb, Cia_v, &
                      Nv_p1, num_g, Ng_lw, gas_name_lw, wn_lw, Nf_lw, &
                      Where,N_cia,GasCk,GasCia,n_type_aer,aer_name_lw, &
                      k_aer_lw, Nwn_Ckcomb, Nwn_Cia
  
  use constituents,    only: cnst_get_ind, cnst_mw, ppcnst
  use ppgrid,          only: pver

  Implicit None
  
!  Arguments:

   real(r8), intent(in) :: mmr_tr(pver,ppcnst)        !mass-mixing ratios(lev,species)
   real(r8), intent(in) :: pint(:)                 !interface pressure, Pa
   real(r8), intent(in) :: tlay(:)                  !layer temperature, K
   real(r8), intent(out) :: tau_tot(Nf_lw, Ng_lw, Nv_p1) !Total optical depth, CK+CIA
   logical  :: debug_rad
    
!  Local variables:

  real(r8), parameter :: Mp = 1.67e-24
  real(r8) :: mixr

  Integer :: TempPres, Wavenum
  Integer :: Nwn
  Integer :: I, J, K, ind(N_gas_ck)
  integer :: i_cia_first  !Index of N2-N2 CIA coef. For Titan, CIA coefs run
                          ! from i_cia_first to N_cia. ***Needs to be changed 
                          ! for gas giants***
                          
  real(r8) :: tauCia(N_cia,Nv_p1,Nwn_Cia)  ! Column optical depth due to CIA
  real(r8) ::  tauCk(Nwn_Ckcomb,Ng_lw,Nv_p1)  ! Column optical depth, CK
  real(r8) ::  taulw_aer(Nf_lw,Nv_p1) ! LW column optical depth, aerosols
  integer  :: tfid  !logical unit for output
  real(r8) :: temp_r, pint_bars(Nv_p1)
  real(r8) :: k_aer_lw_comb(pver,Nf_lw), sum_aero
  integer  :: ind_aer(n_type_aer) 
  
    
!--------------------------------------------------------------
  
  Cia_v%Val=CiaComb%Val
     
  Nwn = Nwn_Ckcomb
  
    tauCk(:,:,:)=0.0
    tauCia(:,:,:)=0.0
    taulw_aer(:,:)=0.0
    tau_tot(:,:,:)=0.0
  
!  Get indices of Ck longwave-active gases:

   do K=1,N_gas_ck
     call cnst_get_ind (Trim(gas_name_lw(K)), ind(K))
   enddo
    

! Multiply the Ck values by the gas mass-mixing ratios:
! (needed before combining into effective k for optical depth calc.  Cks are
!  divided by Mp*molec-weight*1.e20 to convert from cm^2/(molecule)*1.e20 to 
!  cm^2/g. The factor of 1.e20 comes from the original Oxford definition) 

  Do K = 1, N_gas_ck
       
    J=1  !Use mixing ratio of layer 1 for top interface
      
     Ck_m(K)%Val(:, J, :)=Ck(K)%Val(:, J, :)*mmr_tr(J,ind(K)) / &  
                          (1.e20*Mp * cnst_mw(ind(K)))      
    Do J = 2, Nv_p1 - 1
     mixr=0.5*(mmr_tr(J-1,ind(K))+mmr_tr(J,ind(K)))  !avg to interface
     Ck_m(K)%Val(:, J, :)=Ck(K)%Val(:, J, :)*mixr/(1.e20*Mp * cnst_mw(ind(K)))
    Enddo
    
    J=Nv_p1 !Use mixing ratio of layer Nv_p1-1 for lowest interface
         
     Ck_m(K)%Val(:, J, :) = Ck(K)%Val(:, J, :)*mmr_tr(J-1,ind(K)) / &  
                          (1.e20*Mp * cnst_mw(ind(K)))
  End Do
  

! Mix the Ck coefficients of the gases to make one effective Ck stored in CkComb

      CkComb%Val(:,:,:)=0.0

  Do Wavenum = 1, Nwn
     
! Mix gas ck's for gases having a nonzero ck at each wavenumber, one level at a time

      Do TempPres = 1, Nv_p1
        
         Do K = 1, N_gas_ck
         
!!!        if (TempPres == 1) then
!!!         if (masterproc) then
!!!          print*,' k (gas) = ',k
!!!          print*,' Ck wavenumber array = '
!!!          do j=1,size(Ck(K)%Wn)
!!!           print*,'j= ',j,' Wn= ',Ck(K)%Wn(j),' Ck%Val= ',Ck(K)%Val(1,TempPres,j)
!!!          enddo
!!!         endif
!!!        endif
        
         If (Where(Ck(K)%Wn, CkComb%Wn(Wavenum)) /= 0) Then
           Call MixGas(CkComb, Ck_m(K), TempPres, CkComb%Wn(Wavenum))
         End If
        
!!!        if (TempPres == 1) then
!!!         if (masterproc) then
!!!          print*,' k (gas) = ',k
!!!          print*,' CkComb%Wn = ',CkComb%Wn(Wavenum),' Val= ',&
!!!                   CkComb%Val(1,TempPres,Wavenum)         
!!!         endif
!!!        endif
 
         End Do
         
      End Do
  End Do
    
! Multiply Cia coefs by vmr's of collision partners
    
    call process_cias(Cia_v, mmr_tr)
    
! Multiply aerosol absorption coefs by aerosol mass-mixing ratios, and sum to
! get total effective coef as function of pressure and wavenumber
 

    do k=1,n_type_aer
     call cnst_get_ind (Trim(aer_name_lw(k)), ind_aer(k))
    enddo

    do i=1,Nf_lw
     do j=1,Nv_p1-1  !mid-layers only
      sum_aero=0.0
       do k=1,n_type_aer
        sum_aero=sum_aero+k_aer_lw(k)%Val(i)*mmr_tr(j,ind_aer(k))
       enddo
      k_aer_lw_comb(j,i) = sum_aero
     enddo
    enddo
      
! Convert pint from Pa to atms for use with init_(x)_OD_table routines

    pint_bars(:)=pint(:)*1.0e-5
    
! Produce CIA optical depths on grid for column i

    call init_cia_OD_table (tauCia, Cia_v, CkComb, pint_bars, tlay, i_cia_first)
    
! Produce Ck optical depths on grid for column i

    call init_ck_OD_table (tauCk, CkComb, num_g, Ng_lw, pint_bars, tlay)
    
! Produce lw aerosol optical depths on grid for column i

    call init_aer_lw_OD (taulw_aer, k_aer_lw_comb, pint_bars, tlay)
        
! Combine Ck and Cia optical depths with full wavenumber grids.  "tau_tot"
! is total lw optical depth to be passed to radiative transfer routines.

   call combine_cia_ck_aer (tauCia, tauCk, taulw_aer, tau_tot, & 
                        Cia_v, CkComb, wn_lw, Ng_lw, num_g, i_cia_first)
                                                
                                                
  if (debug_rad) then
                        
!  Diagnostic Output (generally keep turned off)

  !--------------------v
  ! Output CIA tau's in Jim's format
  tfid=14
  open(tfid, File = 'CIA_taus', Status = 'replace', Err = 7)
  write(tfid, *, Err = 7) 'TOTAL CIA OPTICAL DEPTH AT EACH LAYER'
  do i=1, size(CiaComb%Wn)
     write(tfid, *, Err=7) 'Wavenumber =', CiaComb%wn(i), 'cm^-1'
     write(tfid, '(A15,A15,A15)', Err=7) 'PRES (B)', 'TAU', '---'
     do j=1, Nv_p1

        temp_r=0
        ! summing only last four gases b/c others negligible on titan
        do k=7, N_cia
           temp_r=temp_r+tauCia(k,j,i)
        end do

        write (tfid,'(E15.5,E15.5,A15)', Err=7) CkComb%P(j), temp_r, '---'  
     end do
  end do
  close(tfid)
  
  !--------------------^
  
  
  !--------------------v
  ! Output CK tau's in Jim's format
  tfid=14
  open(tfid, File ='CK_taus', Status = 'replace', Err = 7)
  write(tfid, *, Err = 7) 'TOTAL CK OPTICAL DEPTH AT EACH LAYER'
  do i=1, size(CkComb%Wn)
     do k=1, num_g(i)
        write(tfid, *, Err=7) 'Wavenumber =', CkComb%wn(i), 'cm^-1'
        write(tfid, *, Err=7) 'G-val =', CkComb%G(k)
        write(tfid, '(A15,A15,A15)', Err=7) 'PRES (B)', 'TAU', '---'
        do j=1, Nv_p1
           write (tfid,'(E15.5,E15.5,A15)', Err=7) CkComb%P(j), tauCk(i,k,j),'---' 
        end do
     end do
  end do
  close(tfid)
  
      
  !--------------------^
  
  !--------------------v
  ! Output aerosol LW tau's in Jim's format
  tfid=14
  open(tfid, File = 'AERLW_taus', Status = 'replace', Err = 7)
  write(tfid, *, Err = 7) 'TOTAL LW AEROSOL OPTICAL DEPTH AT EACH LAYER'
  do i=1, Nf_lw
     write(tfid, *, Err=7) 'Wavenumber =', wn_lw(i), 'cm^-1'
     write(tfid, '(A15,A15,A15)', Err=7) 'PRES (B)', 'TAU', '---'
     do j=1, Nv_p1
        write (tfid,'(E15.5,E15.5,A15)', Err=7) CkComb%P(j), taulw_aer(i,j),'---'  
     end do
  end do
  close(tfid)
  

  !--------------------v
  ! Output combined tau's in Jim's format
  tfid=14
  open(tfid, File ='TOT_taus', Status = 'replace', Err = 7)
  write(tfid, *, Err = 7) 'TOTAL OPTICAL DEPTH AT EACH LAYER'
  do i=1, size(wn_lw)
     do k=1, num_g(i)
        write(tfid, *, Err=7) 'Wavenumber =', wn_lw(i), 'cm^-1'
        write(tfid, *, Err=7) 'G-val =', CkComb%G(k)
        write(tfid, '(A15,A15,A15)', Err=7) 'PRES (B)', 'TAU', '---'
        do j=1, Nv_p1
           write (tfid,'(E15.5,E15.5,A15)', Err=7) CkComb%P(j), tau_tot(i,k,j),'---' 
        end do
     end do
  end do
  close(tfid)
  
  !--------------------^ 
   
  
  endif
  

  return
  
7   print*, 'Failed to open or write to file: <CIA or CK or TOT>_taus'
    stop
End subroutine get_taulw_cam

!###########################################################################



  !MixGas:
  !  Combines the k values for two gases at a shared wavenumber interval. The 
  !  twenty k values and twenty weights are combined to form four hundred of 
  !  each. The four hundred k values are sorted and twenty new k values are 
  !  linearly interpolated at the gaussian points.
  !
  !  Arguments:
  !     Ck1   : First gas Ck
  !     Ck2   : Second gas Ck
  !     TP    : Pressure/Temperature index
  !     WnVal : Wavenumber interval
  
  ! AJF 4/08/08: Changed David Ye's original version.  Must sort combined weights
  !              simultaneously with cks and create new x-axis for sorted
  !              combined-ck set

  Subroutine MixGas(Ck1, Ck2, TP, WnVal)
  
    Use radcnst,  only:  GasCk, wght1d_lw, Where, Sort, Interpolate 
    Implicit None

    Type(GasCk), Intent(InOut) :: Ck1
    Type(GasCk), Intent(In)    :: Ck2
    Integer,     Intent(In)    :: TP
    real(r8),     Intent(In)    :: WnVal

    ! WComb, CkComb : Combined k and weight array
    ! W             : Array of weights
    ! I, J, K       : Loop indicies
    ! Wn1, Wn2      : Wavenumber indicies

    real(r8), Dimension (Size(Ck1%G)**2) :: &
             WComb,       & ! Unsorted combined weights
             CombK,       & ! Unsorted combined CKs
             W_s,         & ! Sorted combined weights (*not* in ascending order)
             CombK_s,     & ! Sorted combined CKs (in ascending order)
             g_s            ! New g values for sorted-Ck x-axis
    real(r8)   :: sumg
    integer, Dimension (Size(Ck1%G)**2) :: indx            
    real(r8), Dimension (Size(Ck1%G))      :: W
    Integer                            :: I, J, K
    Integer                            :: Wn1, Wn2
    Integer                            :: Nm   !Total number of combined wts, cks

!!!    if( masterproc ) print*,'MIXGAS:'

    W = wght1d_lw  !Standard (unchanging) gaussian weights for gas Ck's

    ! Return the index in the wavenumber array that equals the wavenumber
    Wn1 = Where(Ck1%Wn, WnVal)
    Wn2 = Where(Ck2%Wn, WnVal)
    
!!!    if (TP == 1) then
!!!    if (masterproc) then
!!!     print*,' Ck1%Wn= ',Ck1%Wn(Wn1),' Ck2%Wn= ',Ck2%Wn(Wn2)
!!!     print*,' Size(Ck1)= ',size(Ck1%Val,1),size(Ck1%Val,2),size(Ck1%Val,3)
!!!     print*,' Size(Ck2)= ',size(Ck2%Val,1),size(Ck2%Val,2),size(Ck2%Val,3)
!!!    endif
!!!    endif

    ! In the overlapping wavenumber interval, combine the k and weight values of
    ! two gases. Each k value of one gas is combined with all k-values of other
    ! gas. The weights are multiplied.

    K=0
    Do I = 1, Size(W)
    Do J = 1, Size(W)
      K = K + 1
      CombK(K) = Ck1%Val(I, TP, Wn1) + Ck2%Val(J, TP, Wn2)
      WComb(K) = W(I) * W(J)
    End Do
    End Do
    Nm=K
    
!!!    if (TP == 1) then
!!!    if (masterproc) then
!!!     print*,' Nm= ',Nm
!!!     print*,' Size(Wcomb,CombK)= ',size(WComb),size(CombK)
!!!    endif
!!!    endif

    ! Sort the combined k values and weights, with k's in ascending order:
    
     call indexx(Nm,CombK,indx)
     
    ! Rearrange CKs, combined weights with CKs in ascending order
    
     do i=1,Nm
      CombK_s(i) = CombK(indx(i))
      W_s(i)     = WComb(indx(i))
     enddo
     
  ! Find new gaussian abscissa points corresponding to the sorted combined Ck values:
    
     sumg=0.0
     do i=1,Nm
      g_s(i)=W_s(i)+sumg
      sumg=sumg+W_s(i)
     enddo
     
!!!    if (TP == 1) then
!!!    if (masterproc) then
!!!     do i=1,Nm
!!!      print*,'i= ',i,' g_s= ',g_s(i),' CombK_s= ',CombK_s(i),' W_s(i)= ',W_s(i)
!!!     enddo
!!!    endif
!!!    endif    
     
  ! Interpolate at the gaussian points if the k value of Ck1 is not zero.
    Do I = 1, Size(W)
       If (Ck1%Val(I, TP, Wn1) /= 0.0) Then
          Ck1%Val(I, TP, Wn1) = Interpolate(g_s, CombK_s, Ck1%G(I))
       Else
          Ck1%Val(I, TP, Wn1) = Ck2%Val(I, TP, Wn2)
       End If
!!!    if (TP == 1) then
!!!    if (masterproc) then
!!!     print*,' i= ',i,' Ck1%Val= ',Ck1%Val(i,TP,Wn1)
!!!    endif
!!!    endif
    End Do
    
  return
  End Subroutine MixGas
  
!#############################################################################

! Process_cias 
! Author     : Hamik Mukelyan
! Date       : 8/24/07
! AJF, 3-25-08: modified for inclusion in CAM model.  Exponentiation of
!               cia coef data (i.e., of ln[input table values]) moved to
!               interpcia_cam

! Description: Multiplies CIA coefficients by volume mixing ratios.
  
  !--------------------V
  ! Description:  multiply Cia coefficients by volume mixing ratios.
  !
  ! In: cias -- Cia coefficient data structure. See table.f90 
  !     chi  -- Mass mixing ratio of methane
  !
  ! Out: Process cias
  subroutine process_cias (cias, chi)
  
  use radcnst,         only: N_cia, Nv_p1, GasCia
  use constituents,    only: cnst_get_ind, cnst_mw
  use shr_const_mod,   only: SHR_CONST_MWDAIR

  implicit none

  real(r8), parameter :: mole_fraction_h2=1.e-3
  
! Arguments:  
    type(GasCia), intent(inout) :: cias
    real(r8), intent(in)        :: chi(:,:)  !mass-mixing-ratio profiles in column 
    
! Local Variables:
    integer                     :: i,j,k
    real(r8)                    :: temp_r, mw_ch4
    real(r8)                    :: mole_fraction_n2(Nv_p1),vmr_ch4(nv_p1)
    
! Compute mole fractions of N2 and CH4

      call cnst_get_ind('CH4',k)
      mw_ch4=cnst_mw(k)

      j=1     !Assume ch4-vmr at top interface = vmr of topmost layer
      vmr_ch4(j) = chi(j,k)*SHR_CONST_MWDAIR/mw_ch4
      mole_fraction_n2(j) = 1.0-vmr_ch4(j)-mole_fraction_h2
      
      do j=2,Nv_p1-1 ! avg layers to interfaces
        vmr_ch4(j) = 0.5*(chi(j-1,k)+chi(j,k))*SHR_CONST_MWDAIR/mw_ch4
        mole_fraction_n2(j) = 1.0-vmr_ch4(j)-mole_fraction_h2
      enddo
      
      j=Nv_p1  !assume ch4-vmr at surface interface = that of bottom layer
      vmr_ch4(j) = chi(j-1,k)*SHR_CONST_MWDAIR/mw_ch4
      mole_fraction_n2(j) = 1.0-vmr_ch4(j)-mole_fraction_h2
    
! Multiply the cia vals by the corresponding volume mixing ratios. Assume 
! cia file is in this order:
    !  1 equilibrium H2-H2 (eH2-H2)
    !  2 normal H2-H2 (nH2-H2)
    !  3 eH2-He
    !  4 nH2-He
    !  5 eH2-CH4
    !  6 nH2-CH4
    !  7 N2-N2
    !  8 N2-H2
    !  9 N2-CH4
    ! 10 CH4-CH4
    ! Note that we only care about the last 4 combinations since Titan's atmosphere 
    ! contains negligible quantities of the others. 
    
   do k=1, Nv_p1
    do i=7, N_cia
       if (i==7) then
          temp_r = mole_fraction_n2(k) * mole_fraction_n2(k)
       end if
       if (i==8) then
          temp_r = mole_fraction_n2(k) * mole_fraction_h2
       end if   
       if (i==9) then
          temp_r = mole_fraction_n2(k) * vmr_ch4(k)
       end if
       if (i==10) then
          temp_r = vmr_ch4(k)*vmr_ch4(k)
       end if
         do j=1, size(cias%Wn)
          cias%Val(i,k,j) = cias%Val(i,k,j) * temp_r
         enddo
    end do
   end do
   
   return
  end subroutine process_cias
  
!#############################################################################

  !--------------------v  
  ! init_cia_OD_table
  ! Description: Initializes CIA optical depth table.
  ! Originated by Hamik Mukelyan, Aug 2007
  ! AJF, 3/25/08: Modified for CAM model
  !
  ! In: tauCia     -- (collision pair, vertical level, wavenumber). CIA optical depths
  !     cias       -- cia coefficient data structure
  !     ck         -- ck coef data struct (needed for pressures)
  !     pint       -- pressure at interface (bars)
  !     tlay       -- mid-layer temperature (K)
  !     
  ! Out: none
  
  subroutine init_cia_OD_table (tauCia, cias, ck, pint, tlay, i_cia_first)
  
    use radcnst,   only: N_cia, Nv_p1, GasCia, GasCk
    use ppgrid,          only: pver
    
    implicit none
    
!  Arguments
    real(r8), intent(inout) :: tauCia(:,:,:)
    real(r8), intent(in)           :: pint(Nv_p1)
    real(r8), intent(in)           :: tlay(pver)   
    type(GasCia), intent(in)       :: cias
    type(GasCk), intent(in)        :: ck
    integer, intent(out) :: i_cia_first
    
! Local variables:    
    integer                        :: i,j,k
    integer, parameter :: i_cia_firstp=7  ! Index of N2-N2 interaction
                                          ! = first one in list to use
                                          ! for Titan, last is = N_cia

    i_cia_first=i_cia_firstp
       
    do i=i_cia_first, N_cia
       do j=1, size(cias%Wn)
         
          tauCia(i,1,j) = optical_thickness (cias%Val(i,1,j), &
               0.0, pint(1), tlay(1), 'a')  !ajf, 12/13/07, 3-27-08
         
          do k=2, Nv_p1         
             tauCia(i,k,j) = &
                  optical_thickness (0.5*(cias%Val(i,k-1,j)+cias%Val(i,k,j)), & 
                  pint(k-1), pint(k), tlay(k-1), 'a')
             tauCia(i,k,j)=tauCia(i,k-1,j)+tauCia(i,k,j) ! add to layers above
          end do
       end do
    end do

    return    
  end subroutine init_cia_OD_table
  
!###############################################################################
  
  !--------------------v
  ! init_ck_OD_table 
  ! Description: Initializes CK optical depth table.
  !
  ! In: tauCk      -- (wavenumber,gauss point,level). Optical Depths at interfaces
  !     ck         -- ck coefficient data structure
  !     nsa        -- number of gaussian points at each spectral interval
  !     nsm        -- max num g points per spectral interval
  !     pint       -- interface pressures, bars
  !     tlay       -- mid-layer temperature, K
  ! 
  ! Out: none
  subroutine init_ck_OD_table (tauCk, ck, nsa, nsm, pint, tlay)
  
  use radcnst, only: Nv_p1, GasCk
  use ppgrid,          only: pver
  
  implicit none
  
! Arguments:
    real(r8), intent(inout) :: tauCk(:,:,:)
    real(r8), intent(in)           :: pint(Nv_p1)
    real(r8), intent(in)           :: tlay(pver)     
    type(GasCk), intent(in)        :: ck
    integer, intent(in)            :: nsa(:)
    integer, intent(in)            :: nsm
    integer                        :: i,j,k

    
    do i=1, size(ck%Wn)
       do j=1, nsa(i)
       
       
          tauCk(i,j,1)=optical_thickness(ck%Val(j,1,i), &
               0.0, pint(1), tlay(1), 'k')  !optical depth above top interface, ajf 12/13/07
          

          do k=2, Nv_p1 ! first layer to last
             ! get optical thickness of layer k-1 
             tauCk(i,j,k) = &
                  optical_thickness(0.5*(ck%Val(j,k-1,i)+ck%Val(j,k,i)),&
                  pint(k-1), pint(k), tlay(k-1),'k')
             tauCk(i,j,k) = tauCk(i,j,k-1) + tauCk(i,j,k) ! add to layers above
          end do
          
       end do
    end do
    
    return
  end subroutine init_ck_OD_table
  
!################################################################################

 !--------------------v
  ! init_aer_lw_OD 
  ! Description: Initializes lw aerosol optical depth table.
  !
  ! In: 
  !     k_aer_lw_comb -- effective aerosol abs coefficient summed over aer types,
  !                      weighted by type mmr's
  !     pint       -- interface pressures, bars
  !     tlay       -- mid-layer temperature, K
  ! 
  ! Out: taulw_aer      -- aerosol lw optical depths at interfaces (wn,vert lev)
  subroutine init_aer_lw_OD (taulw_aer, k_aer_lw_comb, pint, tlay)
  
  use radcnst, only: Nv_p1, Nf_lw
  use ppgrid,          only: pver
  
  implicit none
  
! Arguments:
    real(r8), intent(inout) :: taulw_aer(:,:)
    real(r8), intent(in)           :: pint(Nv_p1)
    real(r8), intent(in)           :: tlay(pver)     
    real(r8), intent(in)        :: k_aer_lw_comb(:,:)
    integer                        :: i,k

    
    do i=1, Nf_lw
          !(assumes mmraero above top interface = that of first model layer)
          taulw_aer(i,1)=optical_thickness(k_aer_lw_comb(1,i), &
               0.0, pint(1), tlay(1), 'k')  !optical depth above top interface, ajf 12/13/07
          

          do k=2, Nv_p1 ! first layer to last
             ! get optical thickness of layer k-1 
             taulw_aer(i,k) = &
                  optical_thickness(k_aer_lw_comb(k-1,i),&
                  pint(k-1), pint(k), tlay(k-1),'k')
             taulw_aer(i,k) = taulw_aer(i,k-1) + taulw_aer(i,k) ! add to layers above
          end do
    end do
    
    return
  end subroutine init_aer_lw_OD

!################################################################################

  !--------------------v  
  ! Combine_cia_ck_aer
  ! Author: Hamik Mukelyan, August 2007
  ! AJF, 3-25-08: Modified for CAM, also to include aerosols
  ! Description: Combines LW aerosol, CK and CIA optical depths.
  !
  ! In: tau_cia      -- (coll. pair, vertical level, wavenumber). CIA optical depths
  !     tau_ck       -- (wavenumber,gauss point,level). CK OD's.
  !     taulw_aer    -- (wavenumber,level). LW aerosol od's
  !     tau_tot      -- (wavenumber,gauss point, vertical level). Combined OD's
  !                     (Note: tau_tot dimensioning changed from Mukelyan's)
  !     ck           -- ck coefficient data structure
  !     cia          -- cia coef data struct
  !     comb_wn_grid -- grid of combined ck and cia wavenumbers
  !     nsm          -- max number gaussian points per wavenumber interval
  !     nsa          -- array of actual number of g points per wn interval
  !     i_cia_first  -- index of N2-N2 interaction, first index for *Titan*
  !
  ! Out: tau_tot
  subroutine combine_cia_ck_aer (tau_cia, tau_ck, taulw_aer, tau_tot, & 
                             cia, ck, comb_wn_grid, nsm, nsa, i_cia_first)
                             
    use radcnst, only: GasCia, GasCk, Nv_p1, N_cia, Nf_lw, arr_contains

!  Arguments

    real(r8), intent(inout) :: tau_tot(:,:,:)
    real(r8), intent(in)     :: tau_ck(:,:,:), tau_cia(:,:,:), taulw_aer(:,:)
    type(GasCk), intent(in)        :: ck
    type(GasCia), intent(in)       :: cia
    real(r8), intent(in)           :: comb_wn_grid(:)
    integer, intent(in)            :: nsm, nsa(:)
    integer                        :: i_cia_first

!  Local Variables

    integer                        :: i, j, k
    real(r8)                       :: sum_gas_taus
    
    ! initialize tau_tot array to vals in tau_ck, 0's for wn's where there
    ! is only cia and aerosol opacity
    
    do i=1, Nf_lw

       if(arr_contains(ck%Wn, comb_wn_grid(i))) then

          do j=1, Nv_p1
             do k=1, nsa(i)
                tau_tot(i,k,j)=tau_ck(where_arr(ck%Wn, comb_wn_grid(i)),k,j)
             end do
          end do

       else

          do j=1, Nv_p1
             do k=1, nsa(i)
                tau_tot(i, k, j) = 0.0
             end do
          end do

       end if
    end do

    ! now add the cia and aerosol contributions
    
    do i=1, Nf_lw
     
       do j=1, Nv_p1

          ! if this wn interval has cia aborption...
          if (arr_contains(cia%Wn, comb_wn_grid(i))) then

          ! sum the cia taus over the collision pair contributions. 
          ! summing only last four gases b/c others negligible on *titan*             
             sum_gas_taus = 0.0
             do k=i_cia_first, N_cia
                sum_gas_taus = sum_gas_taus + &
                     tau_cia(k, j, where_arr(cia%Wn, comb_wn_grid(i)))
             end do

             ! add the cia contribution to the taus at each gaussian point
             do k=1, nsa(i)
                tau_tot(i,k,j) = tau_tot(i,k,j) + sum_gas_taus
             end do

          end if
          
          ! Add the aerosol contribution
            tau_tot(i,1:nsa(i),j) = tau_tot(i,1:nsa(i),j) + taulw_aer(i,j)
          
       end do
    end do
    
  end subroutine combine_cia_ck_aer
  
!###############################################################################

  !--------------------v  
  ! optical_thickness
  ! Description: Calculates optical thickness using either a cia or ck 
  !              coefficient. The routine must know which coefficient to 
  !              expect since optical thickness has a different form for each.
  !              Units are cgs throughout (although pressure is passed in
  !              in units of atms).
  !
  ! In: p1        -- pressure at upper interface (p1 < p2), bars
  !     p2        -- pressure at lower interface
  !     cia_or_ck -- character which tells function to use either cia or ck
  !                  coefficient; 'a' for cia and 'k' for ck
  !     k         -- ck or cia coefficient (units = cm^2/g)
  !     T         -- temperature (K)
  !
  ! Out: The optical thickness.
  ! 
  ! Note: The ck 'k' comes in multiplied by (mass mixing ratio) / 
  !       [1.e20*(proton mass)*(molecular weight of gas)].  The original
  !        tabular values are in units of 1.e20*cm^2/molecule to follow
  !        Oxford convention.
  
  real(r8) function optical_thickness (k, p1, p2, T, cia_or_ck)
  
    use shr_const_mod, only:  SHR_CONST_BOLTZ, SHR_CONST_MWDAIR, SHR_CONST_G
     
    implicit none
    
! Arguments

    real(r8),  intent(in)      :: p1, p2, k, T
    character, intent(in)      :: cia_or_ck
    
! Local Variables
    real(r8)                   :: cfi
    real(r8), parameter :: n0 = 2.687e19  !Loschmidt's number in cm^-3
    real(r8), parameter :: protmas = 1.672e-24  !proton mass in g
    real(r8)            :: boltz     !Boltzmann's constant in cgs units
    real(r8)            :: grav_accl !in cgs units
    
    grav_accl=SHR_CONST_G*100.0  !convert to cm/s^2
        
    if (cia_or_ck == 'a') then
       boltz=SHR_CONST_BOLTZ*1.e7   !convert to cgs
       cfi = 1. / (n0 * boltz * SHR_CONST_MWDAIR * protmas * n0)
       optical_thickness = 0.5 * k * (p2*p2 - p1*p1) * cfi / T * & 
            1.e12 / grav_accl   ! factor of 1.e12 because pressure in bars --> dyne/cm^2
    else if (cia_or_ck == 'k') then
       optical_thickness = k * (p2 - p1) / grav_accl * 1.e6 ! factor of 1.e6 because pressure in bars
    else
       write(*,*) 'error in func optical_thickness: cia_or_ck must = a or k'
       stop
    end if
    
  end function optical_thickness
  
!###############################################################################

! interpck:
!   Interpolates the ck tables for a given temperature and pressure grid and
!   saves results in table format.
!
! author:
!   David Ye
!   AJF, 3-12-08: modified for use inside CAM code

subroutine interpck_cam(state)
  use shr_kind_mod,    only: r8 => shr_kind_r8

  use physics_types,   only: physics_state
  use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
  use phys_gmean,      only: gmean
  use radcnst,         only: N_gas_ck, CkPt, & !Number of lw-active gases, Ck table
                             Ng_ck, Nwn_ck, Ck, Ck_m, & ! Ck table on CAM P-T grid
                             CkComb, BilinearInterpolate
  Implicit None
  
!  Arguments

  type(physics_state), intent(in), dimension(begchunk:endchunk) :: state 
  
!  Local Variables

  Integer I, J, K, n, lchnk, ncol
  Integer Ngg, Nwn
  Integer Np           ! Number of pressure interfaces (pverp)
! ajf, 3-17-08: variables for computing global mean P-T profile:
  real(r8) :: te(pcols,begchunk:endchunk,1)  !temp storage of P,T flds at one vert-lev
  real(r8) :: te_glob(1)              ! global means of P,T at vertical level
  real(r8) :: pmean(pverp)            ! glob-mean pressure on interfaces
  real(r8) :: tmean(pverp)            !  "    "   temperature on interfaces
  real(r8) :: tm_lay(pver)            ! mean temperature in layers


!  Find global mean pressure, temperature fields at first time or restart step:

   do j=1,pverp
   
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
! interface pressure, Pascals, on interface j:
       te(:ncol,lchnk,1) = state(lchnk)%pint(:ncol,j)
    enddo
    
    call gmean(te, te_glob, 1)
    
    pmean(j)=te_glob(1) !mean pressure on interfaces
    
   enddo

! temperature, layer j
       
   do j=1,pver
   
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
! layer temperature, layer j:
       te(:ncol,lchnk,1) = state(lchnk)%t(:ncol,j)
    enddo
    
    call gmean(te, te_glob, 1)
    
    tm_lay(j)=te_glob(1) !mean temperature in layer j
    
   enddo
   
! Interpolate/Extrapolate global mean temperatures to interfaces

    j=1
     tmean(j)=tm_lay(j)
    do j=2,pver
     tmean(j)=0.5*(tm_lay(j-1)+tm_lay(j))
    enddo
    j=pverp
     tmean(j)=tm_lay(j-1) !mean surface T set equal to surf-layer T

!ajf, 3/12/08:
    Np=pverp
    do n=1,N_gas_ck
    do j=1,Np
     Ck(n)%P(j)=pmean(j)*1.0e-5  !convert pressure from Pa to bars
     Ck_m(n)%P(j)=Ck(n)%P(j)
    enddo 
    do j=1,Np
     Ck(n)%T(j)=tmean(j)
     Ck_m(n)%T(j)=Ck(n)%T(j)
    enddo
    enddo

!!!  print*,'FROM INTERPCK:'
!!!  print*,'Pressure grid for Ck values: '
!!!  print*, Ck(1)%P

!!!  print*,'Temperature grid for Ck values: '
!!!  print*, Ck(1)%T

!!!  Write(*,*) "Input ck file in form (gas)-(collisionpartner).ck:"
!!!  Read(*,*) infil
!!!  Write(*,*) "Output ck file in form (gas)-(collisionpartner)-(version name).ck:"
!!!  Read(*,*) outfil
  

!  Interpolate grand table CkPt onto CAM P-T grid table Ck for each gas

  do n=1,N_gas_ck

  Nwn = Nwn_ck(n)
  Ngg =  Ng_ck(n)

! Interpolate at pressure and temperature coordinates specified
        
  Do I = 1, Ngg
  Do J = 1, Np
  Do K = 1, Nwn
  
!!!    if( k == Nwn .and. i == 1 .and. j == 1) then
!!!     if(masterproc) then
!!!      print*,'INTERPCK_CAM:'
!!!      print*,'n= ',n,' k= ',k,' i= ',i,' j= ',j
!!!      print*,' Nwn= ',Nwn,' Ngg= ',Ngg
!!!      print*,' Ck(n)%P,T(j)= ',Ck(n)%P(j),Ck(n)%T(j)
!!!      print*,' CkPt(n)%T(8,9)= ',CkPt(n)%T(8),CkPt(n)%T(9)
!!!      print*,' CkPt(n)%P(4,5)= ',CkPt(n)%P(4),CkPt(n)%P(5)
!!!      print*,' CkPt(n)%Val(1,4,8,Nwn)= ',CkPt(n)%Val(1,4,8,Nwn)
!!!      print*,' CkPt(n)%Val(1,5,8,Nwn)= ',CkPt(n)%Val(1,5,8,Nwn)
!!!      print*,' CkPt(n)%Val(1,4,9,Nwn)= ',CkPt(n)%Val(1,4,9,Nwn)
!!!      print*,' CkPt(n)%Val(1,5,9,Nwn)= ',CkPt(n)%Val(1,5,9,Nwn)
!!!     endif
!!!    endif     
          
     Ck(n)%Val(I, J, K) = BilinearInterpolate(log(CkPt(n)%P), CkPt(n)%T, &
           CkPt(n)%Val(I, :, :, K), log(Ck(n)%P(J)), Ck(n)%T(J))
           
!!!    if( k == Nwn .and. i == 1 .and. j == 1) then
!!!     if(masterproc) then
!!!      print*,'Ck(n)%Val(I, J, K)= ',Ck(n)%Val(I, J, K)
!!!    endif
!!!    endif           
    
  End Do
  End Do
  End Do
  
  
  
  enddo !end loop over the n_gas_ck gases
  
! Finish initialization of CkComb to retain the P, T arrays of CAM grid

  CkComb%P = Ck(1)%P
  CkComb%T = Ck(1)%T

end subroutine interpck_cam

!###############################################################################

! interpcia_cam:
!   Creates table of Cia values interpolated to model reference  
!   pressure/temperature profile (which can be periodically updated).
!
! author:
!   David Ye
!   AJF, 3-21-08: Extensive modifications for CAM

Subroutine interpcia_cam
  
  use shr_kind_mod, only: r8 => shr_kind_r8
  
  use radcnst,  only: N_cia, CkComb, BilinearInterpolate, &
                      Cia, CiaComb, Cia_v, Nv_p1

  Implicit None
  
!  Local Variables

  Integer            :: Nwn
  Integer            :: I, J, K


! Pass the "current" temperature profile into CiaComb
! (Wavenumber grid already established in radcnst_initialize)

  CiaComb%T = CkComb%T
  Cia_v%T=CiaComb%T

! Bilinearly interpolate to return values at the required temperature and
! wavenumber.

  Nwn=Size(CiaComb%Wn)
  
  Do K = 1, N_cia
  Do J = 1, Nwn
  Do I = 1, Nv_p1
  
     CiaComb%Val(K, I, J) = BilinearInterpolate(Log(Cia%T), Cia%Wn, &
                            Cia%Val(K, :, :), Log(CiaComb%T(I)), CiaComb%Wn(J))
                            
     ! The Cia vals have been natural log'd, so exponentiate CiaComb                       
     CiaComb%Val(K, I, J) = exp(CiaComb%Val(K,I,J))
                                 
  End Do
  End Do
  End Do

End subroutine interpcia_cam
  
!###############################################################################
  
  !--------------------v  
  ! where_arr
  ! Description: Returns index of first instance of 'num' in array argument.
  !              -1 if not found.
  !
  ! In: arr  -- array
  !     num  -- ...
  !
  ! Out: See description.
  integer function where_arr (arr, num)

!Arguments

    real(r8), intent(in) :: arr(:), num
    
!Local

    integer :: i
    
    where_arr = -1
    do i = 1, size(arr)
       if (arr(i) == num) then
          where_arr = i
          return
       end if
    end do
  end function where_arr
  
!###############################################################################

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END SUBROUTINE INDEXX
  

End Module taulw_cam