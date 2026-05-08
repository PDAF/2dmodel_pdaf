!> Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the offline mode of PDAF.
!!
!! This routine is generic. However, it assumes a constant observation
!! error for each observation type (rms_obs_OBSTYPE).
!!
!! __Revision history:__
!! * 2008-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine init_pdaf_offline()

  use PDAF, &                          ! PDAF
       only: PDAF3_init, PDAF_set_iparam, PDAF_set_rparam, &
       PDAFomi_set_domain_limits, PDAF_abort, &
       PDAF_DA_NETF, PDAF_DA_LNETF, PDAF_DA_PF, PDAF_DA_LKNETF
  use parallel_pdaf_mod, &             ! Parallelization variables
       only: myproc_ens, myproc_assim
  use assim_pdaf_mod, &                ! Variables for assimilation
       only: screen, dim_state_p, dim_state, dim_ens, filtertype, subtype, delt_obs, &
       step_offline, type_forget, forget, locweight, cradius, sradius, coords_p, &
       type_trans, type_sqrt, observe_ens, type_obs_init, &
       type_winf, limit_winf, pf_res_type, pf_noise_type, pf_noise_amp, &
       type_hyb, hyb_gamma, hyb_kappa, type_obs_init, type_ens_init, file_covar
  use statevector_pdaf_mod, &          ! State vector variables and init routine
       only: setup_statevector, n_fields

  ! Specific for model
  use model_pdaf_mod, &                ! Model variables
       only: n_dim !nx_p, ny, coords_x_p, coords_y_p

  ! Specific for each observation type
  use obs_OBSTYPE_pdafomi, &           ! Variables for observation OBSTYPE
       only: assim_OBSTYPE, rms_obs_OBSTYPE

  implicit none

! *** Local variables ***
  integer :: pdaf_param_i(2)           ! Integer parameter array for filter
  real    :: pdaf_param_r(1)           ! Real parameter array for filter
  integer :: status_pdaf               ! PDAF status flag
  real    :: lim_coords(2,2)           ! limiting coordinates of process sub-domain

! *** External subroutines ***
  external :: init_ens_cb_pdaf         ! Ensemble initialization
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  if (myproc_ens == 0) then
     write (*,'(/a,1x,a)') 'model-PDAF', 'INITIALIZE PDAF - OFFLINE MODE'
  end if


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF3_init     ***
! **********************************************************

! *** Ensemble settings ***
  dim_ens = 6             ! Size of ensemble (arbitrary value; the actual value is set at runtime)

! *** Model time step to process
  step_offline = 1        ! This determines which observation are read (can be set on command line)

! *** Options for DA method

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++ For available options see ASSIM_PDAF_MOD and INIT_PDAF_PARSE +++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ************************************************************
! *** Settings used in call-back routines                  ***
! *** Here we specify value that deviate from the defaults ***
! ************************************************************

! *** Localization settings - usually set at run-time
  cradius = 5.0           ! Cut-off radius in local filters (in units of model coordinate)

!+++ Specific variables for observations - see defaults in each observation module

! *** Which observation type to assimilate
  assim_OBSTYPE = .true.

!+++ End of specific variables for observations


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options - optional, but useful
! *** One could also use a namelist file

  call init_pdaf_parse()


! ***************************
! *** Define state vector ***
! ***************************

  call setup_statevector(dim_state, dim_state_p, screen)


! ********************************************************
! *** Call PDAF initialization routine - all processes ***
! ***                                                  ***
! *** For all filters, PDAF3_init is first called      ***
! *** specifying only the required parameters.         ***
! *** Further settings are done afterwards using       ***
! *** calls to PDAF_set_iparam & PDAF_set_rparam.      ***
! ********************************************************

  ! *** Here we specify only the required integer and real parameters
  ! *** Other parameters are set using calls to PDAF_set_iparam/PDAF_set_rparam
  pdaf_param_i(1) = dim_state_p ! State dimension
  pdaf_param_i(2) = dim_ens     ! Size of ensemble
  pdaf_param_r(1) = forget      ! Forgetting factor

  call PDAF3_init(filtertype, subtype, step_offline, &
       pdaf_param_i, 2,&
       pdaf_param_r, 1, &
       init_ens_cb_pdaf, screen, status_pdaf)

  ! *** Additional parameter specifications ***
  ! *** -- These are all optional --        ***

  ! Generic settings
  call PDAF_set_iparam(5, type_forget, status_pdaf)      ! Type of forgetting factor
  call PDAF_set_iparam(6, type_trans, status_pdaf)       ! Type of ensemble transformation
  call PDAF_set_iparam(7, type_sqrt, status_pdaf)        ! Type of transform square-root (SEIK-sub4/ESTKF)
  call PDAF_set_iparam(8, observe_ens, status_pdaf)      ! Whether to apply observation operator to ensemble mean
  call PDAF_set_iparam(9, type_obs_init, status_pdaf)    ! Initialize observation before or after call to prepoststep

  ! Settings for NETF and LNETF
  if (filtertype==PDAF_DA_NETF .or. filtertype==PDAF_DA_LNETF) then
     call PDAF_set_iparam(4, pf_noise_type, status_pdaf) ! Perturbation type
     call PDAF_set_iparam(7, type_winf, status_pdaf)     ! Type of weights inflation
     call PDAF_set_rparam(2, limit_winf, status_pdaf)    ! Limit for weights inflation
     call PDAF_set_rparam(3, pf_noise_amp, status_pdaf)  ! Noise amplitude
  end if

  ! Settings for particle filter PF
  if (filtertype==PDAF_DA_PF) then
     call PDAF_set_iparam(4, pf_noise_type, status_pdaf) ! Perturbation type
     call PDAF_set_iparam(6, pf_res_type, status_pdaf)   ! Resampling type
     call PDAF_set_iparam(7, type_winf, status_pdaf)     ! Type of weights inflation
     call PDAF_set_rparam(2, limit_winf, status_pdaf)    ! Limit for weights inflation
     call PDAF_set_rparam(3, pf_noise_amp, status_pdaf)  ! Noise amplitude
  end if

  ! Settings for hybrid filter LKNETF
  if (filtertype==PDAF_DA_LKNETF) then
     call PDAF_set_iparam(4, type_hyb, status_pdaf)      ! Choice of hybrid rule
     call PDAF_set_rparam(2, hyb_gamma, status_pdaf)     ! Hybrid filter weight for state
     call PDAF_set_rparam(3, hyb_kappa, status_pdaf)     ! Normalization factor for hybrid weight 
  end if

! *** Check whether initialization of PDAF was successful ***
  if (status_pdaf /= 0) then
     write (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, ' in initialization of PDAF - stopping! (Process ', myproc_ens,')'
     call PDAF_abort(1)
  end if


! ***************************************************
! *** Set coordinates of elements in state vector ***
! *** (used for localization in EnKF/ENSRF)       ***
! ***************************************************

  ! For localization in EnKF and EnSRF/EAKF, PDAFomi_set_localize_covar
  ! is called in the observation modules. This routine requires a 
  ! coordinate array corresponding to the state vector.

  allocate(coords_p(n_dim, dim_state_p))

  ! +++ This can be filled using coords_x_p and coords_y_p
!  coords_p(i,j) = ?


! *************************************************************************
! *** Set sub-domain coordinate limits                                  ***
! *** Used when the PDAF-OMI option thisobs%use_global_obs is set to 0  ***
! *************************************************************************
  
! +++ This is optional

!   lim_coords(1,1) = ?     ! West
!   lim_coords(1,2) = ?     ! East
!   lim_coords(2,1) = ?     ! North
!   lim_coords(2,2) = ?     ! South

  call PDAFomi_set_domain_limits(lim_coords)

end subroutine init_pdaf_offline
