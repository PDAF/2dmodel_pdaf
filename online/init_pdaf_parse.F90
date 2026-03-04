!>  Parse command line options for PDAF
!!
!! This routine calls the command line parser to initialize
!! variables for the data assimilation with PDAF.
!!
!! Using the parser is optional and shows one possibility
!! to modify the variables of the compiled program. An 
!! alternative to this might be Fortran namelist files.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_pdaf_parse()

  use parser, &                ! Parser function
       only: parse
  use assimilation_pdaf_mod, & ! Variables for assimilation
       only: screen, filtertype, subtype, dim_ens, delt_obs, &
       step_offline, twin_experiment, &
       model_error, model_err_amp, type_forget, forget, &
       type_iau, steps_iau, rank_ana_enkf, &
       locweight, cradius, sradius, &
       type_trans, type_sqrt, dim_lag, type_hyb, &
       hyb_gamma, hyb_kappa, type_winf, limit_winf, &
       pf_res_type, pf_noise_type, pf_noise_amp, &
       observe_ens, type_obs_init, do_omi_obsstats, &
       type_ens_init, file_covar
  use io_pdaf_mod, &           ! File input/output control
       only: write_state, write_ens, write_var

  ! Specific for 2D tutorial model
  use obs_A_pdafomi, &         ! Variables for observation type A
       only: assim_A, rms_obs_A, file_obs_A
  use obs_B_pdafomi, &         ! Variables for observation type B
       only: assim_B, rms_obs_B, file_obs_B

  implicit none

! *** Local variables ***
  character(len=32) :: handle  ! handle for command line parser


! **********************************
! *** Parse command line options ***
! **********************************

!+++ Specific part for 2D tutorial model

  ! Observation settings - particular for the implemented observation modules
  handle = 'assim_A'                 ! Whether to assimilation observation type A
  call parse(handle, assim_A)
  handle = 'assim_B'                 ! Whether to assimilation observation type B
  call parse(handle, assim_B)
  handle = 'rms_obs_A'               ! Assumed uniform RMS error of the observations type A
  call parse(handle, rms_obs_A)
  handle = 'rms_obs_B'               ! Assumed uniform RMS error of the observations type B
  call parse(handle, rms_obs_B)
  handle = 'file_obs_A'              ! Path and name of observation file type A
  call parse(handle, file_obs_A)
  handle = 'file_obs_B'              ! Path and name of observation file type B
  call parse(handle, file_obs_B)

  ! Settings for ensemble initialization
  handle = 'type_ens_init'           ! Type of ensemble initialization
  call parse(handle, type_ens_init)
  handle = 'file_covar'              ! Path and name of covariance matrix file
  call parse(handle, file_covar)

  ! Setting controlling file output
  handle = 'write_state'             ! Whether to write ensemble mean fields
  call parse(handle, write_state)
  handle = 'write_ens'               ! Whether to write ensemble files
  call parse(handle, write_ens)
  handle = 'write_var'               ! Whether to write ensemble variance files
  call parse(handle, write_var)

!+++ End of specific part

!------------------------------------------------------------------------------
! The remaining parse commands should be generic; usually no change necessary

  ! Observation settings
  handle = 'delt_obs'                ! Time step interval between filter analyses
  call parse(handle, delt_obs)
  handle = 'observe_ens'             ! (0) apply H also to ensemble mean; (1) apply H only to ensemble states
  call parse(handle, observe_ens)
  handle = 'type_obs_init'           ! init obs. (0) before or (1) after call to prepostsstep
  call parse(handle, type_obs_init)
  handle = 'do_omi_obsstats'         ! Whether to let PDAF-OMI compute observation statistics
  call parse(handle, do_omi_obsstats)
  handle = 'twin_experiment'         ! T: perform twin experiment with synthetic observations
  call parse(handle, twin_experiment)


  ! Model step to process in offline mode
  handle = 'step'                    ! Time step interval between filter analyses
  call parse(handle, step_offline)

  ! Settings for model and time stepping
  handle = 'model_error'             ! Control application of model error
  call parse(handle, model_error)
  handle = 'model_err_amp'           ! Amplitude of model error
  call parse(handle, model_err_amp)

  ! General settings for PDAF
  handle = 'screen'                  ! set verbosity of PDAF
  call parse(handle, screen)
  handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
  call parse(handle, dim_ens)
  handle = 'filtertype'              ! Choose filter algorithm
  call parse(handle, filtertype)
  handle = 'subtype'                 ! Set subtype of filter
  call parse(handle, subtype)

  ! Control IAU
  handle = 'type_iau'                ! Set whether to use incremental updating
  call parse(handle, type_iau)
  handle = 'steps_iau'               ! Number of time steps over which IAU is applied
  call parse(handle, steps_iau)

  ! Settings for smoother
  handle = 'dim_lag'                 ! Size of lag in smoother
  call parse(handle, dim_lag)

  ! Filter-specific settings
  handle = 'forget'                  ! Set forgetting factor
  call parse(handle,forget)
  handle = 'type_forget'             ! Set type of forgetting factor
  call parse(handle, type_forget)
  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  call parse(handle, type_trans)
  handle = 'type_sqrt'               ! Set type of transformation square-root (SEIK-sub4, ESTKF)
  call parse(handle, type_sqrt)
  handle = 'rank_ana_enkf'           ! Set rank for pseudo inverse in EnKF
  call parse(handle, rank_ana_enkf)

  ! Settings for localization in LSEIK/LETKF
  handle = 'cradius'                 ! Set cut-off radius in grid points for observation domain
  call parse(handle, cradius)
  handle = 'locweight'               ! Set type of localizating weighting
  call parse(handle, locweight)
  sradius = cradius                  ! By default use cradius as support radius
  handle = 'sradius'                 ! Set support radius for 5th-order polynomial
                                     !  or radius for 1/e in exponential weighting
  call parse(handle, sradius)

  ! Settings for nonlinear filters
  handle = 'pf_res_type'             ! Resampling type for particle filter
  call parse(handle, pf_res_type)        
  handle = 'pf_noise_type'           ! Type of perturbing noise in PF
  call parse(handle, pf_noise_type)        
  handle = 'pf_noise_amp'            ! Amplitude of perturbing noise in PF
  call parse(handle, pf_noise_amp)        
  handle = 'type_winf'               ! Set type of weights inflation in NETF/LNETF
  call parse(handle, type_winf)
  handle = 'limit_winf'              ! Set limit for weights inflation
  call parse(handle, limit_winf)

  ! Hybrid weights for LKNETF
  handle = 'type_hyb'                ! Set type of hybrid weight
  call parse(handle, type_hyb)
  handle = 'hyb_gamma'               ! Set hybrid filter weight for state (1.0 LETKF, 0.0 LNETF)
  call parse(handle, hyb_gamma)
  handle = 'hyb_kappa'               ! Set hybrid norm (>1.0)
  call parse(handle, hyb_kappa)

end subroutine init_pdaf_parse
