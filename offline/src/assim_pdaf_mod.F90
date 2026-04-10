!>  Module for assimilation variables
!!
!! This module provides variables needed by the user code for the
!! assimilation. For simplicity, all assimilation-related variables
!! are stored here, even if they are only used in the main program
!! for the filter initialization.
!! When using the template init_pdaf_parse, most variables can be
!! specified as a command line argument.
!!
!! Implementation for the 2D online example with two fields.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
module assim_pdaf_mod

  implicit none
  save

!+++ Specific part for 2D tutorial model

  real :: coords_l(2)                 !< Coordinates of local analysis domain (only size is case-specific)
  integer :: type_ens_init=1          !< Type of ensemble init:
                                      !< (1) read model files; use mean state from these files
                                      !< (2) read model files; use mean state from covariance file 
                                      !< (3) generate ensemble from covariance matrix
  character(len=100) :: file_covar    !< Path and name of covariance matrix file

!+++ End of specific part

! Do not remove the following line, it is relevant for OpenMP parallelization
!$OMP THREADPRIVATE(coords_l)


! -----------------------------------------------------------------
! --- Below are the generic variables used for configuring PDAF ---
! --- Default values are set here, and deviation in init_pdaf   ---

! Settings for state vector size
  integer :: dim_state          !< Global model state dimension
  integer :: dim_state_p        !< Model state dimension for Process-local domain

! Settings for time stepping - available as command line options
  logical :: model_error        !< Control application of model error
  real    :: model_err_amp      !< Amplitude for model error
  integer :: step_offline       !< Model step to process in offline mode

! Settings for observations - available as command line options
  integer :: delt_obs=1         !< time step interval between assimilation steps
  integer :: observe_ens=0      !< (0) apply H also to ensemble mean; (1) apply H only to ensemble states
  integer :: type_obs_init=0    !< init obs. (0) before or (1) after call to prepostsstep
  logical :: twin_experiment=.false. !< Whether to run an twin experiment with synthetic observations
  logical :: do_omi_obsstats=.false. !< Whether to let OMI compute observation statistics

! General control of PDAF - available as command line options
  integer :: screen=2           !< Control verbosity of PDAF
                                !< * (0) no outputs
                                !< * (1) progress info
                                !< * (2) add timings
                                !< * (3) debugging output
  integer :: dim_ens            !< Size of ensemble
  integer :: filtertype=6       !< Select filter algorithm:
                                !<   * SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                                !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                                !<   LKNETF (11), PF (12), EnSRF/EAKF (13)
                                !< GENOBS (100), 3DVAR (200)
  integer :: subtype=0          !< Subtype of filter algorithm
                                !<   * SEIK:
                                !<     (0) ensemble forecast; new formulation
                                !<     (1) ensemble forecast; old formulation
                                !<     (2) SEIK with ensemble transformation
                                !<     (10) fixed error space basis
                                !<     (11) fixed state covariance matrix
                                !<   * LSEIK:
                                !<     (0) ensemble forecast;
                                !<     (2) LSEIK with ensemble transformation
                                !<     (10) fixed error space basis
                                !<     (11) fixed state covariance matrix
                                !<   * ETKF:
                                !<     (0) ETKF using T-matrix like SEIK
                                !<     (1) ETKF following Hunt et al. (2007)
                                !<       There are no fixed basis/covariance cases, as
                                !<       these are equivalent to SEIK subtypes 2/3
                                !<   * LETKF:
                                !<     (0) LETKF using T-matrix like SEIK
                                !<     (1) LETKF following Hunt et al. (2007)
                                !<       There are no fixed basis/covariance cases, as
                                !<       these are equivalent to LSEIK subtypes 2/3
                                !<   * EnKF:
                                !<     (0) analysis for large observation dimension
                                !<     (1) analysis for small observation dimension
                                !<   * LEnKF:
                                !<     (0) standard analysis
                                !<   * ESTKF:
                                !<     (0) standard ESTKF 
                                !<       There are no fixed basis/covariance cases, as
                                !<       these are equivalent to SEIK subtypes 2/3
                                !<   * LESTKF:
                                !<     (0) standard LESTKF 
                                !<       There are no fixed basis/covariance cases, as
                                !<       these are equivalent to LSEIK subtypes 2/3
                                !<   * NETF:
                                !<     (0) standard NETF 
                                !<   * LNETF:
                                !<     (0) standard LNETF 
                                !<   * LKNETF:
                                !<     (0) HNK: 2-step LKNETF with NETF before LETKF
                                !<     (1) HKN: 2-step LKNETF with LETKF before NETF
                                !<     (2) HSync: LKNETF synchronous
                                !<   * PF:
                                !<     (0) standard PF 
                                !<   * ENSRF/EAKF:
                                !<     (0) ENSRF with serial observation processing
                                !<     (1) EAKF with local least square regression
                                !<   * 3D-Var:
                                !<     (0) parameterized 3D-Var
                                !<     (1) 3D Ensemble Var using LESTKF for ensemble update
                                !<     (2) 3D Ensemble Var using ESTKF for ensemble update
                                !<     (3) hybrid 3D-Var using LESTKF for ensemble update
                                !<     (4) hybrid 3D-Var using ESTKF for ensemble update
  integer :: type_iau=0         !< Type of incremental updating:
                                !<     (0) no IAU
                                !<     (1) constant IAU weight
                                !<     (2) linear increase/decrease with maimum in middle of period
                                !<     (3) Null IAU: initialize increments arrays, but do not add increment
  integer :: steps_iau=1        !< Number of time steps over which IAU is applied
  integer :: dim_lag=0          !< Number of time instances for smoother

! Filter settings - available as command line options
!    ! General
  integer :: type_forget=0      !< Type of forgetting factor
                                !<  SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                                !<   (0) fixed
                                !<   (1) global adaptive
                                !<   (2) local adaptive for LSEIK/LETKF/LESTKF
                                !<  NETF/LNETF/PF
                                !<   (0) apply inflation on forecast ensemble
                                !<   (2) apply inflation on analysis ensemble
  real    :: forget=1.0         !< Forgetting factor for inflation at filter analysis
  integer :: dim_bias=0         !< dimension of bias vector
!    ! All localized filters
  real    :: cradius=0.0        !< Cut-off radius for local observation domain
  integer :: locweight=2        !< * Type of localizing weighting of observations
                                !<   (0) constant weight of 1
                                !<   (1) exponentially decreasing with SRADIUS
                                !<   (2) use 5th-order polynomial
                                !<   (3) regulated localization of R with mean error variance
                                !<   (4) regulated localization of R with single-point error variance
  real    :: sradius=0.0        !< Support radius for 5th order polynomial
                                !<   or radius for 1/e for exponential weighting
!    ! ENKF
  integer :: rank_ana_enkf=0    !< Rank to be considered for inversion of HPH in analysis of EnKF
                                !<  (0) for analysis w/o eigendecomposition
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF/NETF/LNETF/LKNETF
  integer :: type_trans=0       !< Type of ensemble transformation 
                                !< * SEIK/LSEIK: 
                                !< (0) use deterministic omega
                                !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                                !< (2) use product of (0) with random orthonormal matrix with
                                !<     eigenvector (1,...,1)^T 
                                !< * ETKF/LETKF with subtype=4: 
                                !< (0) use deterministic symmetric transformation
                                !< (2) use product of (0) with random orthonormal matrix with
                                !<     eigenvector (1,...,1)^T 
                                !< * ESTKF/LESTKF:
                                !< (0) use deterministic omega
                                !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                                !< (2) use product of (0) with random orthonormal matrix with
                                !<     eigenvector (1,...,1)^T
                                !< * NETF/LNETF:
                                !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                                !< (1) use identity transformation
                                !< * LKNETF:
                                !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                                !< (1) use identity transformation
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  integer :: type_sqrt=0        !< * Type of the transform matrix square-root 
                                !<   (0) symmetric square root
                                !<   (1) Cholesky decomposition
!    ! NETF/LNETF/PF
  integer :: type_winf=0        !< Set weights inflation: 
                                !<   (0) no weights inflation
                                !<   (1) use N_eff/N>limit_winf
  real    :: limit_winf=0.5     !< Limit for weights inflation: N_eff/N>limit_winf
!    ! hybrid LKNETF
  integer :: type_hyb=3         !< * Type of hybrid weight:
                                !<   (0) use fixed hybrid weight hyb_gamma
                                !<   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                                !<   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
                                !<   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                                !<   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  real    :: hyb_gamma=0.8      !< Hybrid filter weight for state (1.0: LETKF, 0.0 LNETF)
  real    :: hyb_kappa          !< Hybrid norm for using skewness and kurtosis (recommended default: dim_ens)
!    ! Particle filter
  integer :: pf_res_type=1      !< * Resampling type for PF
                                !<   (1) probabilistic resampling
                                !<   (2) stochastic universal resampling
                                !<   (3) residual resampling        
  integer :: pf_noise_type=0    !< * Resampling type for PF
                                !<   (0) no perturbations, (1) constant stddev, 
                                !<   (2) amplitude of stddev relative of ensemble variance
  real :: pf_noise_amp          !< Noise amplitude (>=0.0, only used if pf_noise_type>0)

!    ! 3D-Var
  integer :: type_opt=1         !< * Type of minimizer for 3DVar
                                !<   (1) LBFGS (default)
                                !<   (2) CG+
                                !<   (3) plain CG
                                !<   (12) CG+ parallelized
                                !<   (13) plain CG parallelized
  integer :: dim_cvec=0         !< Size of control vector (parameterized part; for subtypes 0,1)
  integer :: dim_cvec_ens=0     !< Size of control vector (ensemble part; for subtypes 1,2)
  integer :: mcols_cvec_ens=1   !< Multiplication factor for number of columns for ensemble control vector
  real :: beta_3dvar = 0.5      !< Hybrid weight for hybrid 3D-Var
  integer :: solver_iparam1=2   !< Solver specific parameter
                                !<  LBFGS: parameter m (default=5)
                                !<       Number of corrections used in limited memory matrix; 3<=m<=20
                                !<  CG+: parameter method (default=2)
                                !<       (1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere
                                !<  CG: maximum number of iterations (default=200)
  integer :: solver_iparam2=1   !< Solver specific parameter
                                !<  LBFGS: - not used - 
                                !<  CG+: parameter irest (default=1)
                                !<       (0) no restarts; (n>0) restart every n steps
                                !<  CG: - not used -
  real :: solver_rparam1=1.0e-6 !< Solver specific parameter
                                !<  LBFGS: limit for stopping iterations 'pgtol' (default=1.0e-5)
                                !<  CG+: convergence parameter 'eps' (default=1.0e-5)
                                !<  CG: conpergence parameter 'eps' (default=1.0e-6)
  real :: solver_rparam2=1.0e+7 !< Solver specific parameter
                                !<  LBFGS: tolerance in termination test 'factr' (default=1.0e+7) 
                                !<  CG+: - not used -
                                !<  CG: - not used -

!    ! Other variables - _NOT_ available as command line options!
  real    :: time               !< model time
  real, allocatable :: coords_p(:,:)    !< Coordinates of process-local state vector entries
                                        !< needed to intialize localization for LEnKF/ENSRF

end module assim_pdaf_mod
