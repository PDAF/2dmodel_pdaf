!>  Initialize communicators for PDAF
!!
!! Parallelization routine for a model with attached PDAF. The subroutine is
!! called in the main program subsequently to the initialization of MPI. It
!! initializes MPI communicators for the model tasks, filter task and the
!! coupling between model and filter tasks. In addition some other variables 
!! for the parallelization are initialized.
!! The communicators and variables are handed over to PDAF in the call to 
!! PDAF_set_parallel toward the end of this routine.
!!
!! 3 Communicators are generated:
!! * _COMM_filter_: Communicator in which the filter itself operates
!! * _COMM_model_: Communicators for parallel model forecasts
!! * _COMM_couple_: Communicator for coupling between models and filter
!!
!! In addition there is the main communicator
!! * _COMM_ensemble_: The main communicator in which PDAF operates
!! COMM_ensemble is set to the communicator in which all model integration 
!! are computed. Typically, this is MPI_COMM_WORLD, but it can be defined
!! differently if the model only operators on a subset to MPI_COMM_WORLD.
!! This happens, e.g. if some processes are separated to operate an
!! I/O server or a model coupler for coupled model systems.
!!
!! Other variables that have to be initialized are:
!! * _filterpe_ - Logical: Does the Process execute the filter?
!! * _task_id_ - Integer: Index identifying the model task
!! * _my_ensemble_ - Integer: The index of the Process's model task
!! * _local_npes_model_ - Integer array holding numbers of Processs per model task
!!
!! For COMM_filter and COMM_model also the size of the communicators
!! (npes_filter and npes_model) and the rank of each process  (mype_filter,
!! mype_model) are initialized. These variables can be used in the model part 
!! of the program, but are not handed over to PDAF.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_parallel_pdaf(screen, model_comm, model_comm_rank, model_comm_size)

  use mpi                         ! MPI
  use PDAF, &                     ! PDAF routines
       only: PDAF3_set_parallel
  use parallel_pdaf_mod, &        ! PDAF parallelization variables
       only: mype_ens, npes_ens, COMM_ensemble, &
       mype_model, npes_model, COMM_model, &
       mype_couple, npes_couple, COMM_couple, &
       mype_filter, npes_filter, COMM_filter, filterpe, &
       n_modeltasks, task_id, local_npes_model, MPIerr
  use parser, &                   ! Command line parser
       only: parse

  implicit none

! *** Arguments ***
  integer, intent(in)    :: screen           !< Whether screen information is shown

  ! Model variables for parallelization (one can keep these generic names)
  integer, intent(inout) :: model_comm       !< Model MPI communicator for model tasks
  integer, intent(inout) :: model_comm_size  !< Number of Processs in model_comm
  integer, intent(inout) :: model_comm_rank  !< Process rank in model_comm

! *** Local variables ***
  integer :: i, j                     ! Counters
  integer :: pe_index                 ! Index of Process
  integer :: my_color, color_couple   ! Variables for communicator-splitting 
  integer :: flag                     ! Status flag
  integer :: dim_ens                  ! Ensemble size / number of model tasks
  character(len=32) :: handle         ! Handle for command line parser
  logical :: online_coupling          ! Whether to ru nonline coupled DA


  ! Specify online of offline coupling
  online_coupling = .true.

  ! *** Define ensemble communicator.                   ***
  ! *** This is the communicator in which PDAF operates ***

  COMM_ensemble = model_comm


!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the online mode of PDAF.
!!
!! This routine is generic. However, it assumes a constant observation
!! error (rms_obs). Further, with parallelization the local state
!! dimension dim_state_p is used.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_pdaf()

  use PDAF, &                          ! PDAF
       only: PDAF3_init, PDAF_set_iparam, PDAF_init_forecast, &
       PDAFomi_set_domain_limits, PDAF_iau_init
  use parallel_pdaf_mod, &             ! Parallelization variables
       only: mype_ens, mype_filter, n_modeltasks, abort_parallel
  use assimilation_pdaf_mod, &         ! Variables for assimilation
       only: dim_state_p, dim_state, dim_ens, &
       screen, filtertype, subtype, delt_obs, step_offline, &
       type_iau, steps_iau, type_forget, forget, &
       locweight, cradius, sradius, coords_p, &
       type_obs_init, type_ens_init, file_covar
  use statevector_pdaf_mod, &          ! State vector variables and init routine
       only: setup_statevector, n_fields
  use io_pdaf_mod, &                   ! File input/output control
       only: write_state, write_ens, write_var

  ! Specific for 2D tutorial model
  use model_pdaf_mod, &                ! Model variables
       only: nx_p, ny, n_dim, coords_x_p, coords_y_p
  use obs_A_pdafomi, &                 ! Variables for observation type A
       only: assim_A, rms_obs_A, file_obs_A
  use obs_B_pdafomi, &                 ! Variables for observation type B
       only: assim_B, rms_obs_B, file_obs_B

  implicit none

! *** Local variables ***
  integer :: i, j, k, s, off_nx        ! Counters
  integer :: pdaf_param_i(2)           ! Integer parameter array for filter
  real    :: pdaf_param_r(1)           ! Real parameter array for filter
  integer :: status_pdaf               ! PDAF status flag
  real    :: lim_coords(2,2)           ! limiting coordinates of process sub-domain

! *** External subroutines ***
  external :: init_ens_pdaf            ! Ensemble initialization
  external :: next_observation_pdaf, & ! Provide time step of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_pdaf                ! User supplied pre/poststep routine
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  if (mype_ens == 0) then
     write (*,'(/a,1x,a)') 'model-PDAF', 'INITIALIZE PDAF - ONLINE MODE'
  end if


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen = 2              ! Write screen output (1) for output, (2) add timings

! *** Ensemble settings ***
  dim_ens = n_modeltasks  ! Size of ensemble
                          !   We use n_modeltasks here, initialized in init_parallel_pdaf

! *** Options for DA method

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++ For available options see ASSIMILATION_PDAF_MOD +++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++

  filtertype = 6     ! Type of filter
  subtype = 0        ! Subtype of filter

  forget  = 1.0      ! Forgetting factor value for inflation
  type_forget = 0    ! Type of forgetting factor


! **********************************************************
! ***   Settings used in call-back routines              ***
! **********************************************************

! *** Forecast length ***
  delt_obs = 2            ! Number of time steps between analysis steps

! *** IO options
  write_state = .true.    ! Write ensemble mean fields
  write_ens = .true.      ! Write all ensemble states
  write_var = .false.     ! Write ensemble variance fields

! *** Specifications for ensemble initialization
  type_ens_init = 1       ! Type of ensemble initialization
  file_covar = '../generate_covar/covar.nc'  ! Path and name of covariance matrix file
                                             ! (generated by program generate_covar)

! *** Whether to initialize observations before prepoststep
  type_obs_init = 0       ! (0) before, (1) after

! *** Incremental updating (IAU)
  type_iau = 0            ! Type of incremental updating
  steps_iau = 1           ! Number of time steps over which IAU is applied

! *** Localization settings
  locweight = 2           ! Type of localizating weighting
  cradius = 5.0           ! Cut-off radius in local filters (in units of model coordinate)
  sradius = cradius       ! Support radius for 5th-order polynomial
                          ! or radius for 1/e for exponential weighting

!+++ Specific variables for observations

! *** Which observation type to assimilate
  assim_A = .true.        ! Whether to assimilation obervations of type A
  assim_B = .false.       ! Whether to assimilation obervations of type B

! *** Observation errors
  rms_obs_A = 0.5         ! Observation error standard deviation for observation A
  rms_obs_B = 0.25        ! Observation error standard deviation for observation B

! *** Obervation files
  file_obs_A = '../inputs_2fields/obsA.nc'  ! Observation file A
  file_obs_B = '../inputs_2fields/obsB.nc'  ! Observation file B

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


! *******************************************************
! *** Call PDAF initialization routine (all processes ***
! ***                                                 ***
! *** For all filters, PDAF_init is first called      ***
! *** specifying only the required parameters.        ***
! *** Further settings are done afterwards using      ***
! *** calls to PDAF_set_iparam & PDAF_set_rparam.     ***
! *******************************************************

  ! *** Here we specify only the required integer and real parameters
  ! *** Other parameters are set using calls to PDAF_set_iparam/PDAF_set_rparam
  pdaf_param_i(1) = dim_state_p ! State dimension
  pdaf_param_i(2) = dim_ens     ! Size of ensemble
  pdaf_param_r(1) = forget      ! Forgetting factor

  call PDAF3_init(filtertype, subtype, 0, &
       pdaf_param_i, 2,&
       pdaf_param_r, 1, &
       init_ens_pdaf, screen, status_pdaf)

  ! *** Additional parameter specifications ***
  ! *** -- These are all optional --        ***

  ! Generic settings
  call PDAF_set_iparam(5, type_forget, status_pdaf)      ! Type of forgetting factor
  call PDAF_set_iparam(9, type_obs_init, status_pdaf)    ! Initialize observation before or after call to prepoststep


! *** Check whether initialization of PDAF was successful ***
  if (status_pdaf /= 0) then
     write (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (Process ', mype_ens,')'
     call abort_parallel()
  end if


! **********************
! *** Initialize IAU ***
! **********************
  
  CALL PDAF_iau_init(type_iau, steps_iau, status_pdaf)


! **********************************
! *** Prepare ensemble forecasts ***
! **********************************

  call PDAF_init_forecast(next_observation_pdaf, distribute_state_pdaf, &
       prepoststep_pdaf, status_pdaf)


! ***************************************************
! *** Set coordinates of elements in state vector ***
! *** (used for localization in EnKF/ENSRF)       ***
! ***************************************************

!+++ Specific initialization for 2D tutorial model

  allocate(coords_p(n_dim, dim_state_p))

  ! For localization in EnKF and EnSRF/EAKF, PDAFomi_set_localize_covar
  ! is called in the observation modules. This routine requires a 
  ! coordinate array corresponding to the state vector.

  s = 0
  do k = 1, n_fields
     do i = 1, nx_p
        do j = 1, ny
           s = s + 1
           coords_p(1, s) = coords_x_p(i)
           coords_p(2, s) = coords_y_p(j)
        end do
     end do
  end do


! ************************************************************************
! *** Set domain coordinate limits (for use with OMI's use_global_obs) ***
! ************************************************************************
  
!+++ Specific initialization for 2D tutorial model

  ! Get offset of local domain in global domain in x-direction
  off_nx = nx_p*mype_filter

  lim_coords(1,1) = real(off_nx + 1)     ! West
  lim_coords(1,2) = real(off_nx + nx_p)  ! East
  lim_coords(2,1) = real(ny)             ! North
  lim_coords(2,2) = 1.0                  ! South

  call PDAFomi_set_domain_limits(lim_coords)

end subroutine init_pdaf
  ! *** Parse number of model tasks ***
  ! *** The module variable is N_MODELTASKS. Since it has to be equal
  ! *** to the ensemble size we parse dim_ens from the command line.

  handle = 'dim_ens'
  call parse(handle, dim_ens)
  if (online_coupling) n_modeltasks = dim_ens

  ! *** Fix number or model tasks for offline DA ***
  if (.not. online_coupling) n_modeltasks = 1

  ! *** Get rank and size of COMM_ensemble ***

  call MPI_Comm_Size(COMM_ensemble, npes_ens, MPIerr)
  call MPI_Comm_Rank(COMM_ensemble, mype_ens, MPIerr)

  ! *** Initialize communicators for ensemble evaluations ***
  if (mype_ens == 0) &
       write (*, '(/a, 2x, a)') 'model-PDAF', 'Initialize communicators for assimilation with PDAF'


  ! *** Check consistency of number of parallel ensemble tasks ***
  if (online_coupling) then
     consist1: if (n_modeltasks > npes_ens) then
        ! *** # parallel tasks is set larger than available Processs ***
        n_modeltasks = npes_ens
        if (mype_ens == 0) write (*, '(a, 3x, a)') &
             'model-PDAF', '!!! Resetting number of parallel ensemble tasks to total number of Processs!'
     end if consist1
     if (dim_ens > 0) then
        ! Check consistency with ensemble size
        consist2: if (n_modeltasks > dim_ens) then
           ! # parallel ensemble tasks is set larger than ensemble size
           n_modeltasks = dim_ens
           if (mype_ens == 0) write (*, '(a, 5x, a)') &
                'model-PDAF', '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
        end if consist2
     end if
  end if


  ! *** Store # Processs per ensemble                 ***
  ! *** used for info on Process 0 and for generation ***
  ! *** of model communicators on other Pes      ***
  allocate(local_npes_model(n_modeltasks))

  local_npes_model = floor(real(npes_ens) / real(n_modeltasks))
  do i = 1, (npes_ens - n_modeltasks * local_npes_model(1))
     local_npes_model(i) = local_npes_model(i) + 1
  end do


  ! ***              COMM_MODEL               ***
  ! *** Generate communicators for model runs ***
  ! *** (Split COMM_ENSEMBLE)                 ***
  pe_index = 0
  doens1: do i = 1, n_modeltasks
     do j = 1, local_npes_model(i)
        if (mype_ens == pe_index) then
           task_id = i
           exit doens1
        end if
        pe_index = pe_index + 1
     end do
  end do doens1


  call MPI_Comm_split(COMM_ensemble, task_id, mype_ens, &
       COMM_model, MPIerr)
  
  ! *** Re-initialize Process informations   ***
  ! *** according to model communicator ***
  call MPI_Comm_Size(COMM_model, npes_model, MPIerr)
  call MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

  if (screen > 1) then
    write (*,*) 'model-PDAF: mype(w)= ', mype_ens, '; model task: ', task_id, &
         '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
  end if


  ! Init flag FILTERProcess (all Processs of model task 1)
  if (task_id == 1) then
     filterpe = .true.
  else
     filterpe = .false.
  end if

  ! ***         COMM_FILTER                 ***
  ! *** Generate communicator for filter    ***
  ! *** For simplicity equal to COMM_couple ***

  if (filterpe) then
     my_color = task_id
  else
     my_color = MPI_UNDEFINED
  endif

  call MPI_Comm_split(COMM_ensemble, my_color, mype_ens, &
       COMM_filter, MPIerr)

  ! *** Initialize Process informations         ***
  ! *** according to coupling communicator ***
  if (filterpe) then
     call MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
     call MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)
  endif


  ! ***              COMM_COUPLE                 ***
  ! *** Generate communicators for communication ***
  ! *** between model and filter Processs             ***
  ! *** (Split COMM_ENSEMBLE)                    ***

  color_couple = mype_model + 1

  call MPI_Comm_split(COMM_ensemble, color_couple, mype_ens, &
       COMM_couple, MPIerr)

  ! *** Initialize Process informations         ***
  ! *** according to coupling communicator ***
  call MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
  call MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

  if (screen > 0) then
     if (mype_ens == 0) then
        write (*, '(/a, 2x, a)') 'model-PDAF Pconf', 'Process configuration:'
        write (*, '(a, 2x, a6, a9, a10, a14, a13, /a, 2x, a5, a9, a7, a7, a7, a7, a7, /a, 2x, a)') &
             'model-PDAF Pconf', 'world', 'filter', 'model', 'couple', 'filterPE', &
             'model-PDAF Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
             'model-PDAF Pconf', '----------------------------------------------------------'
     end if
     call MPI_Barrier(COMM_ensemble, MPIerr)
     if (task_id == 1) then
        write (*, '(a, 2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
             'model-PDAF Pconf', mype_ens, mype_filter, task_id, mype_model, color_couple, &
             mype_couple, filterpe
     endif
     if (task_id > 1) then
        write (*,'(a, 2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
         'model-PDAF Pconf', mype_ens, task_id, mype_model, color_couple, mype_couple, filterpe
     end if
     call MPI_Barrier(COMM_ensemble, MPIerr)

     if (mype_ens == 0) write (*, '(/a)') ''

  end if


! ***************************************************
! *** Provide parallelization information to PDAF ***
! ***************************************************

  call PDAF3_set_parallel(COMM_ensemble, COMM_model, COMM_filter, COMM_couple, &
       task_id, n_modeltasks, filterpe, flag)


! ******************************************************************************
! *** Initialize model equivalents to COMM_model, npes_model, and mype_model ***
! ******************************************************************************

  ! We use here generic names as defined as arguments. One could adapt them
  ! to the specific names using in the model code, but this is not required.

  model_comm = COMM_model
  model_comm_rank = mype_model
  model_comm_size = npes_model

end subroutine init_parallel_pdaf

