!>  Interface routine to call initialization of parallelization for PDAF
!!
!! This routine stores the parallelization information provided by the
!! model (which after the communicator initialization for PDAF holds
!! the joint communicator for all processes running model integrations). 
!! Afterwards, it call the parallelization routine for PDAF, which 
!! initializes the communicators for handling the ensemble and 
!! the analysis of the data assimilation.
!!
!! The parallelization variables returned from the PDAF initialization
!! routines are stored in the module parallel_pdaf_mod so that they can
!! be used in the different user-provided routines.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_parallel_pdaf(screen, model_comm, model_comm_rank, model_comm_size)

  use PDAF, &                     ! Command line parser
       only: PDAF_parse
  use parallel_pdaf_mod, &        ! PDAF parallelization variables
       only: n_modeltasks, task_id, mype_ens, npes_ens, COMM_ensemble, &
       mype_model, npes_model, COMM_model, mype_assim, npes_assim, COMM_assim

  implicit none

! *** Arguments ***
  integer, intent(in)    :: screen           !< Whether screen information is shown

  ! Model parallelization variables (one can keep these generic names)
  integer, intent(inout) :: model_comm       !< Model MPI communicator for model tasks
  integer, intent(inout) :: model_comm_size  !< Number of processes in model_comm
  integer, intent(inout) :: model_comm_rank  !< Process rank in model_comm

! *** Local variables ***
  integer :: dim_ens                         ! Ensemble size
  character(len=32) :: handle                ! Handle for command line parser


  ! Parse ensemble size
  handle = 'dim_ens'
  call PDAF_parse(handle, dim_ens)

  ! Set number of model tasks for fully-parallel mode
  n_modeltasks = dim_ens

  ! Store the parallelization variables provided by the model
  ! At this point, they describe all processes doing model integrations
  COMM_ensemble = model_comm
  mype_ens = model_comm_rank
  npes_ens = model_comm_size

  ! Initialize ensemble parallelization
  call init_parallel_pdaf_doinit(screen, 0, 1, n_modeltasks, dim_ens, &
     COMM_model, mype_model, npes_model, &
     COMM_assim, mype_assim, npes_assim, &
     task_id)

  ! Update parallelization variables that are returned to model
  model_comm = COMM_model
  model_comm_rank = mype_model
  model_comm_size = npes_model

end subroutine init_parallel_pdaf


!>  Initialize communicators for PDAF
!!
!! Parallelization routine for a model with attached PDAF. The subroutine is
!! called in the main program subsequently to the initialization of MPI. It
!! initializes MPI communicators for the model tasks, assimilation task and the
!! coupling between model and assimilation tasks. In addition some other variables 
!! for the parallelization are initialized.
!! The communicators and variables are handed over to PDAF in the call to 
!! PDAF_set_parallel toward the end of this routine.
!!
!! 3 Communicators are generated:
!! * _COMM_assim_: Communicator in which the assimilation analysis is computed
!! * _COMM_model_: Communicators for parallel model forecasts
!! * _COMM_couple_: Communicator for coupling between model and assi. processes
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
!! * _assimpe_ - Logical: Does the Process execute the analysis step?
!! * _task_id_ - Integer: Index identifying the model task
!! * _my_ensemble_ - Integer: The index of the Process's model task
!! * _local_npes_model_ - Integer array holding numbers of Processs per model task
!!
!! For COMM_assim and COMM_model also the size of the communicators
!! (npes_assim and npes_model) and the rank of each process  (mype_assim,
!! mype_model) are initialized. These variables can be used in the model part 
!! of the program, but are not handed over to PDAF.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_parallel_pdaf_doinit(screen, type_parallel, online_coupling, n_modeltasks, dim_ens, &
     COMM_model, mype_model, npes_model, COMM_assim, mype_assim, npes_assim, &
     task_id)

  use mpi                         ! MPI
  use PDAF, &                     ! PDAF routines
       only: PDAF3_set_parallel

  implicit none

! *** Arguments ***
  integer, intent(in)    :: screen           !< Whether screen information is shown

  ! Model variables for parallelization (one can keep these generic names)
  integer, intent(in) :: type_parallel       !< Type of parallelization
  integer, intent(in) :: online_coupling     !< 1: online DA coupling, 0: offline DA coupling
  integer, intent(inout) :: n_modeltasks     !< Number of model tasks
  integer, intent(inout) :: dim_ens          !< Ensemble size / number of model tasks
  integer, intent(out) :: COMM_model         !< Model MPI communicator for model tasks
  integer, intent(out) :: npes_model         !< Number of Processs in COMM_model
  integer, intent(out) :: mype_model         !< Process rank in COMM_model
  integer, intent(out) :: COMM_assim         !< MPI communicator for assimilation processes 
  integer, intent(out) :: npes_assim         !< Number of processes in COMM_assim
  integer, intent(out) :: mype_assim         !< Process rank in COMM_assim
  integer, intent(out) :: task_id            !< Index of my model task (1,...,n_modeltasks)

! *** Local variables ***
  integer :: i, j                     ! Counters
  integer :: pe_index                 ! Index of Process
  integer :: my_color, color_couple   ! Variables for communicator-splitting 
  integer :: flag                     ! Status flag
  character(len=32) :: handle         ! Handle for command line parser
  integer :: MPIerr                   ! Error flag for MPI
  integer, allocatable :: local_npes_model(:) ! Number of processes per ensemble
  integer :: COMM_couple              ! MPI communicator for coupling assim and model
  integer :: mype_couple              ! Rank in COMM_couple
  integer :: npes_couple              ! Size in COMM_couple
  integer :: COMM_ensemble            ! Communicator for entire ensemble
  integer :: mype_ens                 ! Rank in COMM_ensemble
  integer :: npes_ens                 ! Size of COMM_ensemble
  logical :: assimpe                 ! Whether we are on a Process in a COMM_assim


  ! *** Define ensemble communicator.                   ***
  ! *** This is the communicator in which PDAF operates ***

  COMM_ensemble = COMM_model

  ! *** Fix number or model tasks for offline DA ***
  if (online_coupling==0) n_modeltasks = 1

  ! *** Get rank and size of COMM_ensemble ***

  call MPI_Comm_Size(COMM_ensemble, npes_ens, MPIerr)
  call MPI_Comm_Rank(COMM_ensemble, mype_ens, MPIerr)

  ! *** Initialize communicators for ensemble evaluations ***
  if (mype_ens == 0) &
       write (*, '(/a, 2x, a)') 'PDAF', 'Initialize communicators for assimilation with PDAF'


  ! *** Check consistency of number of parallel ensemble tasks ***
  if (online_coupling==1) then
     consist1: if (n_modeltasks > npes_ens) then
        ! *** # parallel tasks is set larger than available Processs ***
        n_modeltasks = npes_ens
        if (mype_ens == 0) write (*, '(a, 3x, a)') &
             'PDAF', '!!! Resetting number of parallel ensemble tasks to total number of Processs!'
     end if consist1
     if (dim_ens > 0) then
        ! Check consistency with ensemble size
        consist2: if (n_modeltasks > dim_ens) then
           ! # parallel ensemble tasks is set larger than ensemble size
           n_modeltasks = dim_ens
           if (mype_ens == 0) write (*, '(a, 5x, a)') &
                'PDAF', '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
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
    write (*,*) 'PDAF: mype(w)= ', mype_ens, '; model task: ', task_id, &
         '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
  end if


  ! Init flag for assim processes (all processes of model task 1)
  if (task_id == 1) then
     assimpe = .true.
  else
     assimpe = .false.
  end if

  ! ***         COMM_ASSIM                  ***
  ! *** Generate communicator for analysis  ***
  ! *** For simplicity equal to COMM_couple ***

  if (assimpe) then
     my_color = task_id
  else
     my_color = MPI_UNDEFINED
  endif

  call MPI_Comm_split(COMM_ensemble, my_color, mype_ens, &
       COMM_assim, MPIerr)

  ! *** Initialize Process informations         ***
  ! *** according to coupling communicator ***
  if (assimpe) then
     call MPI_Comm_Size(COMM_assim, npes_assim, MPIerr)
     call MPI_Comm_Rank(COMM_assim, mype_assim, MPIerr)
  endif


  ! ***              COMM_COUPLE                 ***
  ! *** Generate communicators for communication ***
  ! *** between model and assim processes        ***
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
        write (*, '(/a, 2x, a)') 'PDAF Pconf', 'Process configuration:'
        write (*, '(a, 2x, a6, a9, a10, a14, a13, /a, 2x, a5, a9, a7, a7, a7, a7, a7, /a, 2x, a)') &
             'PDAF Pconf', 'world', 'assim', 'model', 'couple', 'assimPE', &
             'PDAF Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
             'PDAF Pconf', '----------------------------------------------------------'
     end if
     call MPI_Barrier(COMM_ensemble, MPIerr)
     if (task_id == 1) then
        write (*, '(a, 2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
             'PDAF Pconf', mype_ens, mype_assim, task_id, mype_model, color_couple, &
             mype_couple, assimpe
     endif
     if (task_id > 1) then
        write (*,'(a, 2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
         'PDAF Pconf', mype_ens, task_id, mype_model, color_couple, mype_couple, assimpe
     end if
     call MPI_Barrier(COMM_ensemble, MPIerr)

     if (mype_ens == 0) write (*, '(/a)') ''

  end if


! ***************************************************
! *** Provide parallelization information to PDAF ***
! ***************************************************

  call PDAF3_set_parallel(COMM_ensemble, COMM_model, COMM_assim, COMM_couple, &
       task_id, n_modeltasks, assimpe, flag)

end subroutine init_parallel_pdaf_doinit

