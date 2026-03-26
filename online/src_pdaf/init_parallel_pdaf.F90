!>  Interface routine to call initialization of parallelization for PDAF
!!
!! This routine stores the parallelization information provided by the
!! model. Afterwards, it calls the parallelization routine for PDAF, which 
!! initializes the communicators for handling the ensemble and the
!! the analysis of the data assimilation. This overwrites the communicator
!! provided by the model by the ensemble configuration.
!!
!! The parallelization variables returned from the PDAF initialization
!! routines are stored in variables from the module parallel_pdaf_mod
!! so that they can be used in the different user-provided routines.
!!
!! The routine is generic, but has to be compiled with the user code
!! because it uses the module parallel_pdaf_mod.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_parallel_pdaf(screen, model_comm, model_comm_rank, model_comm_size)

  use PDAF, &                     ! Command line parser
       only: PDAF_parse, PDAF3_init_parallel
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
  call PDAF3_init_parallel(screen, 0, 1, n_modeltasks, dim_ens, &
     COMM_model, mype_model, npes_model, &
     COMM_assim, mype_assim, npes_assim, &
     task_id)

  ! Update parallelization variables that are returned to model
  model_comm = COMM_model
  model_comm_rank = mype_model
  model_comm_size = npes_model

end subroutine init_parallel_pdaf
