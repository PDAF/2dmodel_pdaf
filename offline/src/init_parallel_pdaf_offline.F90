!>  Interface routine to call initialization of parallelization for PDAF
!!
!! This routine calls the parallelization routine for PDAF, which 
!! initializes the communicators for handling the analysis step
!! of the data assimilation.
!!
!! The parallelization variables returned from the PDAF initialization
!! routines are stored in the module parallel_pdaf_mod so that they can
!! be used in the different user-provided routines.
!!
!! The routine is generic, but has to be compiled with the user code
!! because it uses the module parallel_pdaf_mod.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_parallel_pdaf_offline(screen)

  use mpi
  use PDAF, &                     ! Command line parser
       only: PDAF_parse, PDAF3_init_parallel
  use parallel_pdaf_mod, &        ! PDAF parallelization variables
       only: n_modeltasks, task_id, mype_ens, npes_ens, COMM_ensemble, &
       mype_model, npes_model, COMM_model, mype_assim, npes_assim, COMM_assim

  implicit none

! *** Arguments ***
  integer, intent(in)    :: screen           !< Whether screen information is shown

! *** Local variables ***
  integer :: dim_ens                         ! Ensemble size
  character(len=32) :: handle                ! Handle for command line parser


  ! Parse ensemble size
  handle = 'dim_ens'
  call PDAF_parse(handle, dim_ens)

  ! Set number of model tasks for offline mode
  n_modeltasks = 1

  ! Initialize variables for calling initialization
  COMM_model = MPI_COMM_WORLD
  mype_model = 0
  npes_model = 0

  ! Initialize ensemble parallelization
  call PDAF3_init_parallel(screen, 0, 0, n_modeltasks, dim_ens, &
     COMM_model, mype_model, npes_model, &
     COMM_assim, mype_assim, npes_assim, &
     task_id)

  ! Initialize variables for all processes as they are used in some routines
  COMM_ensemble = COMM_model
  mype_ens = mype_model
  npes_ens = npes_model

end subroutine init_parallel_pdaf_offline
