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
!! The routine is generic, but has to be part of the user code
!! because it uses the module parallel_pdaf_mod.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * 2026-02 - Lars Nerger - Revision for using PDAF3_init_forecast 
!! * Later revisions - see repository log
!!
subroutine init_parallel_pdaf_offline(screen)

  use mpi
  use PDAF, &                     ! Command line parser
       only: PDAF3_init_parallel
  use parallel_pdaf_mod, &        ! PDAF parallelization variables
       only: n_modeltasks, task_id, myproc_ens, nproc_ens, COMM_ens, &
       myproc_model, nproc_model, COMM_model, myproc_assim, nproc_assim, COMM_assim

  implicit none

! *** Arguments ***
  integer, intent(in)    :: screen           !< Whether screen information is shown

  ! Set number of model tasks for offline mode
  n_modeltasks = 1

  ! Initialize ensemble parallelization
  call PDAF3_init_parallel(screen, 0, 0, 0, n_modeltasks, &
     COMM_model, myproc_model, nproc_model, &
     COMM_assim, myproc_assim, nproc_assim, &
     task_id)

  ! Initialize variables for all processes as they are used in some routines
  COMM_ens   = COMM_model
  myproc_ens = myproc_model
  nproc_ens  = nproc_model

end subroutine init_parallel_pdaf_offline
