!> Module for ensemble parallelization
!!
!! This module provides variables for the MPI parallelization
!! to be shared between model-related routines. these are the
!! variables typically used in the user code, while a parallelized
!! model has its own set of variables.
!!
!! In addition, methods to initialize and finalize MPI are provided.
!! Then can be used, e.g. when coupling PDAF for a model this is not
!! yet parallelized. The initialization routine (init_parallel) only
!! perform the overall initialization. The initialization of the
!! ensemble parallelization for PDAF is performed in init_parallel_pdaf.
!!
!! The module is generic, but has to be part of the user code
!! to provide the freedom to adapt it.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
module parallel_pdaf_mod

  use mpi

  implicit none
  save 

  ! Parallelization variables that can be used in the user code

  ! Variables for each model task
  integer :: COMM_model=0         !< MPI communicator for model tasks
  integer :: myproc_model=0       !< Process rank in COMM_model
  integer :: nproc_model=1        !< Number of Processses in COMM_model

  ! Variables describing all processes involved in model integrations
  integer :: COMM_ens=0           !< Jont Communicator for entire ensemble
  integer :: myproc_ens=0         !< Rank in COMM_ens
  integer :: nproc_ens=1          !< Size of COMM_ens

  ! Variables describing the processes involved in the analysis step
  integer :: COMM_assim=0         !< MPI communicator processes in analysis step
  integer :: myproc_assim=1       !< Process rank in COMM_da
  integer :: nproc_assim=0        !< Number of processes in COMM_da

  ! Additional variables for use with PDAF
  integer :: n_modeltasks=1       !< Number of parallel model tasks
  integer :: task_id=1            !< Index of my model task (1,...,n_modeltasks)

end module parallel_pdaf_mod
