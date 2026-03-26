!> Module for ensemble parallelization
!!
!! This module provides variables for the MPI parallelization
!! to be shared between model-related routines. The are variables
!! that are used in the model, even without PDAF and additional
!! variables that are only used, if data assimilation with PDAF
!! is performed.
!!
!! In addition methods to initialize and finalize MPI are provided.
!! The initialization routine is only for the model itself, the 
!! more complex initialization of communicators for execution with
!! PDAF is peformed in init_parallel_pdaf.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module parallel_pdaf_mod

  use mpi

  implicit none
  save 

  ! Basic variables for model state integrations
  integer :: COMM_model         !< MPI communicator for model tasks
  integer :: mype_model         !< Number of Processs in COMM_model
  integer :: npes_model         !< Process rank in COMM_model

  integer :: COMM_ensemble      !< Communicator for entire ensemble
  integer :: mype_ens           !< Rank in COMM_ensemble
  integer :: npes_ens           !< Size of COMM_ensemble

  ! Additional variables for use with PDAF
  integer :: n_modeltasks = 1   !< Number of parallel model tasks

  integer :: COMM_filter        !< MPI communicator for filter Processs 
  integer :: npes_filter        !< Number of processes in COMM_filter
  integer :: mype_filter        !< Process rank in COMM_filter

  integer :: task_id            !< Index of my model task (1,...,n_modeltasks)

  integer :: MPIerr             !< Error flag for MPI

contains
!-------------------------------------------------------------------------------
!> Initialize MPI
!!
!! Routine to initialize MPI, the number of processes
!! (npes_model) and the rank of a process (mype_model).
!! The model is executed within the scope of the
!! communicator Comm_model.
!!
  subroutine init_parallel()

    implicit none

    integer :: i
  
    call MPI_INIT(i);
    call MPI_Comm_Size(MPI_COMM_WORLD,npes_model,i)
    call MPI_Comm_Rank(MPI_COMM_WORLD,mype_model,i)

    ! Initialize model communicator
    ! Here the same as for MPI_COMM_WORLD
    Comm_model = MPI_COMM_WORLD
   
  end subroutine init_parallel
!-------------------------------------------------------------------------------
!> Finalize MPI
!!
!! Routine to finalize MPI
!!
  subroutine finalize_parallel()

    implicit none
    
    call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
    call  MPI_Finalize(MPIerr)

  end subroutine finalize_parallel
!-------------------------------------------------------------------------------
!> Abort MPI
!!
!! Routine for abort MPI program.
!!
  subroutine abort_parallel()

    implicit none
    
    call  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  end subroutine abort_parallel

end module parallel_pdaf_mod
