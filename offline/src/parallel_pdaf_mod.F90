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
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module parallel_pdaf_mod

  use mpi

  implicit none
  save 

  ! Parallelization variables that can be used in the user code

  ! Variables for each model task
  integer :: COMM_model         !< MPI communicator for model tasks
  integer :: mype_model         !< Number of Processs in COMM_model
  integer :: npes_model         !< Process rank in COMM_model

  ! Variables describing all processes involved in model integrations
  integer :: COMM_ensemble      !< Communicator for entire ensemble
  integer :: mype_ens           !< Rank in COMM_ensemble
  integer :: npes_ens           !< Size of COMM_ensemble

  ! Variables describing the processes involved in the analysis step
  integer :: COMM_assim         !< MPI communicator processes in analysis step
  integer :: npes_assim         !< Number of processes in COMM_da
  integer :: mype_assim         !< Process rank in COMM_da

  ! Additional variables for use with PDAF
  integer :: n_modeltasks = 1   !< Number of parallel model tasks
  integer :: task_id            !< Index of my model task (1,...,n_modeltasks)

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
    COMM_model = MPI_COMM_WORLD
   
  end subroutine init_parallel
!-------------------------------------------------------------------------------
!> Finalize MPI
!!
!! Routine to finalize MPI
!!
  subroutine finalize_parallel()

    implicit none
    
    integer :: MPIerr             !< Error flag for MPI

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
    
    integer :: MPIerr             !< Error flag for MPI

    call  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  end subroutine abort_parallel

end module parallel_pdaf_mod
