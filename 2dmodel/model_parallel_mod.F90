!> Module for model parallelization
!!
!! This module provides variables for the MPI parallelization
!! of the tutorial model to be shared between model-related routines. 
!!
!! In addition, methods to initialize and finalize MPI are provided.
!!
!! Revision history:
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_parallel_mod

  use mpi

  implicit none
  save 

  ! Basic variables for model state integrations
  integer :: COMM_2Dmodel  !< MPI communicator for model tasks
  integer :: npes_2Dmodel  !< Number of Processs in COMM_model
  integer :: mype_2Dmodel  !< Process rank in COMM_model
  integer :: npes_world    !< Number of Processs in MPI_COMM_WORLD
  integer :: mype_world    !< Process rank in MPI_COMM_WORLD
  integer :: MPIerr        !< Error flag for MPI
  
contains
!-------------------------------------------------------------------------------
!> Initialize MPI
!!
!! Routine to initialize MPI, the number of Processs
!! (npes_world) and the rank of a Process (mype_world).
!! The model is executed within the scope of the
!! communicator Comm_model. It is also initialized
!! here together with its size (npes_model) and 
!! the rank of a Process (mype_model) within Comm_model.
!!
  subroutine init_parallel()

    implicit none

    integer :: i
  
    call MPI_INIT(i);
    call MPI_Comm_Size(MPI_COMM_WORLD,npes_world,i)
    call MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,i)

    ! Initialize model communicator, its size and the process rank
    ! Here the same as for MPI_COMM_WORLD
    Comm_2Dmodel = MPI_COMM_WORLD
    npes_2Dmodel = npes_world
    mype_2Dmodel = mype_world

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

end module model_parallel_mod
