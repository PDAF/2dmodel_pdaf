!>  Initialize model fields from state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all DA methods.
!!
!! During the forecast phase of the filter this subroutine
!! is called from PDAF_init_forecast or PDAF3_assimilate.
!! supplying a model state, which has to be evolved. 
!! The routine has to initialize the fields of the 
!! model (typically available through a module) from 
!! the state vector of PDAF. With parallelization, 
!! MPI communication might be required to 
!! initialize all subdomains on the model Processs.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine distribute_state_pdaf(dim_p, state_p)

  use model_pdaf_mod, &             ! Model variables
       only: nx_p, ny, fieldA_p

  implicit none
  
! *** Arguments ***
  integer, intent(in) :: dim_p           !< Process-local state dimension
  real, intent(inout) :: state_p(dim_p)  !< Process-local state vector

! *** local variables ***
  integer :: j         ! Counters


! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************

  ! + For the 2D tutorial model the state vector and
  ! + the model field are identical. Hence, the field
  ! + array is directly initialized from an ensemble 
  ! + state vector by each model Process.

  do j = 1, nx_p
     fieldA_p(1:ny, j) = state_p(1 + (j-1)*ny : j*ny)
  end do

end subroutine distribute_state_pdaf
