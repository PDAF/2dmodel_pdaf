!>  Initialize state vector from model fields
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all DA methods.
!!
!! This subroutine is called during the forecast 
!! phase from PDAF_put_state_X or PDAF_assimilate_X
!! after the propagation of each ensemble member. 
!! The supplied state vector has to be initialized
!! from the model fields (typically via a module). 
!! With parallelization, MPI communication might be 
!! required to initialize state vectors for all 
!! subdomains on the model Processs. 
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine collect_state_pdaf(dim_p, state_p)

  use model_pdaf_mod, &             ! Model variables
       only: nx_p, ny, fieldA_p

  implicit none
  
! *** Arguments ***
  integer, intent(in) :: dim_p           !< Process-local state dimension
  real, intent(inout) :: state_p(dim_p)  !< local state vector

! *** local variables ***
  integer :: j         ! Counters
  

! *************************************************
! *** Initialize state vector from model fields ***
! *** for process-local model domain            ***
! *************************************************

  ! + For the 2D tutorial model the state vector and
  ! + the model field are identical. Hence, state vector
  ! + directly initialized from the model field by
  ! + each model Process.

  do j = 1, nx_p
     state_p(1 + (j-1)*ny : j*ny) = fieldA_p(1:ny, j)
  end do

end subroutine collect_state_pdaf
