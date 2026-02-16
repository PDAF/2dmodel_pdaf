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
       only: nx_p, ny, fieldA_p, fieldB_p
  use statevector_pdaf_mod, &       ! State vector variables
       only: id, sfields

  implicit none
  
! *** Arguments ***
  integer, intent(in) :: dim_p           !< Process-local state dimension
  real, intent(inout) :: state_p(dim_p)  !< local state vector

! *** local variables ***
  integer :: i, j, s         ! Counters
  

! *************************************************
! *** Initialize state vector from model fields ***
! *** for process-local model domain            ***
! *************************************************

  ! +++ Note on counter s:
  ! +++ Using the counter s looks primitive, but it
  ! +++ makes the code fail-save because it avoids
  ! +++ index calculations involving nx_p or ny.

  ! FieldA
  s = sfields(id%fieldA)%off
  do j = 1, nx_p
     do i = 1, ny
        s = s + 1
        state_p(s) = fieldA_p(i, j)
     end do
  end do

  ! FieldB
  s = sfields(id%fieldB)%off
  do j = 1, nx_p
     do i = 1, ny
        s = s + 1
        state_p(s) = fieldB_p(i, j)
     end do
  end do

end subroutine collect_state_pdaf
