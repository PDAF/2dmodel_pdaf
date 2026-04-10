!>  Initialize model fields from state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all DA methods.
!!
!! During the forecast phase of the filter this subroutine is called from
!! PDAF_init_forecast or PDAF3_assimilate  supplying a model state, which
!! has to be integrated. The routine has to initialize the fields of the 
!! model (typically available through a module) from the state vector 
!! provided by PDAF.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine distribute_state_pdaf(dim_p, state_p)

  use statevector_pdaf_mod, &            ! State vector variables
       only: id, sfields

  ! Specific for model
  use model_pdaf_mod, &                  ! Model variables
       only: nx_p, ny, fieldA_p, fieldB_p

  implicit none
  
! *** Arguments ***
  integer, intent(in) :: dim_p           !< Process-local state dimension
  real, intent(inout) :: state_p(dim_p)  !< Process-local state vector

! *** local variables ***
  integer :: i, j, s         ! Counters


! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************

  ! +++ Note on counter s:
  ! +++ Using the counter s looks primitive, but it
  ! +++ makes the code fail-save because it avoids
  ! +++ index calculations involving nx_p or ny.

!+++ Specific part for 2D tutorial model

  ! FieldA
  s = sfields(id%fieldA)%off
  do j = 1, nx_p
     do i = 1, ny
        s = s + 1
        fieldA_p(i, j) = state_p(s)
     end do
  end do

  ! FieldB
  s = sfields(id%fieldB)%off
  do j = 1, nx_p
     do i = 1, ny
        s = s + 1
        fieldB_p(i, j) = state_p(s)
     end do
  end do

!+++ End of specific part

end subroutine distribute_state_pdaf
