!>  Initialize state vector from model fields
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used by all DA methods.
!!
!! This subroutine is called during the forecast phase from PDAF
!! after the propagation of each ensemble member. 
!! The supplied state vector has to be initialized from the model fields
!! (typically accessible via a module).
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine coll_state_cb_pdaf(dim_p, state_p)

  use statevector_pdaf_mod, &            ! State vector variables
       only: id, sfields

  ! Include model variables
!   use model_pdaf_mod, &                  ! Model variables
!        only: field_X, field_Y

  implicit none
  
! *** Arguments ***
  integer, intent(in) :: dim_p           !< Process-local state dimension
  real, intent(inout) :: state_p(dim_p)  !< local state vector

! *** local variables ***
  

! *************************************************
! *** Initialize state vector from model fields ***
! *** for process-local model domain            ***
! *************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE collect_state_pdaf.F90: Implement initialization of state vector here!'

!   state_p = field_X ????

end subroutine coll_state_cb_pdaf
