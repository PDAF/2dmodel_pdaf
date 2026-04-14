!>  Initialize model fields from state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all DA methods.
!!
!! During the forecast phase, this subroutine is called from PDAF
!! providing a model state, which has to be integrated.
!! The routine has to initialize the fields of the model (typically
!! available through a module) from the provided state vector.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine distribute_state_pdaf(dim_p, state_p)

  use statevector_pdaf_mod, &            ! State vector variables
       only: id, sfields

  ! Include model variables
!   use model_pdaf_mod, &                  ! Model variables
!        only: field_X, field_Y

  implicit none
  
! *** Arguments ***
  integer, intent(in) :: dim_p           !< Process-local state dimension
  real, intent(inout) :: state_p(dim_p)  !< Process-local state vector

! *** local variables ***


! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE distribute_state_pdaf.F90: Implement initialization of model fields here!'

!  FIELD_X = state_p

end subroutine distribute_state_pdaf
