!>  Set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in domain-localized filters
!!
!! The routine is called during analysis step in the loop over all local
!! analysis domains. It sets the dimension of the local model state on the
!! current analysis domain, its coordiantes and the indices to map between
!! the full state vector and the local state vector.
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine init_dim_l_cb_pdaf(step, domain_p, dim_l)

  use PDAF, &                       ! Routine to provide local indices to PDAF
       only: PDAFlocal_set_indices
  use statevector_pdaf_mod, &       ! State vector variables
       only: n_fields, sfields

  ! Specific for model
  use model_pdaf_mod, &             ! Model variables
       only: coords_l, ny, coords_x_p, coords_y_p

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step      !< Current time step
  integer, intent(in)  :: domain_p  !< Current local analysis domain
  integer, intent(out) :: dim_l     !< Local state dimension

! *** local variables ***
  integer, allocatable :: id_lstate_in_pstate(:) ! Indices of local state vector in process-local global state vector


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! Template reminder - delete when implementing functionality
  write (*,*) 'TEMPLATE init_dim_l_pdaf.F90: Set local state dimension here!'

!  dim_l = ??


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

    ! Template reminder - delete when implementing functionality
  write (*,*) 'TEMPLATE init_dim_l_pdaf.F90: Ensure that coords_l is declared correctly in model_pdaf_mod!'

  ! Global coordinates of local analysis domain
  ! these can be determined from coords_x_p and coords_y_p

!  coords_l(1) = ??
!  coords_l(2) = ??


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  allocate(id_lstate_in_pstate(dim_l))

!  id_lstate_in_pstate = ??

  ! Provide the index vector to PDAF
  call PDAFlocal_set_indices(dim_l, id_lstate_in_pstate)

  ! Deallocate index array
  deallocate(id_lstate_in_pstate)

end subroutine init_dim_l_cb_pdaf
