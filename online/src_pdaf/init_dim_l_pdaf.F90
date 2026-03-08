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
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_dim_l_pdaf(step, domain_p, dim_l)

  use PDAF, &                       ! Routine to provide local indices to PDAF
       only: PDAFlocal_set_indices
  use assimilation_pdaf_mod, &      ! Variables for assimilation
       only: coords_l
  use statevector_pdaf_mod, &       ! State vector variables
       only: n_fields, sfields

  ! Specific for 2D tutorial model
  use model_pdaf_mod, &             ! Model variables
       only: ny, coords_x_p, coords_y_p

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step      !< Current time step
  integer, intent(in)  :: domain_p  !< Current local analysis domain
  integer, intent(out) :: dim_l     !< Local state dimension

! *** local variables ***
  integer :: i                      ! Counters
  integer :: idx(2)                 ! Grid point indices for index domain_p
  integer, allocatable :: id_lstate_in_pstate(:) !< Indices of local state vector in Process-local global state vector


! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
!+++ Specific for 2D tutorial model

  dim_l = n_fields


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

!+++ Specific part for 2D tutorial model

  ! Determine grid point indices from domain_p
  idx(1) = ceiling(real(domain_p)/real(ny))
  idx(2) = domain_p - (idx(1)-1)*ny

  ! Initialize local coordinates
  coords_l(1) = coords_x_p(idx(1))
  coords_l(2) = coords_y_p(idx(2))

!+++ End of specific part

! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  allocate(id_lstate_in_pstate(dim_l))

!+++ Specific part for 2D tutorial model

  ! Here the local domain is a single grid point holding two variables
  ! The variables given by DOMAIN_P + offsets
  DO i=1, n_fields
     id_lstate_in_pstate(i) = domain_p + sfields(i)%off
  END DO

!+++ End of specific part

  ! Provide the index vector to PDAF
  call PDAFlocal_set_indices(dim_l, id_lstate_in_pstate)

  ! Deallocate index array
  deallocate(id_lstate_in_pstate)

end subroutine init_dim_l_pdaf
