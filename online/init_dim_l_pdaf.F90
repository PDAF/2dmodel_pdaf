!>  Set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in domain-localized filters
!!
!! The routine is called during analysis step in the 
!! loop over all local analysis domains. It sets
!! the dimension of the local model state on the
!! current analysis domain, its coordiantes and the
!! indices to map between the full state vector and 
!! the local state vector.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_dim_l_pdaf(step, domain_p, dim_l)

  use PDAF, &                       ! Routine to provide local indices to PDAF
       only: PDAFlocal_set_indices
  use model_pdaf_mod, &             ! Model variables
       only: ny, nx_p
  use assimilation_pdaf_mod, &      ! Variables for assimilation
       only: coords_l
  use parallel_pdaf_mod, &          ! assimilation parallelization variables
       only: mype_filter
  use statevector_pdaf_mod, &       ! State vector variables
       only: n_fields, id, sfields

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step     !< Current time step
  integer, intent(in)  :: domain_p !< Current local analysis domain
  integer, intent(out) :: dim_l    !< Local state dimension

! *** local variables ***
  integer :: i                       ! Counters
  integer :: off_p                   ! Process-local offset in global state vector
  integer, allocatable :: id_lstate_in_pstate(:) !< Indices of local state vector in Process-local global state vector


! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
  dim_l = n_fields


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Global coordinates of local analysis domain
  ! We use grid point indices as coordinates, but could e.g. use meters
  off_p = 0
  do i = 1, mype_filter
     off_p = off_p + nx_p*ny
  end do
  coords_l(1) = real(ceiling(real(domain_p+off_p)/real(ny)))
  coords_l(2) = real(domain_p+off_p) - (coords_l(1)-1)*real(ny)


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  allocate(id_lstate_in_pstate(dim_l))

  ! Here the local domain is a single grid point holding two variables
  ! The variables given by DOMAIN_P + offsets
  DO i=1, n_fields
     id_lstate_in_pstate(i) = domain_p + sfields(i)%off
  END DO

  ! Provide the index vector to PDAF
  call PDAFlocal_set_indices(dim_l, id_lstate_in_pstate)

  ! Deallocate index array
  deallocate(id_lstate_in_pstate)

end subroutine init_dim_l_pdaf
