!> Set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in domain-localized filters and the 3DEnVar
!! and hybrid 3DVAR.
!!
!! The routine is called in PDAF at the beginning of
!! the analysis step before the loop through all local
!! analysis domains. It has to set the number of local
!! analysis domains for the process-local domain.
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine init_n_domains_cb_pdaf(step, n_domains_p)

  ! Specific for model
!   use model_pdaf_mod, &               ! Model-related variables
!        only: nx_p, ny

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step        !< Current time step
  integer, intent(out) :: n_domains_p !< Process-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************

  ! Template reminder - delete when implementing functionality
  write (*,*) 'TEMPLATE init_n_domains_pdaf.F90: Set number of local analysis domains here!'

!  n_domains_p = ?

  ! dummy initialization to allow testing
  n_domains_p = 1

end subroutine init_n_domains_cb_pdaf
