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
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_n_domains_pdaf(step, n_domains_p)

  use assimilation_pdaf_mod, &        ! Assimilation variables
       only: dim_state_p

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step        !< Current time step
  integer, intent(out) :: n_domains_p !< Process-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************
  
  ! Here simply the state dimension
  n_domains_p = dim_state_p

end subroutine init_n_domains_pdaf
