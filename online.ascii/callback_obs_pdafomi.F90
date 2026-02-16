!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!!   When adding an observation type, one has to add one module
!!   obs_OBSTYProcess_pdafomi (based on the template obs_OBSTYProcess_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files. 
!! 
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYProcess.
!!
subroutine init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  use obs_A_pdafomi, only: assim_A, init_dim_obs_A
  use obs_B_pdafomi, only: assim_B, init_dim_obs_B

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step     !< Current time step
  integer, intent(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  integer :: dim_obs_type ! Dimension of one observation type


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_type = 0
  dim_obs = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called

  if (assim_A) then
     call init_dim_obs_A(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if
  if (assim_B) then
     call init_dim_obs_B(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if

end subroutine init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYProcess.
!!
subroutine obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  use obs_A_pdafomi, only: obs_op_A
  use obs_B_pdafomi, only: obs_op_B

  implicit none

! *** Arguments ***
  integer, intent(in) :: step                 !< Current time step
  integer, intent(in) :: dim_p                !< Process-local state dimension
  integer, intent(in) :: dim_obs              !< Dimension of full observed state
  real, intent(in)    :: state_p(dim_p)       !< Process-local model state
  real, intent(inout) :: ostate(dim_obs)      !< Process-local full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi

  call obs_op_A(dim_p, dim_obs, state_p, ostate)
  call obs_op_B(dim_p, dim_obs, state_p, ostate)

end subroutine obs_op_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type.
!!
subroutine init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! Include functions for different observations
  use obs_A_pdafomi, only: init_dim_obs_l_A
  use obs_B_pdafomi, only: init_dim_obs_l_B
  
  implicit none

! *** Arguments ***
  integer, intent(in)  :: domain_p   !< Index of current local analysis domain
  integer, intent(in)  :: step       !< Current time step
  integer, intent(in)  :: dim_obs    !< Full dimension of observation vector
  integer, intent(out) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Call init_dim_obs_l specific for each observation

  call init_dim_obs_l_A(domain_p, step, dim_obs, dim_obs_l)
  call init_dim_obs_l_B(domain_p, step, dim_obs, dim_obs_l)

end subroutine init_dim_obs_l_pdafomi
