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
!!   obs_OBSTYPE_pdafomi (based on the template obs_OBSTYPE_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files. 
!! 
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_OBSTYPE.
!!
subroutine init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  use obs_OBSTYPE_pdafomi, only: assim_OBSTYPE, init_dim_obs_OBSTYPE

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step     !< Current time step
  integer, intent(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  integer :: dim_obs_type ! Dimension of one observation type


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Template reminder - delete when implementing functionality
  write (*, *) 'TEMPLATE callback_obs_pdafomi.F90/init_dim_obs_pdafomi: complete interface to observation modules'

  ! Initialize number of observations
  dim_obs_type = 0
  dim_obs = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called

  if (assim_OBSTYPE) then
     call init_dim_obs_OBSTYPE(step, dim_obs_type)
     dim_obs = dim_obs + dim_obs_type
  end if

end subroutine init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_OBSTYPE.
!!
subroutine obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  use obs_OBSTYPE_pdafomi, only: obs_op_OBSTYPE

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

  ! Template reminder - delete when implementing functionality
  write (*, *) 'TEMPLATE callback_obs_pdafomi.F90/obs_op_pdafomi: complete interface to observation modules'

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi

  call obs_op_OBSTYPE(dim_p, dim_obs, state_p, ostate)

end subroutine obs_op_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_l_OBSTYPE.
!!
subroutine init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! Include functions for different observations
  use obs_OBSTYPE_pdafomi, only: init_dim_obs_l_OBSTYPE
  
  implicit none

! *** Arguments ***
  integer, intent(in)  :: domain_p   !< Index of current local analysis domain
  integer, intent(in)  :: step       !< Current time step
  integer, intent(in)  :: dim_obs    !< Full dimension of observation vector
  integer, intent(out) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Template reminder - delete when implementing functionality
  write (*, *) 'TEMPLATE callback_obs_pdafomi.F90/init_dim_obs_l_pdafomi: complete interface to observation modules'

  ! Call init_dim_obs_l specific for each observation

  call init_dim_obs_l_OBSTYPE(domain_p, step, dim_obs, dim_obs_l)

end subroutine init_dim_obs_l_pdafomi
