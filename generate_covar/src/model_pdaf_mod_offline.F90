!> Module to declare model variables for offline coupling
!!
!! This module declares model-related variables for offline
!! coupled DA. While for online coupling the variables
!! would be included by 'use' statements from model modules,
!! we have to declare the variables explicitly of the 
!! offline coupling.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code from restructuring
!! * Later revisions - see repository log
!!
module model_pdaf_mod

! In the online-coupled mode, we would include this from the model
!   use model_mod, &                ! Model variables
!        ONLY: nx, ny, nx_p, n_dim, offset_x_p, &
!        coords_x_p, coords_y_p, total_steps, fieldA_p, fieldB_p

  implicit none

! *** Declare variables specific for 2D tutorial model ***

  integer :: nx                      !< Size of 2D grid in x-direction
  integer :: ny                      !< Size of 2D grid in y-direction
  integer :: n_dim                   !< Number of model dimensions
  real, allocatable :: coords_x_p(:) !< Process-local coordinates in x-direction
  real, allocatable :: coords_y_p(:) !< Process-local coordinates in y-direction

  integer :: nx_p                    !< Process-local size in x-direction
  integer :: offset_x_p              !< Offset of sub-domain inglobal grid


! *** Additional variables ***

! One might like to declare here additional variables relating
! to the model fields or model grid only used for DA.

end module model_pdaf_mod

 
