!> Module for 2D tutorial model
!!
!! This module provides variables for the 
!! 2-dimensional tutorial model with parallelization.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_mod

  implicit none
  save
  public

! *** Variables specific for 2D tutorial model ***

  integer :: nx                      !< Size of 2D grid in x-direction
  integer :: ny                      !< Size of 2D grid in y-direction
  integer :: total_steps             !< Total number of time steps
  real, allocatable :: fieldA_p(:,:) !< Process-local part of model field A
  real, allocatable :: fieldB_p(:,:) !< Process-local part of model field B

  integer :: nx_p                    !< Process-local size in x-direction

end module model_mod
