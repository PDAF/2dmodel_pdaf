!> Get vector of synthetic observations from PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called when synthetic observations are generated with
!! PDAF. With the call, the user is provided with a generated observation
!! vector. This can then e.g. be written to a file using write_syn_obs.
!!
!! The routine is called by all filter processes.
!!
!! The routine is generic, but has to be part of the user code
!! because it is a call-back routine for PDAF.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code based on template
!! * Later revisions - see repository log
!!
subroutine get_obs_pdaf(step, dim_obs, observation)

  use synobs_pdaf_mod, &
       only: file_synobs, write_syn_obs

  implicit none

  ! *** Arguments
  integer, intent(in) :: step                 !< Current time step
  integer, intent(in) :: dim_obs              !< Dimension of obs. vector
  real, intent(in)    :: observation(dim_obs) !< Observation vector


! *************************
! *** store observation ***
! *************************

  ! Note 'file_synobs' is set in the observation module when
  ! initializing the output file

  call write_syn_obs(step, file_synobs, dim_obs, observation, 1)

end subroutine get_obs_pdaf

