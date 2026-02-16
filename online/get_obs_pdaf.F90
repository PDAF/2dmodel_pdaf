!> Get vector of synthetic observations from PDAF
!!
!! User-supplied routine for PDAF.
!! The routine is called when synthetic observations
!! are generated with PDAF. With the call, the user
!! is provided with a generated observation vetor. 
!! This can then e.g. be written to a file.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code based on template
!! * Later revisions - see repository log
!!
subroutine get_obs_pdaf(step, dim_obs, observation)

  use synobs_pdaf_mod, &
       only: write_syn_obs
  use assimilation_pdaf_mod, &
       only: file_synobs
  use parallel_pdaf_mod, &
       only: mype_filter

  implicit none

  ! *** Arguments
  integer, intent(in) :: step                 !< Current time step
  integer, intent(in) :: dim_obs              !< Dimension of obs. vector
  real, intent(in)    :: observation(dim_obs) !< Observation vector

  ! Local variables
  character(len=100) :: file_syntobs          ! Full name of file for synthetic observations
  character(len=4) :: procstr                 ! 4-digit string for process rank


! *************************
! *** store observation ***
! *************************

  ! Note 'file_synobs' is set in the observation module when
  ! initialize the output file

  call write_syn_obs(step, file_synobs, dim_obs, observation, 1)

end subroutine get_obs_pdaf

