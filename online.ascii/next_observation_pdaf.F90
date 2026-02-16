!>  Initialize information for next forecast phase
!!
!! User-supplied call-back routine for PDAF.
!!
!! The subroutine is called before at the initialization
!! of the forecastss by PDAF_init_forecast and then for
!! each subsequent forecasts phase by PDAF3_assimilate.
!! It has to initialize the number of time steps until
!! the next available observation (nsteps) and the
!! current model time (time). In addition the exit
!! flag (exit) has to be initialized. This indicates for
!! the case of the flexible parallelization variant, if
!! the data assimilation process is completed such that
!! the ensemble loop in the model routine can be exited.
!!
!! The routine is called by all processes
!!         
!! Version for the 2D tutorial model.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine next_observation_pdaf(stepnow, nsteps, doexit, time)

  use assimilation_pdaf_mod, &     ! Assimilation variables
       only: delt_obs
  use parallel_pdaf_mod, &         ! Parallelization variables
       only: mype_ens
  use model_pdaf_mod, &            ! Model variables
       only: total_steps

  implicit none

! *** Arguments ***
  integer, intent(in)  :: stepnow  !< Number of the current time step
  integer, intent(out) :: nsteps   !< Number of time steps until next obs
  integer, intent(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  real, intent(out)    :: time     !< Current model (physical) time


! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************

  time = 0.0          ! Not used in this implementation
  doexit = 0          ! Not used in this implementation

  if (stepnow + nsteps <= total_steps) then
     ! *** During the assimilation process ***
     nsteps = delt_obs   ! This assumes a constant time step interval

     if (mype_ens == 0) write (*, '(a, i7, 3x, a, i7)') &
          'model-PDAF', stepnow, 'Next observation at time step', stepnow + nsteps
  else
     ! *** End of assimilation process ***
     nsteps = 0          ! No more steps

     if (mype_ens == 0) write (*, '(a, i7, 3x, a)') &
          'model-PDAF', stepnow, 'No more observations - end assimilation'
  end if

end subroutine next_observation_pdaf
