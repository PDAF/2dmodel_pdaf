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
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine next_observation_pdaf(stepnow, nsteps, doexit, time)

  use assim_pdaf_mod, &            ! Assimilation variables
       only: delt_obs
  use parallel_pdaf_mod, &         ! Parallelization variables
       only: myproc_ens

  ! Specific for 2D tutorial model
  use model_pdaf_mod, &            ! Model variables
       only: total_steps, time_model => time

  implicit none

! *** Arguments ***
  integer, intent(in)  :: stepnow  !< Number of the current time step
  integer, intent(out) :: nsteps   !< Number of time steps until next obs
  integer, intent(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  real, intent(out)    :: time     !< Current model (physical) time


! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE next_observation_pdaf.F90: Set number of time steps in forecast!'

  nsteps = delt_obs   ! This assumes a constant time step interval

  if (stepnow+nsteps > total_steps) then
     nsteps = total_steps - stepnow
  end if

  if (stepnow == total_steps) then
     nsteps = 0
     doexit = 1
  end if


! *********************************
! *** Set current physical time ***
! *********************************

  time = time_model


! *********************
! *** Set exit flag ***
! *********************

  ! This is not used on the fully parallel mode
!   doexit = ??

end subroutine next_observation_pdaf
