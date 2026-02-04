!>  Main driver for PDAF tutorial
!!
!! This is a simple model program to demonstrate the
!! fully-parallel implementation of the online mode of PDAF. 
!!
!! The simple model has a 2-dimensional grid. The initial state
!! is read from a file. The time stepping consists in shifting
!! the field vertically (in the direction of the first array index)
!! by one grid point per time step. A period boundary condition is
!! applied by inserting the field from the upper boundary into the
!! lower one. 
!!
!! __Revision history:__
!! * 2013-09 - Lars Nerger - Initial code based on dummy model example
!! * Later revisions - see repository log
!!
program main

  use model_parallel_mod, &          ! Model parallelization
       only: init_parallel, finalize_parallel, mype_world
  use model_init_mod, &              ! Model initialization
       only: initialize
  use model_step_mod, &              ! Model integration
       only: stepping

  implicit none

! ********************************
! ***      INITIALIZATION      ***
! ********************************

  ! *** Initialize parallelization ***
  call init_parallel()

  ! *** Initial Screen output ***
  if (mype_world==0) then
     write (*, '(/10x, a/)') '+++++ PDAF tutorial - online mode +++++'
     write (*, '(10x, a)')   '   Decomposed 2D model with 2 fields'
     write (*, '(/)')
  end if

  ! *** Initialize model ***
  call initialize()


! *****************************
! ***      Integration      ***
! *****************************

  ! *** Perform integration ***
  call stepping()


! **************************
! ***      Clean up      ***
! **************************

#ifdef USE_PDAF
  ! End parallelization
  call finalize_pdaf()
#endif

  call finalize_parallel()

end program main
