!!
!! This is the main program for an example implementation of
!! PDAF with domain-decomposition and offline configuration.
!!
!! In the offline mode, we assume that the ensemble
!! integrations are performed in a separate program (model)
!! and the forecasted ensemble can be read from files. After
!! initializing the ensemble information by reading model
!! outputs, a single analysis step is performed. Subsequently,
!! the analysis ensemble can be written to files that can be 
!! used to initialize another ensemble forecast.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
program main_offline

  use parallel_pdaf_mod, &   ! Parallelization
       only: npes_ens, mype_ens, &
       init_parallel, finalize_parallel
  use initialize_grid_mod, &
       only: initialize_grid

  implicit none


! *************************************************
! *** Initialize MPI and communicators for PDAF ***
! *************************************************

  call init_parallel_pdaf_offline(1)


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Initial Screen output ***
  initscreen: if (mype_ens == 0) then

     write (*, '(/8x, a/)') '+++++ PDAF offline mode +++++'
     write (*, '(9x, a)') 'Data assimilation with PDAF'

     if (npes_ens > 1) then
        write (*, '(/8x, a, i4, a/)') 'Running on ', npes_ens, ' MPI processes'
     else
        write (*, '(/12x, a/)') 'Running on 1 process'
     end if
     write (*, '(/)')
     
  end if initscreen


! *** Initialize model information ***
! *** This should only be information on the model grid

  call initialize_grid()


! *******************************
! ***      ASSIMILATION       ***
! *******************************

  ! *** Initialize PDAF ***

  call init_pdaf_offline()


  ! *** Perform analysis ***

  if (mype_ens == 0) &
       write (*, '(/2x, a)') 'PDAF offline mode: START ASSIMILATION'

  call assimilate_pdaf_offline()

! ********************
! *** Finishing up ***
! ********************

! *** Final screen output ***
  if (mype_ens == 0) &
       write (*, '(/1x, a)') 'PDAF offline mode: EXITED ASSIMILATION'

  ! *** Finalize PDAF - print memory and timing information
  call finalize_pdaf(0)

end program MAIN_OFFLINE
