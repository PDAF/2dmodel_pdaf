!>  Main driver for PDAF offline tutorial
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

  use mpi                    ! MPI
  use parallel_pdaf_mod, &   ! Parallelization
       only: MPIerr, npes_world, mype_world, &
       init_parallel, finalize_parallel
  use initialize_grid_mod, &
       only: initialize_grid

  implicit none


! **********************
! *** Initialize MPI ***
! **********************

  call init_parallel() ! initializes MPI


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Initial Screen output ***
  initscreen: if (mype_world == 0) then

     write (*, '(/8x, a/)') '+++++ PDAF offline mode +++++'
     write (*, '(9x, a)') 'Data assimilation with PDAF'

     if (npes_world > 1) then
        write (*, '(/21x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     else
        write (*, '(/21x, a/)') 'Running on 1 PE'
     end if
     write (*, '(/)')
     
  end if initscreen

  
! *** Initialize MPI communicators for PDAF (model and filter) ***
! *** NOTE: It is always n_modeltasks=1 for offline mode       ***

  call init_parallel_pdaf_offline(1)

! *** Initialize model information ***
! *** This should only be information on the model grid

  call initialize_grid()


! *******************************
! ***      ASSIMILATION       ***
! *******************************

  ! *** Initialize PDAF ***

  call init_pdaf_offline()


  ! *** Perform analysis ***

  if (mype_world == 0) &
       write (*, '(/2x, a)') 'PDAF offline mode: START ASSIMILATION'

  call assimilate_pdaf_offline()


  ! Synchronize at barrier for exit
  call MPI_Barrier(MPI_COMM_WORLD, MPIerr) 


! ********************
! *** Finishing up ***
! ********************

! *** Final screen output ***
  if (mype_world == 0) &
       write (*, '(/1x, a)') 'PDAF offline mode: EXITED ASSIMILATION'

  ! *** Finalize PDAF - print memory and timing information
  call finalize_pdaf(0)

  ! End parallelization
  call finalize_parallel()

end program MAIN_OFFLINE
