!>  Main driver for PDAF offline program
!!
!! This is the main program for the PDAF assimilation
!! program in offline coupled mode.
!!
!! In the offline mode, the ensemble integrations are performed
!! in a separate model program and the forecasted ensemble can
!! be read from files. After initializing the ensemble information
!! by reading model outputs, a single analysis step is performed.
!! Subsequently, the analysis ensemble can be written to files
!! that can be used to initialize another ensemble forecast.
!!
!! Using a parallelization with domain-decomposition, the offline
!! mode can be used to perform assimilation with domain-decomposed
!! models. Alternatively one can just split each field in the state
!! vector into parts of approximately equal size.
!!
!! This file is generic and usually does not need any changes when
!! implementing the offline data assimilation for a specific model.
!!
!! __Revision history:__
!! * 2008-07 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
program main_offline

  use parallel_pdaf_mod, &   ! Parallelization
       only: nproc_ens, myproc_ens
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
  initscreen: if (myproc_ens == 0) then

     write (*, '(/8x, a/)') '+++++ PDAF offline mode +++++'
     write (*, '(9x, a)') 'Data assimilation with PDAF'

     if (nproc_ens > 1) then
        write (*, '(/8x, a, i4, a/)') 'Running on ', nproc_ens, ' MPI processes'
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

  if (myproc_ens == 0) &
       write (*, '(/2x, a)') 'PDAF offline mode: START ASSIMILATION'

  call assimilate_pdaf_offline()


! ********************
! *** Finishing up ***
! ********************

! *** Final screen output ***
  if (myproc_ens == 0) &
       write (*, '(/1x, a)') 'PDAF offline mode: EXITED ASSIMILATION'

  ! *** Finalize PDAF - print memory and timing information
  call finalize_pdaf(0)

end program MAIN_OFFLINE
