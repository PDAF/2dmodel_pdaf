!> Main program
!!
!! The program only serves to be able to compile
!! the PDAF online template routines for testing
!! their consistency.
!!
!! This variant is for the 
!!           flexible parallelization
!! variant of PDAF using PDAF3_assimilate routines. 
!! It shows the structure of the required outper loop
!! which enables to integrate an ensemble of model states.
!! 
!! In the online implementation with a real model
!! this driver program is replaced by the actual
!! model code.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code using PDAF3_assimilate
!!
program MAIN

  use mpi                      ! MPI
  use PDAF, &                  ! PDAF
       only: PDAF_get_fcst_info
  use parallel_pdaf_mod, &     ! Parallelization
       only: nproc_model, myproc_model
  use model_mod, &             ! Module provided by model code
       only: total_steps, n_dim, nx, ny, nx_p, offset_x_p, &
       time, dt

  implicit none

! local variables
  integer :: istep             ! Counter
  integer :: mpierr            ! MPI error flag
  integer :: mycomm, myproc, nproc  ! MPI variables for communicator, rank and size
  integer :: nsteps            ! Number of time steps to be performed in current forecast
  integer :: doexit            ! Whether to exit forecasting (1=true)
  real :: timenow              ! Model time
  

! *** Initialize MPI ***

  ! If the model itself is parallelized this step is done by the model

  call MPI_INIT(mpierr);
  call MPI_Comm_Size(MPI_COMM_WORLD, nproc, mpierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myproc, mpierr)

  mycomm = MPI_COMM_WORLD


  if (myproc==0) then
     write (*,*) '**********************************************************************'
     write (*,*) '*   THIS IS A TEST PROGRAM TO CHECK THE TEMPLATE CODE CONSISTENCY    *'
     write (*,*) '*                   Run this program with:                           *'
     write (*,*) '*   ./model_pdaf_flexible -dim_ens NENS                              *'
     write (*,*) '* with ensemble size NENS (=2 is good for testing, =1 does not work) *'
     write (*,*) '* to run a single model task. Or                                     *'
     write (*,*) '*  mpirun -np NP ./model_pdaf_flexible -dim_ens NENS -n_tasks NTSK   *'
     write (*,*) '* with 1<=NTSK<=NENS (NTSK=1 runs a single model task)               *'
     write (*,*) '* to run with ensemble size NENS and NTSK model tasks                *'
     write (*,*) '**********************************************************************'
  end if

  
! *** Initialize MPI communicators for PDAF (model, filter and coupling) ***

  ! This step is always inserted directly after the MPI initialization

  call init_parallel_pdaf(1, mycomm, myproc, nproc)

  
  ! MODEL: Here the model would perform its initialization
  ! MODEL: We just set the local grid dimension and offset

  nx_p = nx / nproc_model
  offset_x_p = nx_p * myproc_model

  write (*, '(/2x, a, i3, a)') &
       '-- Domain decomposition over', nproc_model, ' Processs'
  write (*, '(2x,a,i3,a,i3)') &
       '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny


! *** Initialize PDAF ***

  ! This step is always inserted after the model initialization
  ! is complete and just before the time stepping starts

  call init_pdaf()


! *** Ensemble forecasting and analysis steps ***

  pdaf_modelloop: do  

     ! *** Forecast ensemble state ***

     ! *** PDAF: Get forecast information ***
     call PDAF_get_fcst_info(nsteps, timenow, doexit)

     ! *** Check exit flag ***
     if (doexit==1) exit pdaf_modelloop
 
     ! Initialize current model time
     time = timenow

     ! *** run time stepper ***  

     ! MODEL: In the real model this is the time stepping loop of the model
     timesteps: do istep = 1, nsteps

        ! MODEL: Here the model code would compute the time stepping

        write (*,'(3x, a, f6.2)') 'main.F90: Do stepping, time', time

        ! The model would increment the time information
        time = time + dt  

        ! *** Let PDAF check forecast progress and perform analysis ***

        ! This step is inserted in the time stepping loop
        ! usually just before the 'end do'

        call assimilate_pdaf()

     end do timesteps
  enddo pdaf_modelloop


! *** Finalize PDAF - print memory and timing information ***

  ! This step can be inserted at the end of the model code
  ! before the MPI parallelization is finalized

  call finalize_pdaf(0)


! *** Terminate MPI

  ! If the model itself is parallelized this step is done by the model

  call  MPI_Finalize(MPIerr)

end program MAIN
