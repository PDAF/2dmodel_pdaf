!> Main program
!!
!! The program only serves to be able to compile
!! the PDAF online template routines for testing
!! their consistency.
!!
!! The program shows the setup for the fully-parallel
!! implementation variant of PDAF.
!! 
!! In the online implementation with a real model
!! this driver program is replaced by the actual
!! model code.
!!
!! __Revision history:__
!! * 2021-11 - Lars Nerger - Initial code
!!
program MAIN

  use mpi                      ! MPI
  use parallel_pdaf_mod, &     ! Parallelization
       only: nproc_model, myproc_model
  use model_mod, &             ! Module provided by model code
       only: total_steps, n_dim, nx, ny, nx_p, offset_x_p

  implicit none

! local variables
  integer :: istep       ! Counter
  integer :: mpierr      ! MPI error flag
  integer :: mycomm, myproc, nproc


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
     write (*,*) '*          mpirun -np NENS ./PDAF_online -dim_ens NENS               *'
     write (*,*) '* with ensemble size NENS (=2 is good for testing, =1 does not work) *'
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

  ! MODEL: In the real model this is the time stepping loop of the model
  timesteps: do istep = 1, total_steps

     ! MODEL: Here the model code would compute the time stepping


     ! *** Perform analysis ***

     ! This step is inserted in the time stepping loop
     ! usually just before the 'end do'

     call assimilate_pdaf()

  enddo timesteps


! *** Finalize PDAF - print memory and timing information ***

  ! This step can be inserted at the end of the model code
  ! before the MPI parallelization is finalized

  call finalize_pdaf(0)


! *** Terminate MPI

  ! If the model itself is parallelized this step is done by the model

  call  MPI_Finalize(MPIerr)

end program MAIN
