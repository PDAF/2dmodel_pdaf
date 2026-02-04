!>  Time stepping loop of tutorial model
!!
!! Time stepping for simple 2D tutorial model
!! with domain decomposition.
!!
!! Each time step the field is shifted by one grid 
!! point in the vertical direction (first array index).
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_step_mod

contains

  subroutine stepping()

    use mpi                     ! MPI
    use model_mod, &            ! Model variables
         only: nx, ny, nx_p, fieldA_p, fieldB_p, total_steps
    use model_parallel_mod, &   ! Model parallelization
         only: mype_world, MPIErr, COMM_2Dmodel

    implicit none

! *** local variables ***
    integer :: step, i, j        ! Counters
    character(len=2) :: stepstr  ! String for time step
    real :: store                ! Store single field element
    real, allocatable :: field(:,:) ! Global model field


! ****************
! *** STEPPING ***
! ****************

    if (mype_world==0) write (*, '(1x, a)') 'MODEL INTEGRATION'

    steps: do step = 1 , total_steps

       if (mype_world==0) write (*,*) 'step', step

! *** Time step: Shift fields vertically ***
       do j = 1, nx_p
          ! Field A
          store = fieldA_p(ny, j)

          do i = ny-1,1,-1
             fieldA_p(i+1, j) = fieldA_p(i, j)
          end do

          fieldA_p(1, j) = store

          ! Field B
          store = fieldB_p(ny, j)

          do i = ny-1,1,-1
             fieldB_p(i+1, j) = fieldB_p(i, j)
          end do

          fieldB_p(1, j) = store
       end do

! *** Write new fields into files ***

       ! Gather global field on process 0
       allocate(field(ny, nx))

       call MPI_Gather(fieldA_p, nx_p*ny, MPI_DOUBLE_PRECISION, field, nx_p*ny, &
            MPI_DOUBLE_PRECISION, 0, COMM_2Dmodel, MPIerr)

       ! Write file from process 0
       if (mype_world==0) then
          write (stepstr, '(i2.2)') step
          open(11, file = 'trueA_step'//trim(stepstr)//'.txt', status = 'replace')

          do i = 1, ny
             write (11, *) field(i, :)
          end do

          close(11)     
       end if

       call MPI_Gather(fieldB_p, nx_p*ny, MPI_DOUBLE_PRECISION, field, nx_p*ny, &
            MPI_DOUBLE_PRECISION, 0, COMM_2Dmodel, MPIerr)

       ! Write file from process 0
       if (mype_world==0) then
          write (stepstr, '(i2.2)') step
          open(11, file = 'trueB_step'//trim(stepstr)//'.txt', status = 'replace')

          do i = 1, ny
             write (11, *) field(i, :)
          end do

          close(11)     
       end if

       deallocate(field)

    end do steps

  end subroutine stepping

end module model_step_mod
