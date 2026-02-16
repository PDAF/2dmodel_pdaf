!>  Initialize model
!!
!! Initialization routine for the simple 2D model with
!! parallelization of the model.
!!
!! The routine defines the size of the model grid and
!! reads the initial state from a file. 
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_init_mod

contains

  subroutine initialize()

    use model_mod, &                ! Model variables
         only: nx, ny, nx_p, n_dim, fieldA_p, fieldB_p, &
         coords_x_p, coords_y_p, total_steps
    use model_parallel_mod, &       ! Model parallelzation variables
         only: mype_world, npes_world, mype_2Dmodel, npes_2Dmodel, abort_parallel
    use model_io_mod, &             ! File operations
         only: io_read_sngl

    implicit none

! *** local variables ***
    integer :: i, j                 ! Counters
    character(len=2) :: stepstr     ! String for time step
    character(len=100) :: filename  ! Name of output file


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***
    nx = 36          ! Extent of grid in x-direction
    ny = 18          ! Extent of grid in y-direction
    n_dim = 2        ! Number of model dimensions
    total_steps = 18 ! Number of time steps to perform

! *** Screen output ***
    if (mype_world == 0) then
       write (*, '(1x, a)') 'INITIALIZE PARALLELIZED 2D TUTORIAL MODEL'
       write (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
       write (*, '(10x,a,i4)') 'Time steps', total_steps
    end if

! *** Initialize size of local nx for parallelization ***
    if (npes_2Dmodel==1 .or. npes_2Dmodel==2 .or. npes_2Dmodel==3 .or. npes_2Dmodel==4 .or. &
         npes_2Dmodel==6 .or. npes_2Dmodel==9 .or. npes_2Dmodel==12 .or. npes_2Dmodel==18) then
       ! Split x-direction in chunks of equal size
       nx_p = nx / npes_2Dmodel
    else
       write (*,*) 'ERROR: Invalid number of processes'
       call abort_parallel()
    end if

    if (mype_world == 0 .and. npes_2Dmodel > 1) then
       write (*, '(/2x, a, i3, a)') &
            '-- Domain decomposition over', npes_2Dmodel, ' Processs'
       write (*, '(2x,a,i3,a,i3)') &
            '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny
    end if

    ! allocate memory for process-local part of fields
    allocate(fieldA_p(ny, nx_p))
    allocate(fieldB_p(ny, nx_p))
    allocate(coords_x_p(nx_p))
    allocate(coords_y_p(ny))


! *************************************
! *** Read initial fields from file ***
! *************************************

    filename = '../inputs_2fields/trueA_ini.nc'
    call io_read_sngl(filename, fieldA_p)

    filename = '../inputs_2fields/trueB_ini.nc'
    call io_read_sngl(filename, fieldB_p)


! *************************************
! *** Initialize coordinates        ***
! *************************************

    ! The model coordinates are the grid point indices
    ! stored as real values

    ! Account for decomposition in x-direction
    do i = 1, nx_p
       coords_x_p(i) = real(i + nx_p*mype_2Dmodel)
    end do

    ! We don't use decomposition in y-direction
    do j = 1, ny
       coords_y_p(j) = real(j)
    end do

#ifdef USE_PDAF
    ! Initialize PDAF
    call init_pdaf()
#endif

  end subroutine initialize

end module model_init_mod
