!>  Initialize model
!!
!! Initialization routine for using the tutorial 2D model
!! in offline mode.
!!
!! The routine defines the size of the model grid
!! and the model coordinates.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module initialize_grid_mod

contains

  subroutine initialize_grid()

    use model_pdaf_mod, &            ! Model grid variables
         only: nx, ny, nx_p, n_dim, &
         coords_x_p, coords_y_p
    use parallel_pdaf_mod, &         ! Model parallelzation variables
         only: mype_world, npes_world, abort_parallel

    implicit none

! *** local variables ***
    integer :: i, j                 ! Counters


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***
    nx = 36          ! Extent of grid in x-direction
    ny = 18          ! Extent of grid in y-direction
    n_dim = 2        ! Number of model dimensions

! *** Screen output ***
    if (mype_world == 0) then
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     WRITE (*, '(22x,a)') 'MODEL: 2D-2fields tutorial model'
       write (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
    end if

! *** Initialize size of local nx for parallelization ***
    if (npes_world==1 .or. npes_world==2 .or. npes_world==3 .or. npes_world==4 .or. &
         npes_world==6 .or. npes_world==9 .or. npes_world==12 .or. npes_world==18) then
       ! Split x-direction in chunks of equal size
       nx_p = nx / npes_world
    else
       write (*,*) 'ERROR: Invalid number of processes'
       call abort_parallel()
    end if

    if (mype_world == 0 .and. npes_world > 1) then
       write (*, '(/2x, a, i3, a)') &
            '-- Domain decomposition over', npes_world, ' Processs'
       write (*, '(2x,a,i3,a,i3)') &
            '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny
    end if

    ! allocate memory for process-local part of fields
    allocate(coords_x_p(nx_p))
    allocate(coords_y_p(ny))


! *************************************
! *** Initialize coordinates        ***
! *************************************

    ! The model coordinates are the grid point indices
    ! stored as real values

    ! Account for decomposition in x-direction
    do i = 1, nx_p
       coords_x_p(i) = real(i + nx_p*mype_world)
    end do

    ! We don't use decomposition in y-direction
    do j = 1, ny
       coords_y_p(j) = real(j)
    end do

  end subroutine initialize_grid

end module initialize_grid_mod
