!>  Initialize model grid information
!!
!! Routine to perform initialization of model dimensions
!! for offline coupled mode of PDAF.
!!
!! Initialized are usually
!! * the size of the model grid in each direction
!! * the number of dimensions in the model grid
!! * the coordinates of the model grid points
!! * with parallelization: the process-local grid sizes
!!   and potentially the offset in the global grid
!!
!! The information is used in different call-back routines, e.g. to
!! determine the size of the state vector and for localization.
!!
!! Variant for using the tutorial 2D model in offline mode.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module initialize_grid_mod

contains

  subroutine initialize_grid()

    use parallel_pdaf_mod, &         ! Model parallelzation variables
         only: myproc_ens, nproc_ens, abort_parallel

    ! Specific for model
    use model_pdaf_mod, &            ! Model grid variables
         only: n_dim, nx, ny, nx_p, offset_x_p, coords_x_p, coords_y_p

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
    if (myproc_ens == 0) then
     write (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     write (*, '(5x,a)') 'MODEL: 2D-2fields tutorial model'
       write (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
    end if

! *** Initialize size of local nx for parallelization ***
    if (nproc_ens==1 .or. nproc_ens==2 .or. nproc_ens==3 .or. nproc_ens==4 .or. &
         nproc_ens==6 .or. nproc_ens==9 .or. nproc_ens==12 .or. nproc_ens==18) then
       ! Split x-direction in chunks of equal size
       nx_p = nx / nproc_ens
    else
       write (*,*) 'ERROR: Invalid number of processes'
       call abort_parallel()
    end if

    if (myproc_ens == 0 .and. nproc_ens > 1) then
       write (*, '(/2x, a, i3, a)') &
            '-- Domain decomposition over', nproc_ens, ' Processs'
       write (*, '(2x,a,i3,a,i3)') &
            '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny
    end if

    ! Set offset of process-local grid in global grid
    offset_x_p = nx_p*myproc_ens


! *************************************
! *** Initialize coordinates        ***
! *************************************

    ! The model coordinates are here the grid point indices
    ! stored as real values

    ! allocate memory for process-local part of fields
    allocate(coords_x_p(nx_p))
    allocate(coords_y_p(ny))

    ! Account for decomposition in x-direction
    do i = 1, nx_p
       coords_x_p(i) = real(i + offset_x_p)
    end do

    ! We don't use decomposition in y-direction
    do j = 1, ny
       coords_y_p(j) = real(j)
    end do

  end subroutine initialize_grid

end module initialize_grid_mod
