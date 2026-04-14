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
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
module initialize_grid_mod

contains

  subroutine initialize_grid()

    use parallel_pdaf_mod, &         ! Model parallelzation variables
         only: myproc_ens, nproc_ens

    ! Specific variables for the model
    use model_pdaf_mod, &            ! Model grid variables
         only: n_dim, nx, ny, nx_p, offset_x_p !, coords_x_p, coords_y_p

    implicit none

! *** local variables ***


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model grid specifications ***

    ! Template reminder - delete when implementing functionality
    write (*,*) 'TEMPLATE initialize_grid_mod.F90: Initialize grid dimensions here!'

! +++ Dummy initialization to ensure that the code can be executed.
! +++ One can adapt the names of the variables defining the grid size
    nx = 12          ! Extent of grid in x-direction
    ny = 10          ! Extent of grid in y-direction
    n_dim = 2        ! Number of model dimensions


! *** Screen output ***
    if (myproc_ens == 0) then
     write (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     write (*, '(5x,a)') 'MODEL: +++ TEMPLATE - replace me +++ '
! +++ Adapt this output if other dimension variables are used
     write (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
    end if

! *** Initialize process-local grid size for parallelization ***

! +++ In case of parallelization one usually want to also
! +++ initialize the process-local grid sizes and possible
! +++ their offset in the global grid. Here, we only use
! +++ a variable nx_p for a decomposition in x-direction.
! +++ For a more sophisticated distribution one would usualy
! +++ distribute the grid in x and y direction.

! +++ Dummy initialization to ensure that the code can be executed
    ! Split x-direction in chunks of equal size
    ! This only works correctly if nx is a multiple of nproc_ens
    nx_p = nx / nproc_ens

    ! Set offset of sub-domain in global grid
    offset_x_p = nx_p*myproc_ens

    if (myproc_ens == 0 .and. nproc_ens > 1) then
       write (*, '(/2x, a, i3, a)') &
            '-- Domain decomposition over', nproc_ens, ' Processs'
! +++ Adapt this output if other dimension variables are used
       write (*, '(2x,a,i3,a,i3)') &
            '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny
    end if


! *************************************
! *** Initialize coordinates        ***
! *************************************

! +++ Localization requires the coordinates of the grid points
! +++ These can be initialized here to be used later.

    ! Template reminder - delete when implementing functionality
    write (*,*) 'TEMPLATE initialize_grid_mod.F90: Initialize model grid coordinates here!'

    ! allocate memory for process-local part of fields
!    allocate(coords_x_p(...))
!    allocate(coords_y_p(...))

!    coords_x_p = ?
!    coords_y_p = ?

  end subroutine initialize_grid

end module initialize_grid_mod
