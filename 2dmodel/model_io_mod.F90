!> Module for netcdf operations in 2D tutorial model
!!
!! This module provides functionality to read and write
!! netcdf files for the 2-dimensional tutorial model
!! with parallelization.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_io_mod

  use netcdf

  implicit none
  save
  private

! *** Variables specific for 2D tutorial model ***

  integer :: ncid                       !< ID of output file
  integer :: id_fieldA, id_fieldB       !< file IDs for both fields
  integer :: dimid_x, dimid_y, dimid_t  !< dimension IDs

  public io_write_sngl, io_read_sngl

contains

!-------------------------------------------------------------------------------
!> Write a field into a netCDF output file
!!
!! Routine to write a specific field at one time step.
!! Each call generates a single file.
!!
  subroutine io_write_sngl(step, filename, field_p)

    use mpi
    use model_parallel_mod, &
         only: mype_world, MPIErr, COMM_2Dmodel
    use model_mod, &
         only: nx_p, nx, ny, total_steps

    implicit none

    ! Arguments
    integer, intent(in) :: step                 !< Model time step
    character(len=100), intent(in) :: filename  !< Name of output file
    real, intent(in) :: field_p(:,:)            !< Decomposed model field

    ! Local variables
    integer :: id_field             !< Netcdf field if
    integer :: id_t                 !< Netcdf timestep id
    integer :: dimids(3)            !< Array for netcdf operation
    integer :: countv(3), startv(3) !< Vectors for NC operations
    character(len=6) :: fieldstr    !< String for field name
    real, allocatable :: field(:,:) !< Array for global model field


    ! *** Gather global field on process 0
    allocate(field(ny, nx))

    call MPI_Gather(field_p, nx_p*ny, MPI_DOUBLE_PRECISION, field, nx_p*ny, &
         MPI_DOUBLE_PRECISION, 0, COMM_2Dmodel, MPIerr)

    ! *** Write file on process 0

    if (mype_world==0) then

       ! *** Create file and define dimensions

       call nfcheck( NF90_CREATE(trim(filename), 0, ncid))
       call nfcheck( NF90_DEF_DIM(ncid, 'dim_x', nx, dimid_x))
       call nfcheck( NF90_DEF_DIM(ncid, 'dim_y', ny, dimid_y))
       call nfcheck( NF90_DEF_DIM(ncid, 'dim_steps', 1, dimid_t))

       ! *** Define variables

       ! Time step
       call nfcheck( NF90_DEF_VAR(ncid, 'timestep', NF90_INT, dimid_t, id_t))

       ! Field
       dimids(1) = dimid_y
       dimids(2) = dimid_x
       dimids(3) = dimid_t
       call nfcheck( NF90_DEF_VAR(ncid, 'field', NF90_DOUBLE, dimids(1:3), id_field))

       call nfcheck( NF90_ENDDEF(ncid)) 

       ! *** Write variables

       ! Write timestep
       call nfcheck( NF90_PUT_VAR(ncid, id_t, step))

       ! Write field
       startv(1) = 1
       countv(1) = ny
       startv(2) = 1
       countv(2) = nx
       startv(3) = 1
       countv(3) = 1
       call nfcheck( NF90_PUT_VAR(ncid, id_field, field, start=startv, count=countv))

       call nfcheck( NF90_CLOSE(ncid))

    end if

    deallocate(field)

  end subroutine io_write_sngl

!-------------------------------------------------------------------------------
!> Read field on subdomain from a netCDF output file
!!
!! Routine to read the subdomain-part of a specific field
!! at one time step.
!!
  subroutine io_read_sngl(filename, field_p)

    use mpi
    use model_parallel_mod, &
         only: mype_world, MPIErr, COMM_2Dmodel
    use model_mod, &
         only: nx_p, nx, ny, total_steps

    implicit none

    ! Arguments
    character(len=100), intent(in) :: filename  !< Name of output file
    real, intent(inout) :: field_p(:,:)         !< Decomposed model field

    ! Local variables
    integer :: id_field             ! Netcdf field id
    integer :: countv(3), startv(3) ! Vectors for NC operations
    integer :: off_nx               ! Offset of local grid in global domain in x-direction

    ! *** Read field

     ! offset in x-direction due to domain-decomposition
     off_nx = nx_p*mype_world

    ! Open file and get field ID

    call nfcheck( NF90_OPEN(trim(filename), NF90_NOWRITE, ncid))
    call nfcheck( NF90_INQ_VARID(ncid, 'field', id_field))

    startv(1) = 1
    countv(1) = ny
    startv(2) = 1+off_nx
    countv(2) = nx_p
    startv(3) = 1
    countv(3) = 1

    call nfcheck( NF90_GET_VAR(ncid, id_field, field_p, &
         start=startv(1:3), count=countv(1:3)))

    call nfcheck( NF90_CLOSE(ncid))

  end subroutine io_read_sngl

!-------------------------------------------------------------------------------
!> Check status of netcdf operation
!!
  subroutine nfcheck(status)

    use netcdf
    use model_parallel_mod, &
         only: abort_parallel

! *** Aruments ***
    integer, intent ( in) :: status   ! Reading status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       call abort_parallel()
    end if

  end subroutine nfcheck

end module model_io_mod
