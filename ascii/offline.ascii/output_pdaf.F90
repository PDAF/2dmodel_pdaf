!>  File output for PDAF
!!
!! The routine output_pdaf gathers each field with MPI
!! and then calls the routine write_field_pdaf to do
!! the actual file writing.
!!
!! Implementation for the 2D example with domain decomposition
!!
!! Note: The file ouptut uses text files. While these are 
!! OK for a simple tutorial like this they are not recommended for
!! real high-dimensional cases. It is recommended to use a binary
!! file format, e.g netCDF, instead.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module output_pdaf_mod

contains
  subroutine write_pdaf(step, dim_p, dim_ens, state_p, ens_p, variance_p)

    use mpi                             ! MPI
    use model_pdaf_mod, &               ! Model variables
         only: nx, ny, nx_p
    use parallel_pdaf_mod, &            ! Parallelization variables
         only: COMM_filter, mype_filter, npes_filter, MPIerr, MPIstatus
    use statevector_pdaf_mod, &         ! Statevector variables
         only: id, sfields, n_fields

    implicit none

! *** Arguments ***
    integer, intent(in) :: step         !< Current time step (negative for call after forecast)
    integer, intent(in) :: dim_p        !< Process-local state dimension
    integer, intent(in) :: dim_ens      !< Size of state ensemble
    real, intent(inout) :: state_p(dim_p)        !< Process-local forecast/analysis state
    real, intent(inout) :: ens_p(dim_p, dim_ens) !< Process-local state ensemble
    real, intent(inout) :: variance_p(dim_p)     !< Process-local variance vector


! *** Local variables ***
    integer :: member                   ! Counters
    integer :: fid                      ! field index
    character(len=2) :: ensstr          ! String for ensemble member
    character(len=100) :: filestr       ! Part of file name specific for state or ensemble member
    real, allocatable :: fieldvec(:)    ! global state vector


! *******************
! *** File output ***
! *******************

    ! Allocate arrays
    if (mype_filter==0) then
       allocate(fieldvec(nx*ny))
    else
       allocate(fieldvec(1))
    end if


! *** Write ensemble state fields ***

    do member = 1, dim_ens
        
       ! Set string for ensemble members
       write (ensstr, '(i2.2)') member

       do fid = 1, n_fields
     
          ! Gather single field
          call MPI_Gather(ens_p(sfields(fid)%off+1, member), sfields(fid)%dim, MPI_DOUBLE_PRECISION, &
               fieldvec, nx_p*ny, MPI_DOUBLE_PRECISION, &
               0, COMM_filter, MPIerr)             

          ! Write gathered ensemble member field
          if (mype_filter==0) then
             filestr = 'ens'//trim(sfields(fid)%fname)//'_'//trim(ensstr)

             ! Call generic writing routine
             call write_field_pdaf(step, nx, ny, fieldvec, filestr) 
          end if
       end do
    end do


! *** Write ensemble mean state fields ***

    do fid = 1, n_fields

       ! Gather single field
       call MPI_Gather(state_p(sfields(fid)%off+1), sfields(fid)%dim, MPI_DOUBLE_PRECISION, &
            fieldvec, nx_p*ny, MPI_DOUBLE_PRECISION, &
            0, COMM_filter, MPIerr)             

       ! Write gathered state vector field
       if (mype_filter==0) then
     
          filestr = 'state'//trim(sfields(fid)%fname)

          ! Call generic writing routine
          call write_field_pdaf(step, nx, ny, fieldvec, filestr) 

       end if
    end do


! *** Write ensemble variance fields ***

    do fid = 1, n_fields

       ! Gather single field
       call MPI_Gather(variance_p(sfields(fid)%off+1), sfields(fid)%dim, MPI_DOUBLE_PRECISION, &
            fieldvec, nx_p*ny, MPI_DOUBLE_PRECISION, &
            0, COMM_filter, MPIerr)             

       ! Write gathered state vector field
       if (mype_filter==0) then
     
          filestr = 'variance'//trim(sfields(fid)%fname)

          ! Call generic writing routine
          call write_field_pdaf(step, nx, ny, fieldvec, filestr) 

       end if
    end do

! *** Clean up ***

    deallocate(fieldvec)

  end subroutine write_pdaf

!-------------------------------------------------------------------------------
!> Write a single field
!!
!! This routine write a field from the 2D tutorial model
!! It gets the field in form a a vector, e.g. a part of a
!! state vector and maps this onto a 2D array assuming 
!! a certain ordering consistent with e.g. init_ens_pdaf.
!!
  subroutine write_field_pdaf(step, nx, ny, fieldvec, filestr)

    implicit none

! *** Arguments
    integer, intent(in) :: step               !< Current time step (negative for call after forecast)
    integer, intent(in) :: nx                 !< Domain size in x-direction
    integer, intent(in) :: ny                 !< Domain size in y-direction
    real, intent(in)    :: fieldvec(nx*ny)    !< Vector holding a single model field
    character(len=100), intent(in) :: filestr !< Part of file name specific for state or ensemble member

! *** Local variables ***
    integer :: i, j                     ! Counters
    character(len=2) :: stepstr         ! String for time step
    character(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
    ! Variables for parallelization - global fields
    real, allocatable :: field(:,:)     ! global model field


! **********************
! *** INITIALIZATION ***
! **********************

    if (step==0) then
       anastr = 'ini'
    elseif (step<0) then
       anastr = 'for'
    else
       anastr = 'ana'
    end if

    ! Set string for time step
    if (step>=0) then
       write (stepstr, '(i2.2)') step
    else
       write (stepstr, '(i2.2)') -step
    end if

    ! Allocate arrays
    allocate(field(ny, nx))


! *******************
! *** File output ***
! *******************

    do j = 1, nx
       do i = 1, ny
          field(i, j) = fieldvec(i + (j-1)*ny)
       end do
    end do

    open(11, file = trim(filestr)//'_step'//trim(stepstr)//'_'//trim(anastr)//'.txt', status = 'replace')
 
    do i = 1, ny
       write (11, *) field(i, :)
    end do

    close(11)


! *******************
! *** Clean up    ***
! *******************

    deallocate(field)

  end subroutine write_field_pdaf

!-------------------------------------------------------------------------------
!> Check status of netcdf operation
!!
  subroutine nfcheck(status)

    use netcdf
    use parallel_pdaf_mod, &
         only: abort_parallel

! *** Aruments ***
    integer, intent ( in) :: status   ! Reading status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       call abort_parallel()
    end if

  end subroutine nfcheck

end module output_pdaf_mod
