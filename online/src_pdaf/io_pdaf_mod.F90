!>  File reading and output for PDAF
!!
!! The routine output_pdaf gathers each field with MPI
!! and then calls the routine write_field_pdaf to do
!! the actual file writing.
!! For reading the field of a sub-domain is read.
!!
!! Implementation for the 2D example with domain decomposition
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module io_pdaf_mod

  logical :: write_state = .true.
  logical :: write_ens   = .false.
  logical :: write_var   = .false.

contains

!-------------------------------------------------------------------------------
!> Read field on subdomain from a netCDF output file into state vector
!!
!! Routine to read the subdomain-part of a specific field
!! at one time step.
!!
  subroutine read_state_field_pdaf(filename, state_p)

    use netcdf
    use parallel_pdaf_mod, &
         only: mype_model
    use model_pdaf_mod, &
         only: nx_p, ny

    implicit none

    ! Arguments
    character(len=100), intent(in) :: filename  !< Name of output file
    real, intent(inout) :: state_p(:)           !< Part of process-local state vector

    ! Local variables
    integer :: ncid                 ! ID of output file
    integer :: id_field             ! Netcdf field id
    integer :: countv(3), startv(3) ! Vectors for NC operations
    integer :: off_nx               ! Offset of local grid in global domain in x-direction

    ! *** Read field

     ! offset in x-direction due to domain-decomposition
     off_nx = nx_p*mype_model

    ! Open file and get field ID

    call nfcheck( NF90_OPEN(trim(filename), NF90_NOWRITE, ncid))
    call nfcheck( NF90_INQ_VARID(ncid, 'field', id_field))

    startv(1) = 1
    countv(1) = ny
    startv(2) = 1+off_nx
    countv(2) = nx_p
    startv(3) = 1
    countv(3) = 1

    call nfcheck( NF90_GET_VAR(ncid, id_field, state_p, &
         start=startv(1:3), count=countv(1:3)))

    call nfcheck( NF90_CLOSE(ncid))

  end subroutine read_state_field_pdaf

!-------------------------------------------------------------------------------
!> Read covariance matrix information
!!
!! Routine to read the subdomain-part of the EOFs
!! from a covariance matrix file and the related
!! singular values.
!!
  subroutine read_covar_pdaf(filename, dim_ens, eofs_p, svals, state_p)

    use netcdf
    use parallel_pdaf_mod, &
         only: mype_model, abort_parallel
    use model_pdaf_mod, &
         only: nx_p, ny
    use statevector_pdaf_mod, &     ! State vector variables
         only: sfields, n_fields

    implicit none

    ! Arguments
    character(len=100), intent(in) :: filename  !< Name of output file
    integer, intent(in) :: dim_ens
    real, intent(out) :: eofs_p(:,:)            !< Process-local EOF field
    real, intent(out) :: svals(:)               !< Singular values
    real, intent(out) :: state_p(:)             !< Process local mean state

    ! Local variables
    integer :: j, fid                   ! Counters
    integer :: ncid                     ! ID of output file
    integer :: rank_file                ! Number of EOFs stored in covariance file
    integer :: id_svals, id_eofs        ! IDs for fields
    integer :: id_state                 ! ID for field
    integer :: id_dim                   ! ID for dimension
    integer :: countv(3), startv(3)     ! Vectors for NC operations
    integer :: off_nx                   ! Offset of local grid in global domain in x-direction


    ! offset in x-direction due to domain-decomposition
    off_nx = nx_p*mype_model

    call nfcheck( NF90_OPEN(filename, NF90_NOWRITE, ncid))

    ! Read rank stored in file
    call nfcheck( NF90_INQ_DIMID(ncid, 'rank', id_dim))
    call nfcheck( NF90_Inquire_dimension(ncid, id_dim, len=rank_file))

    ! Check consistency of dimensions
    checkdim: if (rank_file >= dim_ens-1) then

       ! *** Read singular values

       ! Inquire ID
       call nfcheck( NF90_INQ_VARID(ncid, 'sigma', id_svals))

       ! Read array
       startv(1) = 1
       countv(1) = dim_ens-1

       call nfcheck( NF90_GET_VAR(ncid, id_svals, svals, start=startv(1:1), count=countv(1:1)))

       ! Read mean state
       do fid = 1, n_fields
          call nfcheck( NF90_INQ_VARID(ncid, 'mean'//trim(sfields(fid)%fname), id_state))

          startv(2) = 1+off_nx
          countv(2) = nx_p
          startv(1) = 1
          countv(1) = ny

          call nfcheck( NF90_GET_VAR(ncid, id_state, &
               state_p(sfields(fid)%off+1 : sfields(fid)%off+sfields(fid)%dim), &
               start=startv(1:2), count=countv(1:2)))
       end do

       ! Read eofs
       do j = 1, dim_ens-1
          do fid = 1, n_fields
             call nfcheck( NF90_INQ_VARID(ncid, 'u_svd'//trim(sfields(fid)%fname), id_eofs))

             startv(3) = j
             countv(3) = 1
             startv(2) = 1+off_nx
             countv(2) = nx_p
             startv(1) = 1
             countv(1) = ny

             call nfcheck( NF90_GET_VAR(ncid, id_eofs, &
                  eofs_p(sfields(fid)%off+1 : sfields(fid)%off+sfields(fid)%dim, j), &
                  start=startv(1:3), count=countv(1:3)))
          end do
        end do

        call nfcheck( NF90_CLOSE(ncid))

     else
        ! *** Rank stored in file is smaller than requested EOF rank ***
        write(*,*) 'Rank stored in file is smaller than requested EOF rank'

        call nfcheck( NF90_CLOSE(ncid))
        call abort_parallel()

     end if checkdim

  end subroutine read_covar_pdaf

!-------------------------------------------------------------------------------
!> Read read mean state from covariance matrix file
!!
!! Routine to read the mean state variable from the 
!! file holding the covariance matrix.
!!
  subroutine read_mean_covar_pdaf(filename, state_p)

    use netcdf
    use parallel_pdaf_mod, &
         only: mype_model
    use model_pdaf_mod, &
         only: nx_p, ny
    use statevector_pdaf_mod, &     ! State vector variables
         only: sfields, n_fields

    implicit none

    ! Arguments
    character(len=100), intent(in) :: filename  !< Name of output file
    real, intent(out) :: state_p(:)             !< Process local mean state

    ! Local variables
    integer :: fid                      ! Counters
    integer :: ncid                     ! ID of output file
    integer :: id_state                 ! ID for field
    integer :: countv(2), startv(2)     ! Vectors for NC operations
    integer :: off_nx                   ! Offset of local grid in global domain in x-direction


    ! offset in x-direction due to domain-decomposition
    off_nx = nx_p*mype_model

    call nfcheck( NF90_OPEN(filename, NF90_NOWRITE, ncid))

    ! Read mean state
    do fid = 1, n_fields
       call nfcheck( NF90_INQ_VARID(ncid, 'mean'//trim(sfields(fid)%fname), id_state))

       startv(2) = 1+off_nx
       countv(2) = nx_p
       startv(1) = 1
       countv(1) = ny

       call nfcheck( NF90_GET_VAR(ncid, id_state, &
            state_p(sfields(fid)%off+1 : sfields(fid)%off+sfields(fid)%dim), &
            start=startv(1:2), count=countv(1:2)))
    end do

    call nfcheck( NF90_CLOSE(ncid))

  end subroutine read_mean_covar_pdaf



!-------------------------------------------------------------------------------
!> Write fields from state vector
!!
!! Routine collecting the writing of the ensemble mean state, 
!! the ensemble states, and the variance fields. This routine
!! sets the file name and then calls the actual output routine.
!!
  subroutine write_pdaf(step, dim_p, dim_ens, state_p, ens_p, variance_p)

    use mpi                             ! MPI
    use statevector_pdaf_mod, &         ! Statevector variables
         only: sfields, n_fields

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
    character(len=2) :: stepstr         ! String for time step
    character(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
    character(len=2) :: ensstr          ! String for ensemble member
    character(len=100) :: filestr       ! Part of file name specific for state or ensemble member
    character(len=100) :: filename      ! Name of output file


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


! *******************
! *** File output ***
! *******************

! *** Write ensemble mean state fields ***

    if (write_state) then
       do fid = 1, n_fields
          ! Set file name
          filestr = 'state'//trim(sfields(fid)%fname)
          filename = trim(filestr)//'_step'//trim(stepstr)//'_'//trim(anastr)//'.nc'

          ! Call generic writing routine
          call write_field_pdaf(step, filename, &
               state_p(sfields(fid)%off+1: sfields(fid)%off+sfields(fid)%dim)) 
       end do
    end if

! *** Write ensemble state fields ***

    if (write_ens) then
       do member = 1, dim_ens
        
          ! Set string for ensemble members
          write (ensstr, '(i2.2)') member

          do fid = 1, n_fields
             ! Set file name
             filestr = 'ens'//trim(sfields(fid)%fname)//'_'//trim(ensstr)
             filename = trim(filestr)//'_step'//trim(stepstr)//'_'//trim(anastr)//'.nc'

             ! Call generic writing routine
             call write_field_pdaf(step, filename, &
                  ens_p(sfields(fid)%off+1:sfields(fid)%off+sfields(fid)%dim, member)) 
          end do
       end do
    end if

! *** Write ensemble variance fields ***

    if (write_var) then
       do fid = 1, n_fields
          ! Set file name
          filestr = 'variance'//trim(sfields(fid)%fname)
          filename = trim(filestr)//'_step'//trim(stepstr)//'_'//trim(anastr)//'.nc'

          ! Call generic writing routine
          call write_field_pdaf(step, filename, &
               variance_p(sfields(fid)%off+1 : sfields(fid)%off+sfields(fid)%dim)) 
       end do
    end if

  end subroutine write_pdaf


!-------------------------------------------------------------------------------
!> Write a single field
!!
!! This routine write a field from the 2D tutorial model
!! It gets the field in form a a vector, e.g. a part of a
!! state vector and write this into a 2D array assuming 
!! a certain ordering consistent with e.g. init_ens_pdaf.
!!
  subroutine write_field_pdaf(step, filename, fieldvec_p)

    use mpi
    use netcdf
    use parallel_pdaf_mod, &             ! Parallelization variables
         only: COMM_assim, mype_assim
    use model_pdaf_mod, &                ! Model variables
         only: nx, ny, nx_p

    implicit none

! *** Arguments
    integer, intent(in) :: step                !< Current time step (negative for call after forecast)
    real, intent(in)    :: fieldvec_p(:)       !< Vector holding a single model field
    character(len=100), intent(in) :: filename !< Part of file name specific for state or ensemble member

! *** Local variables ***
    integer :: ncid                       ! ID of output file
    integer :: MPIerr                     ! Error flag for MPI
    integer :: id_field                   ! Netcdf field if
    integer :: id_t                       ! Netcdf timestep id
    integer :: dimid_x, dimid_y, dimid_t  ! dimension IDs
    integer :: dimids(3)                  ! Array for netcdf operation
    integer :: countv(3), startv(3)       ! Vectors for NC operations
    real, allocatable :: field(:,:)       ! Array for global model field


! **************************************
! *** Use MPI to obtain global field ***
! **************************************

    if (mype_assim==0) then
       allocate(field(ny, nx))
    else
       allocate(field(1, 1))
    end if

    ! Gather global field
    call MPI_Gather(fieldvec_p, nx_p*ny, MPI_DOUBLE_PRECISION, &
         field, nx_p*ny, MPI_DOUBLE_PRECISION, &
         0, COMM_assim, MPIerr)             


! *******************
! *** File output ***
! *******************

    if (mype_assim==0) then

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

    endif

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

end module io_pdaf_mod
