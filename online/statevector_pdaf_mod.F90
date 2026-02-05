!> Building the state vector
!!
!! This module provides variables & routines for
!! defining the state vector.
!!
!! The module contains three routines
!! - **init_id** - initialize the array `id`
!! - **init_sfields** - initialize the array `sfields`
!! - **setup_statevector** - generic routine controlling the initialization
!!
!! The declarations of **id** and **sfields** as well as the
!! routines **init_id** and **init_sfields** might need to be
!! adapted to a particular modeling case. However, for most
!! parts also the configruation using the namelist is possible.
!!
!! __Revision history__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module statevector_pdaf_mod

  implicit none
  save

! *** Variables to handle multiple fields in the state vector ***

  !< Fortran type holding the indices of model fields in the state vector
  !< This can be extended to any number of fields - it severs to give each field a name
  type field_ids
     integer :: fieldA 
     INTEGER :: fieldB
  end type field_ids

  !< Fortran type storing size and offset of each model field in the state vector
  !< This is generic, but one could extend this type to more variables
  type state_field
     integer :: dim               ! size of field in state vector
     integer :: off               ! offset of field in state vector
     integer :: ndims             ! Number of dimensions
     character(len=10) :: name    ! Name of field variable (optional)
     character(len=10) :: fname   ! Name in output file (optional)     
  end type state_field

  !---- The next variables usually do not need editing -----

  !< Type variable holding field IDs in state vector
  type(field_ids) :: id

  !< number of fields in state vector
  integer :: n_fields                   

  !< Vector of type variable holding dimension and offset of each field
  type(state_field), allocatable :: sfields(:)

contains


! -----------------------------------------------------------------
!> This routine initializes the array `id`
!!
  subroutine init_id(n_fields)

    implicit none

! *** Arguments ***
    integer, intent(out) :: n_fields


! Set total number of fields
    n_fields = 2

! Set field IDs
    id%fieldA = 1
    id%fieldB = 2

  end subroutine init_id

! -----------------------------------------------------------------
!> This initializes the array sfields
!!
!! This routine initializes the sfields array with specifications
!! of the fields in the state vector.
!!
  subroutine init_sfields()

    use model_pdaf_mod, &
         only: nx_p, ny

    implicit none

! *** Local variables ***
    integer :: i           ! Counter


! *** Allocate ***

    allocate(sfields(n_fields))

! *****************************************************
! *** Specify sfields entry for each field variable ***
! *****************************************************

    ! fieldA
    sfields(id%fieldA)%ndims = 2
    sfields(id%fieldA)%name = 'fieldA'
    sfields(id%fieldA)%fname = 'A'

    ! fieldB
    sfields(id%fieldB)%ndims = 2
    sfields(id%fieldB)%name = 'fieldB'
    sfields(id%fieldB)%fname = 'B'


! **************************************
! ***   Set dimensions and offsets   ***
! **************************************

    ! Set field dimensions
    do i = 1, n_fields
       if (sfields(i)%ndims == 2) then
          sfields(i)%dim = nx_p * ny
       end if
    end do

    ! Define field offsets in state vector
    sfields(1)%off = 0
    do i = 2, n_fields
       sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
    end do

  end subroutine init_sfields


! -----------------------------------------------------------------
!> Initialize the state vector
!!
!! This routine is generic. Case-specific adaptions should only
!! by done in the routines init_id and init_sfields.
!!
  subroutine setup_statevector(dim_state, dim_state_p, screen)

    use parallel_pdaf_mod, &
         only: mype_ens, npes_ens, task_id, comm_ensemble, &
         comm_model, MPI_SUM, MPI_INTEGER, MPIerr

    implicit none

! *** Arguments ***
    integer, intent(out) :: dim_state    !< Global dimension of state vector
    integer, intent(out) :: dim_state_p  !< Local dimension of state vector
    integer, intent(in)  :: screen       !< Verbosity flag

! *** Local variables ***
    integer :: i                 ! Counters


! ***********************************
! *** Initialize the state vector ***
! ***********************************

! *** Initialize array `id` ***

    call init_id(n_fields)

! *** Initialize array `sfields` ***

    call init_sfields()

! *** Set state vector dimension ***

    dim_state_p = sum(sfields(:)%dim)

! *** Write information about the state vector ***

    if (mype_ens==0) then
       write (*,'(/a,2x,a)') 'model-PDAF', '*** Setup of state vector ***'
       write (*,'(a,5x,a,i5)') 'model-PDAF', '--- Number of fields in state vector:', n_fields
       write (*,'(a,a4,3x,a2,2x,a8,6x,a3,7x,a6)') &
            'model-PDAF','PE','ID', 'variable', 'dim', 'offset'
    end if

    if ((mype_ens==0 .and. screen<=2) .or. (task_id==1 .and. screen>2)) then
       do i = 1, n_fields
          write (*,'(a, i4, i5,3x,a10,2x,i3,2x,i10,3x,i10,4x,l,4x,l,2x,i4)') 'model-PDAF', &
               mype_ens, i, sfields(i)%name, sfields(i)%dim, sfields(i)%off
       end do
    end if

    if (npes_ens==1) then
       write (*,'(a,2x,a,1x,i10)') 'model-PDAF', 'Full state dimension: ',dim_state_p
       dim_state = dim_state_p
    else
       if (task_id==1) then
          if (screen>2 .or. mype_ens==0) &
               write (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
               'model-PDAF', 'PE', mype_ens, 'process-local full state dimension: ',dim_state_p

          call MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)
          if (mype_ens==0) then
             write (*,'(a,2x,a,1x,i10)') 'model-PDAF', 'Global state dimension: ',dim_state
          end if
       end if
    end if
    call MPI_Barrier(comm_ensemble, MPIerr)

  end subroutine setup_statevector

end module statevector_pdaf_mod
