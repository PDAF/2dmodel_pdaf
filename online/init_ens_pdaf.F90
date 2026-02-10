!> Initialize ensemble
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all DA methods.
!!
!! The routine is called by PDAF when the DA is initialized.
!! It has to initialize an ensemble of dim_ens states.
!!
!! The routine is called by all filter processes and 
!! initializes the ensemble for the process-local domain.
!!
!! Implementation for the 2D online example with parallelization.
!! Here, the ensemble is directly read from files.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  use model_pdaf_mod, &           ! Model variables
       only: nx, ny, nx_p
  use parallel_pdaf_mod, &        ! Assimilation parallelization variables
       only: mype_filter, mype_model
  use statevector_pdaf_mod, &     ! State vector variables
       only: id, sfields, n_fields

  implicit none

! *** Arguments ***
  integer, intent(in) :: filtertype                !< Type of filter to initialize
  integer, intent(in) :: dim_p                     !< process-local state dimension
  integer, intent(in) :: dim_ens                   !< Size of ensemble
  real, intent(inout) :: state_p(dim_p)            !< process-local model state
  !< (It is not necessary to initialize the array 'state_p' for ensemble filters.
  !< It is available here only for convenience and can be used freely.)
  real, intent(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for ensemble filters
  real, intent(out)   :: ens_p(dim_p, dim_ens)     !< process-local state ensemble
  integer, intent(inout) :: flag                   !< PDAF status flag

! *** local variables ***
  integer :: i, j, s, fid, member     ! Counters
  real, allocatable :: field(:,:)     ! global model field
  character(len=2) :: ensstr          ! String for ensemble member


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-Process 0 ***
  if (mype_filter==0) then
     write (*, '(/9x, a)') 'Initialize state ensemble'
     write (*, '(9x, a)') '--- read ensemble from files'
     write (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  end if

  ! allocate memory for temporary fields
  allocate(field(ny, nx))


! ********************************
! *** Read ensemble from files ***
! ********************************

  do member = 1, dim_ens
     write (ensstr, '(i1)') member

     do fid = 1, n_fields

       ! Read field
        open(11, file = '../inputs_online_2fields/ens'//trim(sfields(fid)%fname)//'_'//trim(ensstr)//'.txt', &
             status='old')

        ! Read global field
        do i = 1, ny
           read (11, *) field(i, :)
        end do

        ! +++ Note on counter s:
        ! +++ Using the counter s looks primitive, but it
        ! +++ makes the code fail-save because it avoids
        ! +++ index calculations involving nx_p or ny.

        ! Initialize process-local part of ensemble
        s = sfields(fid)%off
        do j = 1, nx_p
           do i = 1, ny
              s = s + 1
              ens_p(s, member) = field(i, nx_p*mype_model + j)
           end do
        end do

        close(11)
     end do
  end do


! ****************
! *** clean up ***
! ****************

  deallocate(field)

end subroutine init_ens_pdaf
