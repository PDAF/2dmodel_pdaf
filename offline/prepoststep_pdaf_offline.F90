!>  Used-defined Pre/Poststep routine for PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all DA methods.
!! 
!! The routine is called before and after the analysis step.
!!
!! The routine provides full access to the state 
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep 
!! operations can be performed here. For example 
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by 
!! computing the estimated variances.
!!
!! If a user considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!!
!! Implementation for the 2D example with domain decomposition
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine prepoststep_pdaf_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  use mpi                             ! MPI
  use PDAF, &                         ! PDAF diagnostic routines
       only: PDAF_diag_stddev, PDAF_diag_variance, PDAFomi_diag_stats
  use parallel_pdaf_mod, &            ! Parallelization variables
       only: COMM_filter, mype_filter
  use statevector_pdaf_mod, &         ! Statevector variables
       only: sfields, n_fields
  use io_pdaf_mod, &                  ! Output file operations
       only: write_pdaf

  implicit none

! *** Arguments ***
  integer, intent(in) :: step         !< Current time step (negative for call after forecast)
  integer, intent(in) :: dim_p        !< Process-local state dimension
  integer, intent(in) :: dim_ens      !< Size of state ensemble
  integer, intent(in) :: dim_ens_p    !< Process-local size of ensemble
  integer, intent(in) :: dim_obs_p    !< Process-local dimension of observation vector
  real, intent(inout) :: state_p(dim_p) !< Process-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  real, intent(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  real, intent(inout) :: ens_p(dim_p, dim_ens)      !< Process-local state ensemble
  integer, intent(in) :: flag         !< PDAF status flag


! *** local variables ***
  integer :: i, j                     ! Counters
  integer :: nobs                     ! Number of observations
  integer :: verbose                  ! Flag for screen output
  integer :: istart, iend             ! stard and end index of a field in state vector
  integer :: pdaf_status              ! status flag
  real :: stddev_g                    ! Global ensemble standard deviation over all fields
  real, allocatable :: ens_stddev(:)  ! ensemble standard deviation for each field (=estimated RMS errors)
  real, allocatable :: variance_p(:)  ! Ensemble variance state vector
  real, pointer :: obsstats_ptr(:,:)  ! Point for observation statistics
  character(len=3) :: anastr          ! String for call type (initial, forecast, analysis)


! **********************
! *** INITIALIZATION ***
! **********************

  if (mype_filter == 0) then
     if (step==0) then
        write (*, '(a, 5x, a)') 'model-PDAF', 'Analyze initial state ensemble'
        anastr = 'ini'
     elseif (step<0) then
        write (*, '(a, 5x, a)') 'model-PDAF', 'Analyze and write forecasted state ensemble'
        anastr = 'for'
     else
        write (*, '(a, 5x, a)') 'model-PDAF', 'Analyze and write assimilated state ensemble'
        anastr = 'ana'
     end if
  end if


! ************************************************************
! *** Compute ensemble mean and standard deviation         ***
! *** (=RMS errors according to sampled covar matrix)      ***
! ************************************************************

  ! Allocate fields
  allocate(ens_stddev(n_fields))

  ! Compute ensemble deviation and mean separately
  ! for each field in the state vector
  do j = 1, n_fields
     ! Start and end index
     istart = 1 + sfields(j)%off
     iend = sfields(j)%dim + sfields(j)%off

     call PDAF_diag_stddev(sfields(j)%dim, dim_ens, &
          state_p(istart:iend), ens_p(istart:iend,:), &
          ens_stddev(j), 1, COMM_filter, pdaf_status)
  end do

  ! Compute ensemble variance in state vector format

  allocate(variance_p(dim_p))

  call PDAF_diag_variance(dim_p, dim_ens, state_p, ens_p, variance_p, &
     stddev_g, 0, 0, COMM_filter, pdaf_status)

  ! Compute observation diagnostics

  verbose = 0
  if (mype_filter==0) verbose = 1
  call PDAFomi_diag_stats(nobs, obsstats_ptr, verbose)


! *****************
! *** Screen IO ***
! *****************

  ! Output ensemble standard deviations
  if (mype_filter == 0) then
     write (*, '(a,6x,a)') 'model-PDAF', 'Ensemble standard deviation (estimated RMS error)'
     do i = 1, n_fields
        write (*,'(a,4x,a13,4x,a10,2x,es12.4)') &
             'model-PDAF', 'stddev-'//anastr, trim(sfields(i)%name), ens_stddev(i)
     end do
  end if


! *******************
! *** File output ***
! *******************

  call write_pdaf(step, dim_p, dim_ens, state_p, ens_p, variance_p)


! *******************
! *** Clean up    ***
! *******************

  deallocate(variance_p)

end subroutine prepoststep_pdaf_offline
