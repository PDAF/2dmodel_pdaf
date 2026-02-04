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
subroutine prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  use mpi                           ! MPI
  use model_pdaf_mod, &             ! Model variables
       only: nx, ny, nx_p
  use parallel_pdaf_mod, &     ! Parallelization variables
       only: COMM_filter, mype_filter, npes_filter, MPIerr, MPIstatus
  use assimilation_pdaf_mod, &      ! Assimilation variables
       only: dim_state
  use PDAF, &                  ! PDAF diagnostic routine
       only: PDAF_diag_stddev

  implicit none

! *** Arguments ***
  integer, intent(in) :: step        !< Current time step (negative for call after forecast)
  integer, intent(in) :: dim_p       !< Process-local state dimension
  integer, intent(in) :: dim_ens     !< Size of state ensemble
  integer, intent(in) :: dim_ens_p   !< Process-local size of ensemble
  integer, intent(in) :: dim_obs_p   !< Process-local dimension of observation vector
  real, intent(inout) :: state_p(dim_p) !< Process-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  real, intent(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  real, intent(inout) :: ens_p(dim_p, dim_ens)      !< Process-local state ensemble
  integer, intent(in) :: flag        !< PDAF status flag


! *** local variables ***
  integer :: i, j, member             ! Counters
  integer :: pdaf_status              ! status flag
  logical, save :: firsttime = .true. ! Routine is called for first time?
  real :: ens_stddev                  ! estimated RMS error
  real, allocatable :: field(:,:)     ! global model field
  character(len=2) :: ensstr          ! String for ensemble member
  character(len=2) :: stepstr         ! String for time step
  character(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
  ! Variables for parallelization - global fields
  integer :: domain                   ! Counter
  integer :: off_p                    ! Row-offset according to domain decomposition
  real, allocatable :: ens(:,:)       ! global ensemble
  real, allocatable :: state(:)       ! global state vector
  real,allocatable :: ens_p_tmp(:,:)  ! Temporary ensemble for some Process-domain


! **********************
! *** INITIALIZATION ***
! **********************

  if (mype_filter == 0) then
     if (firsttime) then
        write (*, '(a, 5x, a)') 'model-PDAF', 'Analyze initial state ensemble'
        anastr = 'ini'
     else
        if (step<0) then
           write (*, '(a, 5x, a)') 'model-PDAF', 'Analyze and write forecasted state ensemble'
           anastr = 'for'
        else
           write (*, '(a, 5x, a)') 'model-PDAF', 'Analyze and write assimilated state ensemble'
           anastr = 'ana'
        end if
     end if
  end if


! ************************************************************
! *** Compute ensemble mean and standard deviation         ***
! *** (=RMS errors according to sampled covar matrix)      ***
! ************************************************************

  call PDAF_diag_stddev(dim_p, dim_ens, state_p, ens_p, &
        ens_stddev, 1, COMM_filter, pdaf_status)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  if (mype_filter == 0) then
     write (*, '(a, 4x, a, es12.4)') &
          'model-PDAF', 'RMS error according to sampled variance: ', ens_stddev
  end if


! *******************
! *** File output ***
! *******************

  notfirst: if (.not. firsttime) then

     allocate(ens(dim_state, dim_ens))
     allocate(state(dim_state))

     ! Gather full ensemble on process with rank 0 and write file
     mype0b: if (mype_filter /= 0) then

        ! *** Send ensemble substates on filter-Processs with rank > 0 ***

        call MPI_send(ens_p, dim_ens * dim_p, &
             MPI_DOUBLE_PRECISION, 0, 1, COMM_filter, MPIerr)

     else mype0b

        ! *** Initialize and receive sub-states on Process 0 ***

        ! Initialize sub-ensemble for Process 0
        do member = 1, dim_ens
           do i=1, dim_p
              ens(i, member) = ens_p(i, member)
           end do
        end do

        ! Define offset in state vectors
        off_p = dim_p

        do domain = 2, npes_filter
           ! Initialize sub-ensemble for other Processs and send sub-arrays

           ! Allocate temporary buffer array
           allocate(ens_p_tmp(nx_p*ny, dim_ens))

           ! Receive sub-arrays
           call MPI_recv(ens_p_tmp, nx_p*ny * dim_ens, MPI_DOUBLE_PRECISION, &
                domain - 1, 1, COMM_filter, MPIstatus, MPIerr)

           ! Initialize MPI buffer for local ensemble
           do member = 1, dim_ens
              do i = 1, nx_p*ny
                 ens(i + off_p, member) = ens_p_tmp(i, member)
              end do
           end do

           deallocate(ens_p_tmp)

           ! Increment offset
           off_p = off_p + nx_p*ny

        end do


        ! *** Now write analysis ensemble ***

        write (*, '(a, 5x, a)') 'model-PDAF', '--- write ensemble and state estimate'

        ! Set string for time step
        if (step>=0) then
           write (stepstr, '(i2.2)') step
        else
           write (stepstr, '(i2.2)') -step
        end if

        allocate(field(ny, nx))

        do member = 1, dim_ens
           do j = 1, nx
              field(1:ny, j) = ens(1 + (j-1)*ny : j*ny, member)
           end do

           write (ensstr, '(i2.2)') member

           open(11, file = 'ens_'//trim(ensstr)//'_step'//trim(stepstr)//'_'//trim(anastr)//'.txt', status = 'replace')
 
           do i = 1, ny
              write (11, *) field(i, :)
           end do

           close(11)
        end do

     end if mype0b

     ! Gather full state vector on process with rank 0 and write to file
     mype0c: if (mype_filter /= 0) then

        ! *** Send ensemble substates on filter-Processs with rank > 0 ***

        call MPI_send(state_p, dim_p, &
             MPI_DOUBLE_PRECISION, 0, 1, COMM_filter, MPIerr)

     else mype0c

        ! *** Initialize and receive sub-states on Process 0 ***

        ! Initialize sub-state for Process 0
        do i = 1, dim_p
           state(i) = state_p(i)
        end do

        ! Define offset in state vectors
        off_p = dim_p

        do domain = 2, npes_filter
           ! Initialize sub-ensemble for other Processs and send sub-arrays

           ! Receive sub-arrays
           call MPI_recv(state(1+off_p), nx_p*ny, MPI_DOUBLE_PRECISION, &
                domain - 1, 1, COMM_filter, MPIstatus, MPIerr)

           ! Increment offset
           off_p = off_p + nx_p*ny

        end do
     
        ! *** Now write analysis state estimate ***

        do j = 1, nx
           field(1:ny, j) = state(1 + (j-1)*ny : j*ny)
        end do

        open(11, file = 'state_step'//trim(stepstr)//'_'//trim(anastr)//'.txt', status = 'replace')
 
        do i = 1, ny
           write (11, *) field(i, :)
        end do

        close(11)

        deallocate(field)
     end if mype0c

     deallocate(ens, state)

  end if notfirst


! ********************
! *** finishing up ***
! ********************

  firsttime = .false.

end subroutine prepoststep_pdaf
