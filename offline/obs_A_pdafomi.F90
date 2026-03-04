!> PDAF-OMI observation module for type A observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! OBSTYPE = A
!!
!! __Observation type A:__
!! The observation type A in this tutorial are observations of fieldA
!! at specified model grid points.
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF
!! usually by callback_obs_pdafomi.F90
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data types (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_OBSTYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_OBSTYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there is one optional routine, which is required if filters 
!! with localization are used:
!! * init_dim_obs_l_OBSTYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module obs_A_pdafomi

  use PDAF, &
       only: obs_f, obs_l   ! Declaration of observation data types
 
  implicit none
  save

!+++ Specific part for the observation type

  ! Variables which are inputs to the module (usually set in init_pdaf)
  logical :: assim_A = .false.     !< Whether to assimilate this data type
  real    :: rms_obs_A = 0.25      !< Observation error standard deviation (for constant errors)
  character(len=100) :: file_obs_A = 'obsA.nc'     !< Name of observation file
  character(len=100) :: stub_synobs_A='obsA_syn_'  !< Name stub of file for synthetic observations

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.

!+++ End of specific part


! *********************************************************
! *** Data type obs_f defines the full observations     ***
! *** Data type obs_l defines the local observations    ***
! ***  for domain-localized ensemble methods            ***
! *********************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could rename the variables
  type(obs_f), target, public :: thisobs      ! full observation
  type(obs_l), target, public :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

contains

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the Process-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - array with indices of module-type observation in process-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!! * thisobs\%inno_omit   - Omit obs. if squared innovation larger this factor times observation variance
!!                          (default: 0.0, omission is active if >0) 
!! * thisobs\%inno_omit_ivar - Value of inverse variance to omit observation
!!                          (default: 1.0e-12, change this if this value is not small compared to actual obs. error)
!!
!! Further variables are set by PDAF-OMI when the routine PDAFomi_gather_obs is called.
!!
  subroutine init_dim_obs_A(step, dim_obs)

    use netcdf
    use PDAF, &                          ! PDAF
         only: PDAFomi_gather_obs, PDAF_local_type, PDAFomi_set_localize_covar, &
         PDAF_DA_GENOBS
    use parallel_pdaf_mod, &             ! Parallelization variables
         only: mype_filter
    use assimilation_pdaf_mod, &         ! Variables for assimilation
         only: filtertype, cradius, coords_p, dim_state_p, &
         locweight, cradius, sradius, twin_experiment
    use statevector_pdaf_mod, &          ! Statevector variables
         only: id, sfields
    use synobs_pdaf_mod, &               ! Routines for synthetic observations
         only: file_synobs, init_file_syn_obs, read_syn_obs
    use io_pdaf_mod, &                   ! IO operations
         only: nfcheck

    ! Specific for 2D tutorial model
    use model_pdaf_mod, &                ! Model variables
         only: nx, ny, nx_p, n_dim

    implicit none

! *** Arguments ***
    integer, intent(in)    :: step       !< Current time step
    integer, intent(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    integer :: i, j                      ! Counters
    integer :: cnt_p, cnt0_p             ! Counters
    integer :: off_nx                    ! Offset of local grid in global domain in x-direction
    integer :: dim_obs_p                 ! Number of process-local observations
    real, allocatable :: obs_field(:,:)  ! Observation field read from file
    real, allocatable :: obs_p(:)        ! Process-local observation vector
    real, allocatable :: ivar_obs_p(:)   ! Process-local inverse observation error variance
    real, allocatable :: ocoord_p(:,:)   ! Process-local observation coordinates 
    integer :: ncid, id_obs              ! variables for netcdf file reading
    integer :: countv(3), startv(3)      ! Vectors for NC operations
    character(len=4) :: procstr          ! 4-digit string for process rank
    logical, save :: firsttime=.true.    ! Flag for first call


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    if (mype_filter==0) &
         write (*,'(a,5x,a)') 'model-PDAF','Assimilate observations - obs type A'

    ! Store whether to assimilate this observation type (used in routines below)
    if (assim_A) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! Specify the overall domain size
    ! ONLY REQUIRED FOR THISOBS%USE_GLOBAL_OBS = 0
    ! OR THISOBS%DISTTYPE = 1 (periodicity)
    allocate(thisobs%domainsize(2))
    thisobs%domainsize(1) = real(nx)
    thisobs%domainsize(2) = real(ny)
    
    ! Specify whether to (1) use global observations for local filters,
    ! or (0) restrict the full observations to those relevant for a process domain
    thisobs%use_global_obs = 1


! **********************************
! *** Read Process-local observations ***
! **********************************

    ! Read observation field from file
    if (mype_filter==0) &
         write (*,'(a,5x,a, i6)') 'model-PDAF','--- read observation at step', step

    allocate(obs_field(ny, nx))

    call nfcheck( NF90_OPEN(trim(file_obs_A), NF90_NOWRITE, ncid))
    call nfcheck( NF90_INQ_VARID(ncid, 'obs', id_obs))
    startv(3) = step
    countv(3) = 1
    startv(2) = 1
    countv(2) = nx
    startv(1) = 1
    countv(1) = ny
    call nfcheck( NF90_GET_VAR(ncid, id_obs, obs_field, start=startv(1:3), count=countv(1:3)))
    call nfcheck( NF90_CLOSE(ncid))


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    ! Get offset of local domain in global domain in x-direction
    off_nx = nx_p*mype_filter

    ! Count process-local observations
    cnt_p = 0
    do j = 1 + off_nx, nx_p + off_nx
       do i= 1, ny
          if (obs_field(i,j) > -999.0) cnt_p = cnt_p + 1
       end do
    end do

    ! Set number of local observations
    dim_obs_p = cnt_p


    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***

    haveobs: if (dim_obs_p > 0) then

       ! Allocate process-local observation arrays
       allocate(obs_p(dim_obs_p))
       allocate(ivar_obs_p(dim_obs_p))
       allocate(ocoord_p(thisobs%ncoord, dim_obs_p))

       ! Allocate process-local index array
       ! This array has a many rows as required for the observation operator
       ! 1 if observations are at grid points; >1 if interpolation is required
       allocate(thisobs%id_obs_p(1, dim_obs_p))

       cnt_p = 0
       cnt0_p = 0
       do j = 1 + off_nx, nx_p + off_nx
          do i= 1, ny
             cnt0_p = cnt0_p + 1
             if (obs_field(i,j) > -999.0) then
                cnt_p = cnt_p + 1
                thisobs%id_obs_p(1, cnt_p) = cnt0_p + sfields(id%fieldA)%off
                obs_p(cnt_p) = obs_field(i, j)
                ocoord_p(1, cnt_p) = real(j)
                ocoord_p(2, cnt_p) = real(i)
             end if
          end do
       end do


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

       ! *** Set inverse observation error variances ***

       ivar_obs_p(:) = 1.0 / (rms_obs_A*rms_obs_A)

    else haveobs

       ! *** For dim_obs_p=0 allocate arrays with minimum size

       allocate(obs_p(1))
       allocate(ivar_obs_p(1))
       allocate(ocoord_p(2, 1))
       allocate(thisobs%id_obs_p(1, 1))

    end if haveobs


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *** and overwrite real observation values.            ***
! *********************************************************

    IF (twin_experiment .AND. filtertype/=PDAF_DA_GENOBS) THEN

       ! Set file name (separate files for each process)
       write (procstr, '(i4.4)') mype_filter
       file_synobs = trim(stub_synobs_A)//procstr//'.nc'

       ! Read synthetic observations
       CALL read_syn_obs(file_synobs, dim_obs_p, obs_p, step, 1-mype_filter)
    END IF


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, cradius, dim_obs)


! ***************************************************************
! *** Provide localization information for LEnKF, EnsRF, EAKF ***
! ***************************************************************

    IF (PDAF_local_type() > 1) THEN
       CALL PDAFomi_set_localize_covar(thisobs, dim_state_p, n_dim, coords_p, &
            locweight, cradius, sradius)
    END IF


! ***************************************************************
! *** When generating synthetic observations: Initialize file ***
! ***************************************************************

    if (filtertype==PDAF_DA_GENOBS .and. firsttime) then

       ! Set file name (separate files for each process)
       write (procstr, '(i4.4)') mype_filter
       file_synobs = trim(stub_synobs_A)//procstr//'.nc'

       ! Initialize file
       call init_file_syn_obs(nx_p*ny, file_synobs, 1)

       firsttime = .false.
    end if


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    deallocate(obs_field)
    deallocate(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  end subroutine init_dim_obs_A



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  subroutine obs_op_A(dim_p, dim_obs, state_p, ostate)

    use PDAF, &
         only: PDAFomi_obs_op_gridpoint

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p                 !< Process-local state dimension
    integer, intent(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    real, intent(in)    :: state_p(dim_p)        !< Process-local model state
    real, intent(inout) :: ostate(dim_obs)       !< Full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    ! observation operator for observed grid point values
    call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)

  end subroutine obs_op_A



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  subroutine init_dim_obs_l_A(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    use PDAF, only: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    use assimilation_pdaf_mod, &   
         only: coords_l, cradius, locweight, sradius

    implicit none

! *** Arguments ***
    integer, intent(in)  :: domain_p     !< Index of current local analysis domain
    integer, intent(in)  :: step         !< Current time step
    integer, intent(in)  :: dim_obs      !< Full dimension of observation vector
    integer, intent(inout) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, cradius, sradius, dim_obs_l)

  end subroutine init_dim_obs_l_A

end module obs_A_pdafomi
