!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
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
!! **Using this template:**
!!   To be able to distinguish the observation type and the routines in this module,
!!   we recommend to rename the module according to the observation module-type.
!!   Further,we recommend to replace 'OBSTYPE' in the routine names according to the
!!   type of the observation so that they can be identified when calling them from 
!!   the call-back routines.
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
!! In addition, there is one optional routines, which is required if filters 
!! with localization are used:
!! * init_dim_obs_l_OBSTYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
module obs_OBSTYPE_pdafomi

  use PDAF, &
       only: obs_f, obs_l   ! Declaration of observation data types
 
  implicit none
  save

!+++ Specific part for the observation type

  ! Variables which are inputs to the module (usually set in init_pdaf)
  logical :: assim_OBSTYPE        !< Whether to assimilate this data type
  real    :: rms_obs_OBSTYPE      !< Observation error standard deviation (for constant errors)

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
!! The routine is called by each filter process at the beginning of
!! the analysis step before the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the observation type
!! handled in this module according to the current time step for all
!! observations required for the analyses in the loop over all local 
!! analysis domains on the process-local state domain.
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
!! Further variables are set by PDAF-OMI in the routine PDAFomi_gather_obs.
!!
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  subroutine init_dim_obs_OBSTYPE(step, dim_obs)

    use PDAF, &                          ! PDAF 
         only: PDAF_local_type, PDAFomi_gather_obs, PDAFomi_set_localize_covar
    use parallel_pdaf_mod, &             ! Parallelization variables
         only: myproc_assim
    use assim_pdaf_mod, &                ! Variables for assimilation
         only: filtertype, cradius, coords_p, dim_state_p, &
         locweight, cradius, sradius, twin_experiment
    use statevector_pdaf_mod, &          ! Statevector variables
         only: id, sfields
!    use synobs_pdaf_mod, &               ! Routines for synthetic observations
!         only: file_synobs, init_file_syn_obs, read_syn_obs

    ! Specific for the model
!    use model_pdaf_mod, &                ! Model variables
!         only: nx, ny, nx_p, n_dim

    implicit none

! *** Arguments ***
    integer, intent(in)    :: step       !< Current time step
    integer, intent(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    integer :: i, j                      ! Counters
    integer :: dim_obs_p                 ! Number of process-local observations
    real, allocatable :: obs_p(:)        ! Process-local observation vector
    real, allocatable :: ivar_obs_p(:)   ! Process-local inverse observation error variance
    real, allocatable :: ocoord_p(:,:)   ! Process-local observation coordinates 
    logical, save :: firsttime=.true.    ! Flag for first call


    ! Template reminder - delete when implementing functionality
    write (*,*) 'TEMPLATE init_OBSTYPE_pdafomi_TEMPLATE.F90: Initialize observations'

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    if (myproc_assim==0) &
         write (*,'(8x,a)') 'Assimilate observations - OBS_OBSTYPE'

    ! Store whether to assimilate this observation type (used in routines below)
    if (assim_OBSTYPE) thisobs%doassim = 1

! +++ Adapt according to the chosen coordinate system (Cartesian or geographic)

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

! +++ Adapt according to how many coordinate direction are used to compute distances

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! Specify the overall domain size
    ! Onle required for thisobs%use_global_obs = 0
    ! or thisobs%disttype=1 (periodicity)
     allocate(thisobs%domainsize(2))
!     thisobs%domainsize(1) = real(nx)
!     thisobs%domainsize(2) = real(ny)
    
! +++ using thisobs%use_global_obs=0 can improve the computing performance
! +++ for domain-localized filters. It requires that PDAFomi_set_domain_limits
! +++ is called once before (see init_pdaf_offline.F90).

    ! Specify whether to (1) use global observations for local filters,
    ! or (0) restrict the full observations to those relevant for a process domain
    thisobs%use_global_obs = 1


! +++++ This is a dummy implementation for a single observation
! +++++ Its only purpose is to let the program run with segmentation fault

    if (myproc_assim==0) &
         write (*,'(8x,a)') 'Dummy implementation of a single observation'

    dim_obs_p = 1   ! This needs to be determined by counting observations

    ! Allocate process-local observation arrays
    allocate(obs_p(dim_obs_p))
    allocate(ivar_obs_p(dim_obs_p))
    allocate(ocoord_p(thisobs%ncoord, dim_obs_p))

    ! Allocate process-local index array
    ! This array has a many rows as required for the observation operator
    ! 1 if observations are at grid points; >1 if interpolation is required
    allocate(thisobs%id_obs_p(1 , dim_obs_p))

    ! Dummy values
    obs_p = 1.0
    ocoord_p = 1.0
    thisobs%id_obs_p = 1
    ivar_obs_p = 1.0
! +++++ END OF DUMMY OBSERVATION


! ***************************************
! *** Read Process-local observations ***
! ***************************************

  ! read observation values and their coordinates
  ! also read observation error information if available


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

!    dim_obs_p = ...
    

    haveobs: if (dim_obs_p > 0) then
    ! *** Initialize vector of observations on the process sub-domain ***

!      ALLOCATE(obs_p(dim_obs_p))

!      obs_p = ....


    ! *** Initialize coordinate array of observations on the process sub-domain ***

  ! The coordinates are only used in case of the local filters or to compute
  ! interpolation coefficients (see below). In any case, the array must be
  ! allocated because it's an argument of PDAFomi_gather_obs called below.

!      ALLOCATE(ocoord_p(thisobs%ncoord, dim_obs_p))

!      ocoord_p = ....

    ! *** Initialize process local index array                         ***
    ! *** This array holds the information which elements of the state ***
    ! *** vector are used in the observation operator.                 ***
    ! *** It has a many rows as required for the observation operator, ***
    ! *** i.e. 1 if observations are at grid points; >1 if             ***
    ! *** interpolation is required                                    ***

  ! The initialization is done locally for each process sub-domain and later
  ! used in the observation operator. 
  ! Examples:
  ! 1. If the observations are model fields located at grid points, one should
  !   initialize the index array thisobs%id_obs_p with one row so that it contains 
  !   the indices of the observed field values in the process-local state vector
  !   (state_p). Then one can use the observation operator OBS_OP_GRIDPOINT 
  !   provided by the module PDAFomi.
  ! 2. If the observations are the average of model fields located at grid points,
  !   one should initialize the index array thisobs%id_obs_p with as many rows as 
  !   values to be averaged. Each column of the arrays then contains the indices of
  !   the elements of the process-local state vector that have to be averaged. With
  !   this index array one can then use the observation operator OBS_OP_GRIDAVG
  !   provided by the module PDAFomi.
  ! 3. If model values need to be interpolated to the observation location
  !   one should initialize the index array thisobs%id_obs_p with as many rows as 
  !   values are required in the interpolationto be averaged. Each column of the 
  !   array then contains the indices of elements of the process-local state vector 
  !   that are used in the interpolation.
  ! Below, you need to replace NROWS by the number of required rows
  ! Note: This array is only used in the observation operator routines. If you
  !   replace the observation operator routines provided by PDAF-OMI by a custom
  !   observation operator, you might not need this array.

!      ALLOCATE(thisobs%id_obs_p( NROWS , dim_obs_p))

!      thisobs%id_obs_p = ...


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

  ! This initialization is only required if an observation operator
  ! with interpolation is used. The coefficients should be determined
  ! here instead of the observation operator, because the operator is
  ! called for each ensemble member while init_dim_obs is only called
  ! once.

  ! Allocate array of interpolation coefficients. As thisobs%id_obs_p, the number
  ! of rows corresponds to the number of grid points using the the interpolation

!      ALLOCATE(thisobs%icoeff_p( NROWS , dim_obs_p))

  ! Ensure that the order of the coefficients is consistent with the
  ! indexing in thisobs%id_obs_p. Further ensure that the order is consistent
  ! with the assumptions used in the observation operator.

!      thisobs%icoeff_p = ...


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

!      ALLOCATE(ivar_obs_p(dim_obs_p))
    
!      ivar_obs_p = ...


    else haveobs

       ! *** For dim_obs_p=0 allocate arrays with minimum size

       ! NOTE FOR DIM_OBS_P=0
       ! For the call to PDAFomi_gather_obs_f
       ! The arrays obs_p, ivar_obs_p, ocoord_p, and thisobs%id_obs_p
       ! have to to be allocated for the call to PDAFomi_gather_obs.
       ! Thus, if dim_obs_p=0 can happen in your application you should
       ! explicitly handle this case as done with the if-statement
       ! checking dim_obs_p

       allocate(obs_p(1))
       allocate(ocoord_p(thisobs%ncoord, 1))
       allocate(thisobs%id_obs_p(1, 1))
       allocate(ivar_obs_p(1))

    end if haveobs


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *** and overwrite real observation values.            ***
! *********************************************************

! +++ Outcommented because the routine needs netcdf
!     IF (twin_experiment .AND. filtertype/=PDAF_DA_GENOBS) THEN
!        CALL read_syn_obs(file_synobs, dim_obs_p, obs_p, step, 1-myproc_assim)
!     END IF


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    ! This routine is generic for the case that only the observations, 
    ! inverse variances and observation coordinates are gathered

    call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, cradius, dim_obs)


! ***************************************************************
! *** Provide localization information for LEnKF, EnsRF, EAKF ***
! ***************************************************************

    ! If one uses the LEnKF or ENSRF methods, one has to set the
    ! localization information. The three localization variables
    ! (locweight, cradius, sradius) can be different for each
    ! observation type.
    ! One also need to initialize of coordinate array (coords_p)
    ! for the state vector (usually in init_pdaf).

    ! For cradius and sradius:
    ! If these are defined as scalar values, isotropic localization is used.
    ! If these are vectors, nonisotropic localization is used
    !   (their length has to be equal to thisobs%ncoord)


!     IF (PDAF_local_type() > 1) THEN
!        CALL PDAFomi_set_localize_covar(thisobs, dim_state_p, ndim, coords_p, &
!             locweight, cradius, sradius)
!     END IF


! ***************************************************************
! *** When generating synthetic observations: Initialize file ***
! ***************************************************************

! +++ Outcommented because the routine needs netcdf
!     if (filtertype==PDAF_DA_GENOBS .and. firsttime) then
! 
!        ! Initialize file
!        call init_file_syn_obs(maxsize, file_syntobs_OBSTYPE, 1)
! 
!        firsttime = .false.
!     end if


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
!    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Note: Arrays in THISOBS are deallocated by PDAF-OMI

  end subroutine init_dim_obs_OBSTYPE



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
  subroutine obs_op_OBSTYPE(dim_p, dim_obs, state_p, ostate)

    use PDAF, &                ! Include PDAF-OMI routine
         only: PDAFomi_obs_op_gridpoint

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p                 !< Process-local state dimension
    integer, intent(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    real, intent(in)    :: state_p(dim_p)        !< Process-local model state
    real, intent(inout) :: ostate(dim_obs)       !< Full observed state


    ! Template reminder - delete when implementing functionality
    write (*,*) 'TEMPLATE init_OBSTYPE_pdafomi_TEMPLATE.F90: Apply observation operator'

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    !+++  Choose suitable observation operator from the
    !+++  module PDAFomi_obs_op or implement your own

    ! Example: Observation operator for observed grid point values
    call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)

  end subroutine obs_op_OBSTYPE



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
  subroutine init_dim_obs_l_OBSTYPE(domain_p, step, dim_obs, dim_obs_l)

    use PDAF, &                ! Include PDAF-OMI routine
         only: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    ! one can also set observation-specific values for the localization.
    use assim_pdaf_mod, &   
         only: coords_l, cradius, locweight, sradius

    implicit none

! *** Arguments ***
    integer, intent(in)  :: domain_p     !< Index of current local analysis domain
    integer, intent(in)  :: step         !< Current time step
    integer, intent(in)  :: dim_obs      !< Full dimension of observation vector
    integer, intent(inout) :: dim_obs_l  !< Local dimension of observation vector


    ! Template reminder - delete when implementing functionality
    write (*,*) 'TEMPLATE init_OBSTYPE_pdafomi_TEMPLATE.F90: Initialize local observations'

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Here one has to specify the coordinates of the local analysis domain
    ! (coords_l) and the localization variables, which can be different for
    ! each observation type and can be made dependent on the index DOMAIN_P.
    ! coords_l should be set in the call-back routine init_dim_l.

    ! For cradius and sradius:
    ! If these are defined as scalar values, isotropic localization is used.
    ! If these are vectors, nonisotropic localization is used
    !   (their length has to be equal to thisobs%ncoord)

    call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, cradius, sradius, dim_obs_l)

  end subroutine init_dim_obs_l_OBSTYPE

end module obs_OBSTYPE_pdafomi
