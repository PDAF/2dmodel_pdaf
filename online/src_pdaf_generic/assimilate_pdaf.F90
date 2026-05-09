!>  Routine to call PDAF for analysis step in fully-parallel mode
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the PDAF assimilation routine PDAF3_assimilate,
!! which checks whether the forecast phase is completed. If so, the
!! analysis step is computed inside PDAF.
!!
!! In this routine, the real names of most of the 
!! user-supplied routines for PDAF are specified (see below).
!!
!! The routine is generic, but has to be part of the user code
!! because one might want to adapt the names of call-back routines.
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code
!! * 2026-02 - Lars Nerger - Revision for advanced tutorial revising tutorial case
!! * Other revisions - see repository log
!!
subroutine assimilate_pdaf()

  use PDAF, &                     ! PDAF interface definitions
       only: PDAF3_assimilate, PDAF3_generate_obs, PDAF_abort, PDAF_DA_GENOBS
  use parallel_pdaf_mod, &        ! Parallelization variables
       only: myproc_ens
  use assim_pdaf_mod, &           ! Variables for assimilation
       only: filtertype

  implicit none

! *** Local variables ***
  integer :: status_pdaf          ! PDAF status flag


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the calls to PDAF3_assimilate.
! This allows the user to specify the actual name of a routine.  
! The PDAF-internal name of a subroutine can be different from the external name!

  ! Interface between model and PDAF, and prepoststep
  external :: coll_state_cb_pdaf, &   ! Collect a state vector from model fields
       dist_state_cb_pdaf, &          ! Distribute a state vector to model fields
       next_obs_cb_pdaf, &            ! Provide time step of next observation
       prepost_cb_pdaf                ! User supplied pre/poststep routine
  ! Localization of state vector
  external :: init_dim_l_cb_pdaf, &   ! Initialize state dimension for local analysis domain
       init_n_domains_cb_pdaf         ! Provide number of local analysis domains
       
  ! Interface to PDAF-OMI for local and global filters
  external :: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for Process-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for Process-local domain
       init_dim_obs_l_pdafomi         ! Get dimension of obs. vector for local analysis domain
  ! Subroutine used for generating observations
  external :: get_obs_cb_pdaf         ! Get vector of synthetic observations from PDAF


! *********************************
! *** Call assimilation routine ***
! *********************************

! +++ Note: The universal routine PDAF3_assimilate can be used to
! +++ execute all filter methods. The specified routines for localization
! +++ are only executed if a local filter is used. If one uses
! +++ exclusively global filters or the LEnKF, one can use the specific
! +++ routine PDAF3_assimilate_global which does not include the
! +++ arguments for localization. This would avoid to include routines
! +++ that are never called for global filters. 

  ! Call universal PDAF3 ensemble assimilation routine
  if (filtertype /= PDAF_DA_GENOBS) then
     call PDAF3_assimilate(coll_state_cb_pdaf, dist_state_cb_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          init_n_domains_cb_pdaf, init_dim_l_cb_pdaf, init_dim_obs_l_pdafomi, &
          prepost_cb_pdaf, next_obs_cb_pdaf, status_pdaf)
  else
     ! Observation generation has its own OMI interface routine
     call PDAF3_generate_obs(coll_state_cb_pdaf, dist_state_cb_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, get_obs_cb_pdaf, &
          prepost_cb_pdaf, next_obs_cb_pdaf, status_pdaf)
  end if


! ************************
! *** Check error flag ***
! ************************

  if (status_pdaf /= 0) then
     write (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF3_assimilate - stopping! (Process ', myproc_ens,')'
     call PDAF_abort(1)
  end if

end subroutine assimilate_pdaf
