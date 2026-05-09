!> Interface routine to call PDAF analysis step
!!
!! This routine performs a single analysis step in the
!! offline implementation of PDAF. For this, it calls the
!! filter-specific assimilation routine of PDAF. For the
!! offline implementation this is PDAF3_assim_offline.
!!
!! In this routine, the real names of most of the 
!! user-supplied routines for PDAF are specified (see below).
!!
!! The routine is generic, but has to be part of the user code
!! because one might want to adapt the names of call-back routines.
!!
!! __Revision history:__
!! * 2009-11 - Lars Nerger - Initial code by restructuring
!! * 2026-02 - Lars Nerger - Revision for advanced tutorial revising tutorial case
!! * Other revisions - see repository log
!!
subroutine assimilate_pdaf_offline()

  use PDAF, &                     ! PDAF interface definitions
       only: PDAF3_assim_offline, PDAF3_generate_obs_offline, PDAF_abort, PDAF_DA_GENOBS
  use parallel_pdaf_mod, &        ! Parallelization variables
       only: myproc_ens
  use assim_pdaf_mod, &           ! Variables for assimilation
       only: filtertype

  implicit none

! *** Local variables ***
  integer :: status_pdaf          ! PDAF status flag


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the call to PDAF3_assim_offline.
! This allows the user to specify the actual name of a routine.  
! The PDAF-internal name of a subroutine might be different from the external name!

  ! Interface between model and PDAF, and prepoststep
  external :: prepost_cb_pdaf           ! User supplied pre/poststep routine
  ! Localization of state vector
  external :: init_n_domains_cb_pdaf, & ! Provide number of local analysis domains
       init_dim_l_cb_pdaf               ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  external :: init_dim_obs_pdafomi, &   ! Get dimension of full obs. vector for Process-local domain
       obs_op_pdafomi, &                ! Obs. operator for full obs. vector for Process-local domain
       init_dim_obs_l_pdafomi           ! Get dimension of obs. vector for local analysis domain
  ! Subroutine used for generating observations
  external :: get_obs_cb_pdaf           ! Get vector of synthetic observations from PDAF


! *********************************
! *** Call assimilation routine ***
! *********************************

! +++ Note: The universal routine PDAF3_assim_offline can be used to
! +++ execute all filter methods. The specified routines for localization
! +++ are only executed if a local filter is used. If one uses
! +++ exclusively global filters, the LEnKF, EnsRF or EAKF, one can use
! +++ the specific routine PDAF3_assim_offline_global which does not
! +++ include the arguments for localization. This would avoid to include
! +++ routines that are never called for global filters. 

  ! Call universal PDAF3 ensemble assimilation routine
  if (filtertype /= PDAF_DA_GENOBS) then
     call PDAF3_assim_offline( &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          init_n_domains_cb_pdaf, init_dim_l_cb_pdaf, init_dim_obs_l_pdafomi, &
          prepost_cb_pdaf, status_pdaf)
  else
     ! Observation generation has its own OMI interface routine
     call PDAF3_generate_obs_offline( &
          init_dim_obs_pdafomi, obs_op_pdafomi, get_obs_cb_pdaf, &
          prepost_cb_pdaf, status_pdaf)
  end if


! ************************
! *** Check error flag ***
! ************************

  if (status_pdaf /= 0) then
     write (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF3_assim_offline - stopping! (Process ', myproc_ens,')'
     call PDAF_abort(1)
  end if

end subroutine assimilate_pdaf_offline
