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
!! * Later revisions - see repository log
!!
subroutine assimilate_pdaf_offline()

  use PDAF, &                     ! PDAF interface definitions
       only: PDAF3_assim_offline, PDAF_abort
  use parallel_pdaf_mod, &        ! Parallelization variables
       only: myproc_ens

  implicit none

! *** Local variables ***
  integer :: status_pdaf          ! PDAF status flag


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the call to PDAF3_assim_offline.
! This allows the user to specify the actual name of a routine.  
! The PDAF-internal name of a subroutine might be different from the external name!

  ! Interface between model and PDAF, and prepoststep
  external :: prepoststep_pdaf_offline ! User supplied pre/poststep routine
  ! Localization of state vector
  external :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf                 ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  external :: init_dim_obs_pdafomi, &  ! Get dimension of full obs. vector for Process-local domain
       obs_op_pdafomi, &               ! Obs. operator for full obs. vector for Process-local domain
       init_dim_obs_l_pdafomi          ! Get dimension of obs. vector for local analysis domain


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
  call PDAF3_assim_offline( &
       init_dim_obs_pdafomi, obs_op_pdafomi, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
       prepoststep_pdaf_offline, status_pdaf)


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
