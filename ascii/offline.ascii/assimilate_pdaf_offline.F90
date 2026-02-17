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
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine assimilate_pdaf_offline()

  use PDAF, &                     ! PDAF interface definitions
       only: PDAF3_assim_offline
  use parallel_pdaf_mod, &        ! Parallelization variables
       only: mype_ens, abort_parallel

  implicit none

! *** Local variables ***
  integer :: status_pdaf          ! PDAF status flag


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the calls to PDAF3_assimilate.
! This allows the user to specify the actual name of a routine.  
! The PDAF-internal name of a subroutine can be different from the external name!

  ! Interface between model and PDAF, and prepoststep
  external :: prepoststep_pdaf        ! User supplied pre/poststep routine
  ! Localization of state vector
  external :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf                ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  external :: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for Process-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for Process-local domain
       init_dim_obs_l_pdafomi         ! Get dimension of obs. vector for local analysis domain


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
  call PDAF3_assim_offline( &
       init_dim_obs_pdafomi, obs_op_pdafomi, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
       prepoststep_pdaf, status_pdaf)


! ************************
! *** Check error flag ***
! ************************

  if (status_pdaf /= 0) then
     write (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF3_assimilate - stopping! (Process ', mype_ens,')'
     call abort_parallel()
  end if

end subroutine assimilate_pdaf_offline
