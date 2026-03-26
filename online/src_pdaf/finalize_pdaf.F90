!>  Finalize PDAF
!!
!! This routine calls routines for output on timing
!! and memory use, to deallocate PDAF-internal arrays,
!! and to finalize MPI of model is not parallelized.
!!
!! The routine is generic, but has to be part of the user code
!! because one might to use different outputs.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine finalize_pdaf()

  use PDAF, &                     ! PDAF
       only: PDAF_print_info, PDAF_deallocate
  use parallel_pdaf_mod, &        ! Parallelization
       only: mype_ens

  implicit none
  
! *** Show globally allocated memory for PDAF ***
  call PDAF_print_info(11)

! *** Print PDAF timings onto screen ***
  if (mype_ens==0) call PDAF_print_info(3)

! *** Deallocate PDAF arrays ***
  call PDAF_deallocate()

end subroutine finalize_pdaf
