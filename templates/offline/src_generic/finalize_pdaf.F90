!> Interface routine to PDAF for finalization
!!
!! This routine calls routines for outputs on timing and
!! memory use, to deallocate PDAF-internal arrays, and to
!! finalize MPI if MPI was initialized by PDAF, as is
!! usually done in the offline coupled model and in the
!! online coupled model if the model is not parallelized.
!!
!! The routine is generic, but has to be part of the user
!! code because it includes myproc_ens from parallel_pdaf_mod. 
!! Further, one might to use different outputs.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine finalize_pdaf()

  use PDAF, &                     ! PDAF
       only: PDAF_print_info, PDAF_finalize
  use parallel_pdaf_mod, &        ! Parallelization
       only: myproc_ens

  implicit none
  
! *** Show globally allocated memory for PDAF ***
  call PDAF_print_info(11)

! *** Print PDAF timings onto screen ***
  if (myproc_ens==0) call PDAF_print_info(3)

! *** Deallocate PDAF arrays and finalize MPI ***
  call PDAF_finalize()

end subroutine finalize_pdaf
