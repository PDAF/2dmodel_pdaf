!> Module to include variables from 2D model
!!
!! This module includes variables from the 2D tutorial model.
!! For PDAF in its online coupling it is the single point which
!! directly links to the 2D-model. All other routines include
!! from this module. For the offline case the variables are
!! declared here and separately initialized.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_pdaf_mod

  use model_mod, &                ! Model variables
       ONLY: nx, ny, nx_p, n_dim, coords_x_p, coords_y_p, &
       total_steps, fieldA_p, fieldB_p

  implicit none


! *** Additional variables ***

! One might like to declare here additional variables relating
! to the model fields or model grid only used for DA.

end module model_pdaf_mod

 
