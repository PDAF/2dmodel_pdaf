!>  Initialize communicators for PDAF
!!
!! Parallelization routine for a model with attached PDAF. The subroutine is
!! called in the main program subsequently to the initialization of MPI. It
!! initializes MPI communicators for the model tasks, filter task and the
!! coupling between model and filter tasks. In addition some other variables 
!! for the parallelization are initialized.
!! The communicators and variables are handed over to PDAF in the call to 
!! PDAF_set_parallel toward the end of this routine.
!!
!! 3 Communicators are generated:
!! * _COMM_filter_: Communicator in which the filter itself operates
!! * _COMM_model_: Communicators for parallel model forecasts
!! * _COMM_couple_: Communicator for coupling between models and filter
!!
!! In addition there is the main communicator
!! * _COMM_ensemble_: The main communicator in which PDAF operates
!! COMM_ensemble is set to the communicator in which all model integration 
!! are computed. Typically, this is MPI_COMM_WORLD, but it can be defined
!! differently if the model only operators on a subset to MPI_COMM_WORLD.
!! This happens, e.g. if some processes are separated to operate an
!! I/O server or a model coupler for coupled model systems.
!!
!! Other variables that have to be initialized are:
!! * _filterpe_ - Logical: Does the Process execute the filter?
!! * _task_id_ - Integer: Index identifying the model task
!! * _my_ensemble_ - Integer: The index of the Process's model task
!! * _local_npes_model_ - Integer array holding numbers of Processs per model task
!!
!! For COMM_filter and COMM_model also the size of the communicators
!! (npes_filter and npes_model) and the rank of each process  (mype_filter,
!! mype_model) are initialized. These variables can be used in the model part 
!! of the program, but are not handed over to PDAF.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
subroutine init_parallel_pdaf_offline(screen)

  use mpi                         ! MPI
  use PDAF, &                     ! PDAF routines
       only: PDAF3_set_parallel
  use parallel_pdaf_mod, &        ! PDAF parallelization variables
       only: mype_ens, npes_ens, COMM_ensemble, &
       mype_model, npes_model, COMM_model, &
       mype_couple, npes_couple, COMM_couple, &
       mype_filter, npes_filter, COMM_filter, filterpe, &
       n_modeltasks, task_id, local_npes_model, MPIerr
  use parser, &                   ! Command line parser
       only: parse

  implicit none

! *** Arguments ***
  integer, intent(in)    :: screen           !< Whether screen information is shown

! *** Local variables ***
  integer :: i, j                     ! Counters
  integer :: pe_index                 ! Index of Process
  integer :: my_color, color_couple   ! Variables for communicator-splitting 
  integer :: flag                     ! Status flag
  integer :: dim_ens                  ! Ensemble size / number of model tasks
  character(len=32) :: handle         ! Handle for command line parser
  logical :: online_coupling          ! Whether to ru nonline coupled DA


  ! Specify online of offline coupling
  online_coupling = .false.

  ! *** Define ensemble communicator.                   ***
  ! *** This is the communicator in which PDAF operates ***

  COMM_ensemble = MPI_COMM_WORLD


  ! *** Parse number of model tasks ***
  ! *** The module variable is N_MODELTASKS. Since it has to be equal
  ! *** to the ensemble size we parse dim_ens from the command line.

  handle = 'dim_ens'
  call parse(handle, dim_ens)
  if (online_coupling) n_modeltasks = dim_ens

  ! *** Fix number or model tasks for offline DA ***
  if (.not. online_coupling) n_modeltasks = 1

  ! *** Get rank and size of COMM_ensemble ***

  call MPI_Comm_Size(COMM_ensemble, npes_ens, MPIerr)
  call MPI_Comm_Rank(COMM_ensemble, mype_ens, MPIerr)

  ! *** Initialize communicators for ensemble evaluations ***
  if (mype_ens == 0) &
       write (*, '(/a, 2x, a)') 'model-PDAF', 'Initialize communicators for assimilation with PDAF'


  ! *** Check consistency of number of parallel ensemble tasks ***
  if (online_coupling) then
     consist1: if (n_modeltasks > npes_ens) then
        ! *** # parallel tasks is set larger than available Processs ***
        n_modeltasks = npes_ens
        if (mype_ens == 0) write (*, '(a, 3x, a)') &
             'model-PDAF', '!!! Resetting number of parallel ensemble tasks to total number of Processs!'
     end if consist1
     if (dim_ens > 0) then
        ! Check consistency with ensemble size
        consist2: if (n_modeltasks > dim_ens) then
           ! # parallel ensemble tasks is set larger than ensemble size
           n_modeltasks = dim_ens
           if (mype_ens == 0) write (*, '(a, 5x, a)') &
                'model-PDAF', '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
        end if consist2
     end if
  end if


  ! *** Store # Processs per ensemble                 ***
  ! *** used for info on Process 0 and for generation ***
  ! *** of model communicators on other Pes      ***
  allocate(local_npes_model(n_modeltasks))

  local_npes_model = floor(real(npes_ens) / real(n_modeltasks))
  do i = 1, (npes_ens - n_modeltasks * local_npes_model(1))
     local_npes_model(i) = local_npes_model(i) + 1
  end do


  ! ***              COMM_MODEL               ***
  ! *** Generate communicators for model runs ***
  ! *** (Split COMM_ENSEMBLE)                 ***
  pe_index = 0
  doens1: do i = 1, n_modeltasks
     do j = 1, local_npes_model(i)
        if (mype_ens == pe_index) then
           task_id = i
           exit doens1
        end if
        pe_index = pe_index + 1
     end do
  end do doens1


  call MPI_Comm_split(COMM_ensemble, task_id, mype_ens, &
       COMM_model, MPIerr)
  
  ! *** Re-initialize Process informations   ***
  ! *** according to model communicator ***
  call MPI_Comm_Size(COMM_model, npes_model, MPIerr)
  call MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

  if (screen > 1) then
    write (*,*) 'model-PDAF: mype(w)= ', mype_ens, '; model task: ', task_id, &
         '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
  end if


  ! Init flag FILTERProcess (all Processs of model task 1)
  if (task_id == 1) then
     filterpe = .true.
  else
     filterpe = .false.
  end if

  ! ***         COMM_FILTER                 ***
  ! *** Generate communicator for filter    ***
  ! *** For simplicity equal to COMM_couple ***

  if (filterpe) then
     my_color = task_id
  else
     my_color = MPI_UNDEFINED
  endif

  call MPI_Comm_split(COMM_ensemble, my_color, mype_ens, &
       COMM_filter, MPIerr)

  ! *** Initialize Process informations         ***
  ! *** according to coupling communicator ***
  if (filterpe) then
     call MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
     call MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)
  endif


  ! ***              COMM_COUPLE                 ***
  ! *** Generate communicators for communication ***
  ! *** between model and filter Processs             ***
  ! *** (Split COMM_ENSEMBLE)                    ***

  color_couple = mype_model + 1

  call MPI_Comm_split(COMM_ensemble, color_couple, mype_ens, &
       COMM_couple, MPIerr)

  ! *** Initialize Process informations         ***
  ! *** according to coupling communicator ***
  call MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
  call MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

  if (screen > 0) then
     if (mype_ens == 0) then
        write (*, '(/a, 2x, a)') 'model-PDAF Pconf', 'Process configuration:'
        write (*, '(a, 2x, a6, a9, a10, a14, a13, /a, 2x, a5, a9, a7, a7, a7, a7, a7, /a, 2x, a)') &
             'model-PDAF Pconf', 'world', 'filter', 'model', 'couple', 'filterPE', &
             'model-PDAF Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
             'model-PDAF Pconf', '----------------------------------------------------------'
     end if
     call MPI_Barrier(COMM_ensemble, MPIerr)
     if (task_id == 1) then
        write (*, '(a, 2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
             'model-PDAF Pconf', mype_ens, mype_filter, task_id, mype_model, color_couple, &
             mype_couple, filterpe
     endif
     if (task_id > 1) then
        write (*,'(a, 2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
         'model-PDAF Pconf', mype_ens, task_id, mype_model, color_couple, mype_couple, filterpe
     end if
     call MPI_Barrier(COMM_ensemble, MPIerr)

     if (mype_ens == 0) write (*, '(/a)') ''

  end if


! ***************************************************
! *** Provide parallelization information to PDAF ***
! ***************************************************

  call PDAF3_set_parallel(COMM_ensemble, COMM_model, COMM_filter, COMM_couple, &
       task_id, n_modeltasks, filterpe, flag)


! ******************************************************************************
! *** Initialize model equivalents to COMM_model, npes_model, and mype_model ***
! ******************************************************************************

  ! If the names of the variables for COMM_model, npes_model, and 
  ! mype_model are different in the numerical model, the 
  ! model-internal variables should be initialized at this point.

  ! THIS STEP IS NOT DONE WITH OFFLINE COUPLING

end subroutine init_parallel_pdaf_offline

