program main
  use mpi
  use mod_state, only: State
  implicit none

  type(State) :: S
  integer :: ierr

  call MPI_Init(ierr)

  call S%init(MPI_COMM_WORLD)
  call S%build_boundary_only()

  ! At this point we have:
  ! - cfg loaded
  ! - dom built
  ! - chemistry + reactions loaded
  ! - fields allocated/zeroed (from State%init)

  call S%finalize()

  call MPI_Finalize(ierr)
end program main