program main
  use mpi
  use mod_state, only: State
  implicit none

  type(State) :: S
  integer :: ierr

  call MPI_Init(ierr)

  call S%init(MPI_COMM_WORLD)
  call S%build_boundary_only()

  call S%finalize()

  call MPI_Finalize(ierr)
end program main