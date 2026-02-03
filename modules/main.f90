program main
  use mpi
  use mod_state, only: State
  implicit none

  integer :: ierr
  type(State) :: st

  call MPI_Init(ierr)

  call st%init(MPI_COMM_WORLD)
  call st%build_boundary_only()
  call st%finalize()

  call MPI_Finalize(ierr)
end program main
