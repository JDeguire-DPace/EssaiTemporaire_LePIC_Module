program main
  use mpi
  use mod_intro,      only: print_introduction
  use mod_simulation, only: Simulation
  implicit none

  type(Simulation) :: sim
  integer :: ierr

  call print_introduction()
  call MPI_Init(ierr)

  call sim%init(MPI_COMM_WORLD)
  call sim%run(150005)
  call sim%finalize()

  call MPI_Finalize(ierr)
end program main
