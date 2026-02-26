module mod_simulation
  use mod_state
  use mpi
  implicit none
  private
  public :: Simulation

  type :: Simulation
     type(State) :: S
     integer     :: comm = MPI_COMM_WORLD
  contains
     procedure :: init
     procedure :: run
     procedure :: finalize
  end type Simulation

contains

  subroutine init(self, comm_in)
    class(Simulation), intent(inout) :: self
    integer, intent(in) :: comm_in
    self%comm = comm_in

    call self%S%init(self%comm)
    call self%S%build_boundary_only()
  end subroutine init

  subroutine run(self, nsteps)
    class(Simulation), intent(inout) :: self
    integer, intent(in) :: nsteps
    integer :: i
    do i = 1, nsteps
      ! placeholder: later you’ll put deposit/poisson/push here
    end do
  end subroutine run

  subroutine finalize(self)
    class(Simulation), intent(inout) :: self
    call self%S%finalize()
  end subroutine finalize

end module mod_simulation