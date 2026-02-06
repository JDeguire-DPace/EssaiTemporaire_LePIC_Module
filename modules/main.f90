program main
  use mpi
  use iso_fortran_env, only: real64
  use mod_state, only: State
  implicit none
  type(State) :: S
  integer :: ierr

  call MPI_Init(ierr)

  call S%init(MPI_COMM_WORLD)
  write(*,*) "Initialization complete. Building boundary conditions..."
  call S%build_boundary_only()

  ! --- Poisson test (Laplace) ---
  S%pdec%rhs_dom = 0.0_real64  ! vacuum
  ! call pdesolver_local(S%pdec%phi_dom, S%pdec%rhs_dom, S%pdec%bcnd_dom, ...)  ! next
  call S%pdec%gather_to_global(S%fld%phi)

  call S%finalize()
  call MPI_Finalize(ierr)
end program main
