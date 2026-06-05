module mod_density
  !!
  !! Fast density reduction / charge-density builder
  !!
  !! Main performance change vs the previous modular version:
  !!   - No full-size temporary np_work allocation/copy every timestep.
  !!   - Periodic stitching is applied directly to np_thread, exactly where the
  !!     old reduce_species_density already did it.
  !!   - build_rho_from_np now reads np_red directly instead of copying it.
  !!
  !! This keeps the same numerical convention as your current modular code:
  !!   reduce_species_density:
  !!      np_thread -> periodic density stitching -> np_red
  !!   build_rho_from_np:
  !!      rho = - sum_s charge(s) * np_red(s)
  !!
  use iso_fortran_env, only: real64, int32
  use mpi
  use mod_constants, only: qe, eps0

  implicit none
  private

  public :: reduce_species_density
  public :: build_rho_from_np
  public :: build_rho_from_np_thread
  public :: density_max_per_species

contains

  subroutine reduce_species_density(n, bcnd, np_thread, ntype, nproc, mpi_comm, np_red)
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc, mpi_comm
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
    real(real64),   intent(out)   :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)

    integer :: ierr, mpi_size
    integer :: ix, iy, iz, ptype, iproc
    real(real64) :: acc

    call MPI_Comm_size(mpi_comm, mpi_size, ierr)

    ! Important: this is intentionally in-place.
    ! In your timestep, np_thread is rebuilt by deposition before this is called,
    ! so modifying its periodic/ghost planes here is equivalent to the old
    ! np_work = np_thread copy, but avoids the large allocation/copy.
    call apply_periodic_density_bc(n, bcnd, np_thread, ntype, nproc)

    np_red = 0.0_real64

    !$omp parallel do collapse(4) private(iproc,acc) schedule(static) default(shared)
    do ptype = 1, ntype
      do iz = 1, n(3)+1
        do iy = 1, n(2)+1
          do ix = 1, n(1)+1
            acc = 0.0_real64
            do iproc = 1, nproc
              acc = acc + np_thread(ix,iy,iz,ptype,iproc)
            end do
            np_red(ix,iy,iz,ptype) = acc
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (mpi_size > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, np_red, &
          (n(1)+3)*(n(2)+3)*(n(3)+3)*ntype, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ierr)
    end if

  end subroutine reduce_species_density


  subroutine build_rho_from_np(n, np_red, charge, ntype, rho, bcnd, flag_pbc)
    integer(int32), intent(in)  :: n(3)
    integer,        intent(in)  :: ntype
    real(real64),   intent(in)  :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(in)  :: charge(ntype)
    real(real64),   intent(out) :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)  :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)  :: flag_pbc

    integer :: ix, iy, iz, ptype
    real(real64) :: r

    ! Keep current modular behavior: no second periodic correction here.
    ! The previous file had that correction commented out after copying np_red
    ! into np_work. Therefore we read np_red directly.
    rho = 0.0_real64

    !$omp parallel do collapse(3) private(ptype,r) schedule(static) default(shared)
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          r = 0.0_real64
          do ptype = 1, ntype
            r = r - charge(ptype) * np_red(ix,iy,iz,ptype)
          end do
          rho(ix,iy,iz) = r
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine build_rho_from_np


  subroutine build_rho_from_np_thread(n, np_thread, charge, ntype, nproc, rho, bcnd, flag_pbc)
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)
    real(real64),   intent(in)    :: charge(ntype)
    real(real64),   intent(out)   :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)    :: flag_pbc

    integer :: ix, iy, iz, ptype, iproc
    real(real64) :: r

    ! This routine is not used in your current mod_simulation timestep path,
    ! but it is kept as a fast drop-in equivalent.
    if (flag_pbc == 1_int32) then
      call apply_periodic_density_bc(n, bcnd, np_thread, ntype, nproc)
    end if

    rho = 0.0_real64

    !$omp parallel do collapse(3) private(iproc,ptype,r) schedule(static) default(shared)
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          r = 0.0_real64
          do iproc = 1, nproc
            do ptype = 1, ntype
              r = r - charge(ptype) * np_thread(ix,iy,iz,ptype,iproc)
            end do
          end do
          rho(ix,iy,iz) = r
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine build_rho_from_np_thread


  subroutine density_max_per_species(n, bcnd, np_red, ntype, np_mx)
    integer(int32), intent(in)  :: n(3)
    integer,        intent(in)  :: ntype
    integer,        intent(in)  :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(in)  :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(out) :: np_mx(ntype)

    integer :: ix, iy, iz, ptype
    real(real64) :: local_max

    np_mx = 0.0_real64

    do ptype = 1, ntype
      local_max = 0.0_real64

      !$omp parallel do collapse(3) reduction(max:local_max) schedule(static) default(shared)
      do iz = 1, n(3)+1
        do iy = 1, n(2)+1
          do ix = 1, n(1)+1
            local_max = max(local_max, np_red(ix,iy,iz,ptype))
          end do
        end do
      end do
      !$omp end parallel do

      np_mx(ptype) = local_max
    end do

  end subroutine density_max_per_species


  subroutine apply_periodic_density_bc(n, bcnd, np_thread, ntype, nproc)
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype,nproc)

    integer :: ix, iy, iz, iproc, ptype

    ! ------------------------------------------------------------
    ! Legacy periodic density stitching in y.
    ! Written with explicit ptype loop to avoid array-section temporaries.
    ! ------------------------------------------------------------
    !$omp parallel do collapse(3) private(ptype) schedule(static) default(shared)
    do iproc = 1, nproc
      do iz = 1, n(3)+1
        do ix = 1, n(1)+1
          if (bcnd(ix,1,iz) == 0) then
            do ptype = 1, ntype
              np_thread(ix,1,iz,ptype,iproc) = 0.5_real64 * ( &
                   np_thread(ix,1,iz,ptype,iproc) + np_thread(ix,n(2)+1,iz,ptype,iproc) )
              np_thread(ix,0,iz,ptype,iproc) = np_thread(ix,n(2),iz,ptype,iproc)
            end do
          end if

          if (bcnd(ix,n(2)+1,iz) == 0) then
            do ptype = 1, ntype
              np_thread(ix,n(2)+2,iz,ptype,iproc) = np_thread(ix,2,iz,ptype,iproc)
              np_thread(ix,n(2)+1,iz,ptype,iproc) = np_thread(ix,1,iz,ptype,iproc)
            end do
          end if
        end do
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------
    ! Legacy periodic density stitching in z.
    ! Written with explicit ptype loop to avoid array-section temporaries.
    ! ------------------------------------------------------------
    !$omp parallel do collapse(3) private(ptype) schedule(static) default(shared)
    do iproc = 1, nproc
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          if (bcnd(ix,iy,1) == 0) then
            do ptype = 1, ntype
              np_thread(ix,iy,1,ptype,iproc) = 0.5_real64 * ( &
                   np_thread(ix,iy,1,ptype,iproc) + np_thread(ix,iy,n(3)+1,ptype,iproc) )
              np_thread(ix,iy,0,ptype,iproc) = np_thread(ix,iy,n(3),ptype,iproc)
            end do
          end if

          if (bcnd(ix,iy,n(3)+1) == 0) then
            do ptype = 1, ntype
              np_thread(ix,iy,n(3)+2,ptype,iproc) = np_thread(ix,iy,2,ptype,iproc)
              np_thread(ix,iy,n(3)+1,ptype,iproc) = np_thread(ix,iy,1,ptype,iproc)
            end do
          end if
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine apply_periodic_density_bc

end module mod_density
