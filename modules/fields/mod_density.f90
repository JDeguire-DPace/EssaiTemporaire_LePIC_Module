module mod_density
  !!
  !! Modern density reduction / charge-density builder
  !!
  !! Replaces the legacy ideas of:
  !!   - dens_red
  !!   - calc_rho   (partially: rho build only; slab extraction can come later)
  !!
  !! Main routines:
  !!   1) reduce_species_density:
  !!        thread-local np_thread(:,:,:,:,iproc)  ->  np_red(:,:,:,ptype)
  !!        with legacy periodic stitching in y/z before reduction
  !!
  !!   2) build_rho_from_np:
  !!        np_red(:,:,:,ptype) -> rho(:,:,:)
  !!
  !!   3) density_max_per_species:
  !!        compute max density in the simulation domain for each species
  !!
  use iso_fortran_env, only: real64, int32
  use mpi
  implicit none
  private

  public :: reduce_species_density
  public :: build_rho_from_np
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

    call MPI_Comm_size(mpi_comm, mpi_size, ierr)

    np_red = 0.0_real64

    ! ------------------------------------------------------------
    ! Legacy periodic boundary stitching on thread-local densities
    ! (same logic as calc_rho before reduction)
    ! ------------------------------------------------------------
    call apply_periodic_density_bc(n, bcnd, np_thread, ntype, nproc)

    ! print *, 'THREAD nx-1 sum p1/p2 = ', &
    !   sum(np_thread(n(1)-1,:,:,1,:)), sum(np_thread(n(1)-1,:,:,2,:))

    ! print *, 'THREAD nx sum p1/p2 = ', &
    !   sum(np_thread(n(1),:,:,1,:)), sum(np_thread(n(1),:,:,2,:))

    ! print *, 'THREAD nx+1 sum p1/p2 = ', &
    !   sum(np_thread(n(1)+1,:,:,1,:)), sum(np_thread(n(1)+1,:,:,2,:))

    ! ------------------------------------------------------------
    ! OpenMP-thread reduction: np_thread -> np_red
    ! Legacy loops only over 1:n+1 active nodes
    ! ------------------------------------------------------------
    !$omp parallel do collapse(3) private(iproc,ptype) default(shared)
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          do iproc = 1, nproc
            do ptype = 1, ntype
              np_red(ix,iy,iz,ptype) = np_red(ix,iy,iz,ptype) + np_thread(ix,iy,iz,ptype,iproc)
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------
    ! MPI reduction across ranks
    ! ------------------------------------------------------------
    if (mpi_size > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, np_red, &
           (n(1)+3)*(n(2)+3)*(n(3)+3)*ntype, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ierr)
    end if

    ! print *, 'CHECK nx+1 np1 = ', sum(np_red(n(1)+1,:,: ,1))
    ! print *, 'CHECK nx+1 np2 = ', sum(np_red(n(1)+1,:,: ,2))

  end subroutine reduce_species_density


  subroutine build_rho_from_np(n, np_red, charge, ntype, rho)
    integer(int32), intent(in)  :: n(3)
    integer,        intent(in)  :: ntype
    real(real64),   intent(in)  :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(in)  :: charge(ntype)
    real(real64),   intent(out) :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    integer :: ix, iy, iz, ptype

    rho = 0.0_real64

    !$omp parallel do collapse(3) private(ptype) default(shared)
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          do ptype = 1, ntype
            rho(ix,iy,iz) = rho(ix,iy,iz) - charge(ptype) * np_red(ix,iy,iz,ptype)
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine build_rho_from_np


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

      !$omp parallel do collapse(3) reduction(max:local_max) default(shared)
      do iz = 1, n(3)+1
        do iy = 1, n(2)+1
          do ix = 1, n(1)+1
            if (bcnd(ix,iy,iz) == -1) then
              local_max = max(local_max, np_red(ix,iy,iz,ptype))
            end if
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

    integer :: ix, iy, iz, iproc

    ! ------------------------------------------------------------
    ! Legacy periodic boundary logic in y
    ! ------------------------------------------------------------
    !$omp parallel do collapse(2) private(iproc) default(shared)
    do iz = 1, n(3)+1
      do ix = 1, n(1)+1
        do iproc = 1, nproc
          if (bcnd(ix,1,iz) == 0) then
            np_thread(ix,1,iz,:,iproc) = 0.5_real64 * ( &
                 np_thread(ix,1,iz,:,iproc) + np_thread(ix,n(2)+1,iz,:,iproc) )
            np_thread(ix,0,iz,:,iproc) = np_thread(ix,n(2),iz,:,iproc)
          end if

          if (bcnd(ix,n(2)+1,iz) == 0) then
            np_thread(ix,n(2)+2,iz,:,iproc) = np_thread(ix,2,iz,:,iproc)
            np_thread(ix,n(2)+1,iz,:,iproc) = np_thread(ix,1,iz,:,iproc)
          end if
        end do
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------
    ! Legacy periodic boundary logic in z
    ! ------------------------------------------------------------
    !$omp parallel do collapse(2) private(iproc) default(shared)
    do iy = 1, n(2)+1
      do ix = 1, n(1)+1
        do iproc = 1, nproc
          if (bcnd(ix,iy,1) == 0) then
            np_thread(ix,iy,1,:,iproc) = 0.5_real64 * ( &
                 np_thread(ix,iy,1,:,iproc) + np_thread(ix,iy,n(3)+1,:,iproc) )
            np_thread(ix,iy,0,:,iproc) = np_thread(ix,iy,n(3),:,iproc)
          end if

          if (bcnd(ix,iy,n(3)+1) == 0) then
            np_thread(ix,iy,n(3)+2,:,iproc) = np_thread(ix,iy,2,:,iproc)
            np_thread(ix,iy,n(3)+1,:,iproc) = np_thread(ix,iy,1,:,iproc)
          end if
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine apply_periodic_density_bc

end module mod_density