module mod_collision_flat_kernel
  use iso_fortran_env, only: int32, real64, int8
  use mpi
  use mod_RNG, only: ran2
  use mod_utils, only: stop_calculation
  use mod_collision_flat_backend, only: FlatCollisionBuffers

  implicit none
  private
  public :: perform_flat_collisions_step

  integer(int32), parameter :: MAX_COL = 128_int32

contains

  subroutine perform_flat_collisions_step( &
      flat, n, h, mass, p_ncol, sig_list, col_info, sigv_mx, &
      np_mx_all, ns_coll, dt, nu_uplim, iseed, nproc_mpi, mpi_rank)

    type(FlatCollisionBuffers), intent(inout) :: flat
    real(real64), intent(in) :: mass(:)
    integer(int32), intent(in) :: p_ncol(:)
    integer(int32), intent(in) :: sig_list(:,:), col_info(:,:)
    real(real64), intent(in) :: sigv_mx(:,:)
    real(real64), intent(in) :: np_mx_all(:)
    integer(int32), intent(in) :: ns_coll
    real(real64), intent(in) :: dt, nu_uplim(:)
    integer(int32), intent(inout) :: iseed(:)
    integer(int32), intent(in) :: nproc_mpi, mpi_rank

    integer(int32) :: ptype, iproc, icol, ind_col, n_re, ttype
    integer(int32) :: sum_np_tot_local, sum_np_tot_global
    integer(int32) :: Nc_tmp, ierr
		integer(int32) :: ic, ip, chosen_icol
		real(real64) :: csum
    integer(int32) :: accepted
    real(real64) :: nu(MAX_COL), sort_arr(MAX_COL), sum_nu
    integer(int32) :: indx(MAX_COL)
    logical :: flag_coll
    real(real64) :: nu_max, Pmax, dNc, rnd
    real(real64) :: vx1, vy1, vz1
    real(real64) :: vx2, vy2, vz2
    real(real64) :: vr
    integer(int32) :: icell, i0, i1, itarget
    real(real64) :: np_t
		integer(int32), intent(in) :: n(3)
		real(real64), intent(in) :: h(3)
		integer(int32) :: Nc_local

    do ptype = 1, flat%ntype
      if (p_ncol(ptype) == 0_int32) cycle
      if (mass(ptype) <= 0.0_real64) cycle

      nu_max = 0.0_real64

      do icol = 1, min(p_ncol(ptype), MAX_COL)
        ind_col = sig_list(ptype,icol)
        n_re    = col_info(ind_col,1)
        ttype   = col_info(ind_col,2+n_re)

        if (ttype < 1 .or. ttype > size(np_mx_all)) cycle

        nu_max = nu_max + np_mx_all(ttype) * sigv_mx(ptype,ind_col)
      end do

      if (ptype <= size(nu_uplim)) then
        nu_max = min(nu_max, nu_uplim(ptype))
      end if

      sum_np_tot_local = sum(flat%np_tot(ptype,1:flat%nproc))

      if (nproc_mpi > 1_int32) then
        call MPI_Allreduce(sum_np_tot_local, sum_np_tot_global, 1, &
             MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      else
        sum_np_tot_global = sum_np_tot_local
      end if

      Pmax = nu_max * real(ns_coll, real64) * dt

      if (Pmax > 1.0_real64) then
        if (mpi_rank == 0) then
          write(*,*) "flat collision probability > 1"
          write(*,*) "ptype, Pmax = ", ptype, Pmax
        end if
        call stop_calculation
      end if

      dNc = real(sum_np_tot_global, real64) * Pmax / &
            real(max(1_int32,nproc_mpi*flat%nproc), real64)

      Nc_tmp = int(dNc, int32)
      rnd = ran2(iseed(1))
      if (rnd <= dNc - real(Nc_tmp, real64)) Nc_tmp = Nc_tmp + 1_int32

 

        if (mpi_rank == 0) then
        write(*,*) "FLAT COLL ptype=", ptype, &
                    " nu_max=", nu_max, &
                    " Pmax=", Pmax, &
                    " dNc=", dNc, &
                    " Nc_tmp=", Nc_tmp
        end if

        accepted = 0_int32

				do iproc = 1_int32, flat%nproc

					if (flat%np_tot(ptype,iproc) <= 0_int32) cycle

					Nc_local = int( &
							real(flat%np_tot(ptype,iproc),real64) / &
							real(sum_np_tot_local,real64) * &
							real(Nc_tmp,real64), &
							int32)

					do ic = 1_int32, Nc_local

            rnd = ran2(iseed(iproc))
            ip = int(real(flat%np_tot(ptype,iproc), real64) * rnd, int32) + 1_int32
            if (ip > flat%np_tot(ptype,iproc)) ip = flat%np_tot(ptype,iproc)

            if (flat%flag_dead(ip,ptype,iproc) /= 0_int8) cycle

            sum_nu = 0.0_real64
            nu(:)  = 0.0_real64

            ! Projectile velocity
            vx1 = flat%vxp(4,ip,ptype,iproc)
            vy1 = flat%vxp(5,ip,ptype,iproc)
            vz1 = flat%vxp(6,ip,ptype,iproc)

						icell = flat_particle_cell_id( &
						flat%vxp(1,ip,ptype,iproc), &
						flat%vxp(2,ip,ptype,iproc), &
						flat%vxp(3,ip,ptype,iproc), &
						n, h)

            do icol = 1, min(p_ncol(ptype), MAX_COL)

            ind_col = sig_list(ptype,icol)

            n_re  = col_info(ind_col,1)
            ttype = col_info(ind_col,2+n_re)

            ! Select a target particle from the same cell if target is tracked.
            if (ttype <= flat%ntype) then

            i0 = flat%Plist(icell-1,ttype,iproc) + 1_int32
            i1 = flat%Plist(icell  ,ttype,iproc)

            if (i1 < i0) cycle

            rnd = ran2(iseed(iproc))
            itarget = i0 + int(real(i1-i0+1_int32, real64) * rnd, int32)
            if (itarget > i1) itarget = i1

            vx2 = flat%vxp(4,itarget,ttype,iproc)
            vy2 = flat%vxp(5,itarget,ttype,iproc)
            vz2 = flat%vxp(6,itarget,ttype,iproc)

            np_t = np_mx_all(ttype)

            else

            ! Neutral target: for now assume cold background gas.
            vx2 = 0.0_real64
            vy2 = 0.0_real64
            vz2 = 0.0_real64

            np_t = np_mx_all(ttype)

            end if

            vr = sqrt((vx1-vx2)**2 + &
                    (vy1-vy2)**2 + &
                    (vz1-vz2)**2)

            nu(icol) = np_t * sigv_mx(ptype,ind_col)


            

            sum_nu = sum_nu + nu(icol)

            end do

            if (sum_nu <= 0.0_real64) cycle

            rnd = ran2(iseed(iproc)) * nu_max

            if (rnd > sum_nu) cycle

            accepted = accepted + 1_int32

            ! Reaction selection
						csum = 0.0_real64
						do icol = 1, min(p_ncol(ptype), MAX_COL)
							csum = csum + nu(icol)

						if (rnd <= csum) then
							chosen_icol = icol
							ind_col = sig_list(ptype, chosen_icol)

							! For now: only elastic reactions.
							! Elastic = same projectile and same target, no species creation.
							n_re  = col_info(ind_col,1)
							ttype = col_info(ind_col,2+n_re)

							if (ttype <= flat%ntype) then
								! Simple velocity swap / scattering placeholder.
								! This validates that flat%vxp can be modified safely.

								flat%vxp(4,ip,ptype,iproc) = vx2
								flat%vxp(5,ip,ptype,iproc) = vy2
								flat%vxp(6,ip,ptype,iproc) = vz2

								flat%vxp(4,itarget,ttype,iproc) = vx1
								flat%vxp(5,itarget,ttype,iproc) = vy1
								flat%vxp(6,itarget,ttype,iproc) = vz1
							end if

							exit
						end if
						end do

            if (accepted <= 5 .and. mpi_rank == 0) then
                write(*,*) "ptype=", ptype, &
                            " chosen reaction=", chosen_icol
            end if

        end do
        end do

        if (mpi_rank == 0) then
        write(*,*) "FLAT ACCEPT TEST ptype=", ptype, " accepted=", accepted
        end if

    end do

  end subroutine perform_flat_collisions_step

  pure integer(int32) function flat_particle_cell_id(x, y, z, n, h) result(icell)
    real(real64), intent(in) :: x, y, z
    integer(int32), intent(in) :: n(3)
    real(real64), intent(in) :: h(3)

    integer(int32) :: ix, iy, iz

    ix = int(x / h(1), int32) + 1_int32
    iy = int(y / h(2), int32) + 1_int32
    iz = int(z / h(3), int32) + 1_int32

    ix = max(1_int32, min(n(1), ix))
    iy = max(1_int32, min(n(2), iy))
    iz = max(1_int32, min(n(3), iz))

    icell = (ix - 1_int32) + n(1) * ((iy - 1_int32) + n(2) * (iz - 1_int32)) + 1_int32
  end function flat_particle_cell_id

end module mod_collision_flat_kernel