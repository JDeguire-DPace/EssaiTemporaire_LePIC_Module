module mod_electricField
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: calc_Efield_modular

contains

  subroutine calc_Efield_modular(n, h, phi, E, bcnd)
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)
    real(real64),   intent(in) :: phi(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64), intent(inout) :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    integer :: ix, iy, iz

    E = 0.0_real64

    !
    ! Interior points
    !
    !$OMP PARALLEL
    !$OMP DO
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1

          if (bcnd(ix,iy,iz) <= 0) then

            ! Second-order centered differences
            E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
            E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
            E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))

          end if

        end do
      end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !
    ! Boundary conditions
    !
    !$OMP PARALLEL
    !$OMP DO
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1

          if (bcnd(ix,iy,iz) == -1) goto 10

          !
          ! Neumann BCs
          !
          if (bcnd(ix,iy,iz) == -2) then
            ! YZ plane, LHS only

            E(1,1,iy,iz) = 0.0_real64

            if (iy > 1 .and. iy < n(2)+1 .and. &
                iz > 1 .and. iz < n(3)+1) then
              E(2,1,iy,iz) = -(phi(1,iy+1,iz) - phi(1,iy-1,iz)) / (2.0_real64*h(2))
              E(3,1,iy,iz) = -(phi(1,iy,iz+1) - phi(1,iy,iz-1)) / (2.0_real64*h(3))
              goto 10
            end if
          end if

          !
          ! Periodic BCs
          !
          if (bcnd(ix,iy,iz) == 0 .or. bcnd(ix,iy,iz) == -2) then

            ! iy = 0,1,ny+1,ny+2 planes
            if (iy == 1) then
              E(2,ix,1,iz)       = -(phi(ix,2,iz) - phi(ix,n(2),iz)) / (2.0_real64*h(2))
              E(2,ix,0,iz)       = E(2,ix,n(2),iz)
              E(2,ix,n(2)+1,iz)  = E(2,ix,1,iz)
              E(2,ix,n(2)+2,iz)  = E(2,ix,2,iz)

              E(1,ix,0,iz)       = E(1,ix,n(2),iz)
              E(1,ix,n(2)+2,iz)  = E(1,ix,2,iz)

              E(3,ix,0,iz)       = E(3,ix,n(2),iz)
              E(3,ix,n(2)+2,iz)  = E(3,ix,2,iz)
            end if

            ! iz = 0,1,nz+1,nz+2 planes
            if (iz == 1) then
              E(3,ix,iy,1)       = -(phi(ix,iy,2) - phi(ix,iy,n(3))) / (2.0_real64*h(3))
              E(3,ix,iy,0)       = E(3,ix,iy,n(3))
              E(3,ix,iy,n(3)+1)  = E(3,ix,iy,1)
              E(3,ix,iy,n(3)+2)  = E(3,ix,iy,2)

              E(1,ix,iy,0)       = E(1,ix,iy,n(3))
              E(1,ix,iy,n(3)+2)  = E(1,ix,iy,2)

              E(2,ix,iy,0)       = E(2,ix,iy,n(3))
              E(2,ix,iy,n(3)+2)  = E(2,ix,iy,2)
            end if

            goto 10
          end if

          !
          ! Walls
          !
          if (bcnd(ix,iy,iz) >= 1) then

            ! West wall
            if (bcnd(ix+1,iy,iz) <= 0) then
              E(1,ix,iy,iz) = 2.0_real64*E(1,ix+1,iy,iz) - E(1,ix+2,iy,iz)
              E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
              E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))
            end if

            ! East wall
            if (bcnd(ix-1,iy,iz) <= 0) then
              E(1,ix,iy,iz) = 2.0_real64*E(1,ix-1,iy,iz) - E(1,ix-2,iy,iz)
              E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
              E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))
            end if

            ! South wall
            if (bcnd(ix,iy+1,iz) <= 0) then
              E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
              E(2,ix,iy,iz) = 2.0_real64*E(2,ix,iy+1,iz) - E(2,ix,iy+2,iz)
              E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))
            end if

            ! North wall
            if (bcnd(ix,iy-1,iz) <= 0) then
              E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
              E(2,ix,iy,iz) = 2.0_real64*E(2,ix,iy-1,iz) - E(2,ix,iy-2,iz)
              E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))
            end if

            ! Bottom wall
            if (bcnd(ix,iy,iz+1) <= 0) then
              E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
              E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
              E(3,ix,iy,iz) = 2.0_real64*E(3,ix,iy,iz+1) - E(3,ix,iy,iz+2)
            end if

            ! Top wall
            if (bcnd(ix,iy,iz-1) <= 0) then
              E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
              E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
              E(3,ix,iy,iz) = 2.0_real64*E(3,ix,iy,iz-1) - E(3,ix,iy,iz-2)
            end if

          end if

10        continue

        end do
      end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

  end subroutine calc_Efield_modular

end module mod_electricField
