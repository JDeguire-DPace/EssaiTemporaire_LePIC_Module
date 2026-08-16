module mod_sor_timing
  ! Isolates the cost of the coarse-grid "concatenate to every rank"
  ! Allgatherv block inside sor_rb (modules/legacy/sors.f90, flag_c==1)
  ! from the rest of the multigrid solve, to check whether it is what
  ! dominates the "solve(ms)" figure in print_poisson_breakdown.
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: accum_concat_time, print_sor_timing, reset_sor_timing

  real(real64), save :: t_concat = 0.0_real64
  integer(int32), save :: n_concat = 0_int32

contains

  subroutine accum_concat_time(dt)
    real(real64), intent(in) :: dt
    t_concat = t_concat + dt
    n_concat = n_concat + 1_int32
  end subroutine accum_concat_time

  subroutine print_sor_timing(t_count)
    integer(int32), intent(in) :: t_count
    real(real64) :: denom
    if (t_count <= 0_int32) return
    denom = real(t_count, real64)
    write(*,'(a)') " ----- sor_rb coarse-grid concat (Allgatherv) timing -----"
    write(*,'(a,f8.2,a,i6)') &
      "  concat(ms)=", 1000.0_real64*t_concat/denom, &
      "  n_concat_per_call=", n_concat
    write(*,'(a)') " -----------------------------------------------------------"
  end subroutine print_sor_timing

  subroutine reset_sor_timing()
    t_concat = 0.0_real64
    n_concat = 0_int32
  end subroutine reset_sor_timing

end module mod_sor_timing
