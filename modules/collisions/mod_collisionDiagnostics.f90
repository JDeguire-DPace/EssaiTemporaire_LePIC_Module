module mod_collisionDiagnostics
  use iso_fortran_env, only: int32, int64
  implicit none
  private
  public :: init_rxn_counts, count_rxn, print_rxn_counts, reset_rxn_counts, count_mcc, count_bmcc

  integer(int64), allocatable, save :: rxn_count(:)
  integer(int64), save :: n_mcc = 0_int64
  integer(int64), save :: n_bmcc = 0_int64

contains

  subroutine init_rxn_counts(ncol)
    integer(int32), intent(in) :: ncol

    if (.not. allocated(rxn_count)) then
      allocate(rxn_count(ncol))
      rxn_count = 0_int64
    end if
  end subroutine init_rxn_counts

  subroutine count_rxn(c_ind)
    integer(int32), intent(in) :: c_ind
    if (.not. allocated(rxn_count)) return
    if (c_ind < 1 .or. c_ind > size(rxn_count)) return
    !$omp atomic
    rxn_count(c_ind) = rxn_count(c_ind) + 1_int64
  end subroutine count_rxn

  subroutine print_rxn_counts()
    integer :: i
    if (.not. allocated(rxn_count)) return
    write(*,'(a)') " Accepted collision reactions:"
    do i = 1, size(rxn_count)
      if (rxn_count(i) > 0_int64) then
        write(*,'(a,i4,a,i12)') " reaction ", i, " accepted = ", rxn_count(i)
      end if
    end do
  end subroutine print_rxn_counts

  subroutine reset_rxn_counts()
    if (allocated(rxn_count)) rxn_count = 0_int64
  end subroutine reset_rxn_counts

  subroutine count_mcc(c_ind)
    integer(int32), intent(in) :: c_ind
    call count_rxn(c_ind)
    !$omp atomic
    n_mcc = n_mcc + 1_int64
  end subroutine

  subroutine count_bmcc(c_ind)
    integer(int32), intent(in) :: c_ind
    call count_rxn(c_ind)
    !$omp atomic
    n_bmcc = n_bmcc + 1_int64
  end subroutine

end module mod_collisionDiagnostics