module mod_utils
    implicit none
    private
    public :: stop_calculation
contains
    subroutine stop_calculation
        use mpi
        implicit none
        integer ierr

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        STOP
        call MPI_Finalize(ierr)

        return
    end subroutine stop_calculation
end module mod_utils