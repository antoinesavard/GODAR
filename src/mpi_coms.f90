subroutine broadcasting (num_threads)

    use mpi

    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_mpi.h"
    include "CB_variables.h"
    include "CB_bond.h"

    integer, intent(in) :: num_threads

    ! openmp variable
    call mpi_bcast(num_threads, 1, mpi_integer,         &
                    master, mpi_comm_world, ierr)

    ! variable broadcast
    call mpi_bcast(x, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(y, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(r, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(h, n, mpi_double_precision,          &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(theta, n, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(omega, n, mpi_double_precision,      &
                    master, mpi_comm_world, ierr)
    call mpi_bcast(bond, n * n, mpi_double_precision,   &
                    master, mpi_comm_world, ierr)

end subroutine broadcasting