subroutine diag_reduction (diag)

    use mpi

    implicit none

    include "parameter.h"
    include "CB_diagnostics.h"
    include "CB_mpi.h"

    character(*), intent(in) :: diag

    ! for stresses
    if (trim(diag) == "stress") then
        call mpi_allreduce( &
        sigxx, sigxx_r, n, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)
        
        call mpi_allreduce( &
        sigyy, sigyy_r, n, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)

        call mpi_allreduce( &
        sigxy, sigxy_r, n, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)

        call mpi_allreduce( &
        sigyx, sigyx_r, n, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)
    end if

    ! ! same for pressure
    ! if (trim(diag) == "pressure") then
        
    ! end if

end subroutine diag_reduction


subroutine diag_broadcast (diag)

    use mpi

    implicit none

    include "parameter.h"
    include "CB_diagnostics.h"
    include "CB_mpi.h"

    character(*), intent(in) :: diag

    ! for stresses
    if (trim(diag) == "stress") then
        call mpi_allreduce( &
        tsigxx_r, tsigxx, n, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)
        
        call mpi_allreduce( &
        tsigyy_r, tsigyy, n, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)

        call mpi_allreduce( &
        tsigxy_r, tsigxy, n, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)

        call mpi_allreduce( &
        tsigyx_r, tsigyx, n, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)
    end if

    ! ! same for pressure
    ! if (trim(diag) == "pressure") then

    ! end if

end subroutine