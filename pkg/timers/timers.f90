!=======================================================================
!   Module: wall-clock timers of the phases of a time step
!
! timer_on(k) / timer_off(k) accumulate the time spent in phase k over
! the run; timer_report prints the totals at the end (minimum and
! maximum over the MPI ranks). A dozen calls per time step, so their
! cost is negligible. Self-contained (no model module), as all of pkg/.
!=======================================================================

module timers

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_wtime

    implicit none

    ! phases of a time step
    integer, parameter ::       &
        t_integ   = 1,          & ! time integration (Verlet, ridging)
        t_tree    = 2,          & ! kd-tree and neighbour lists
        t_pairs   = 3,          & ! force loop over the particles
        t_forcing = 4,          & ! winds, currents and Coriolis
        t_comm    = 5,          & ! MPI reductions and broadcasts
        t_output  = 6,          & ! outputs (with the bond gathering)
        n_timers  = 6

    character(len=8), parameter :: timer_names(n_timers) = [ &
        character(len=8) :: "integ", "tree", "pairs", "forcing", &
        "comm", "output" ]

    real(dp) :: timer_total(n_timers) = 0.0_dp
    real(dp) :: timer_start(n_timers) = 0.0_dp

contains

    subroutine timer_on(k)

        integer, intent(in) :: k

        timer_start(k) = omp_get_wtime()

    end subroutine timer_on


    subroutine timer_off(k)

        integer, intent(in) :: k

        timer_total(k) = timer_total(k) + omp_get_wtime() - timer_start(k)

    end subroutine timer_off


    subroutine timer_report(loop_time)

        ! prints the time of each phase on rank 0, with the minimum and
        ! maximum over the ranks (a large difference in comm means that
        ! the ranks wait for each other: load imbalance)

        use mpi_f08

        real(dp), intent(in) :: loop_time

        real(dp) :: tmin(n_timers), tmax(n_timers)
        integer :: k, rank, ierr

        call mpi_comm_rank(mpi_comm_world, rank, ierr)

        call mpi_reduce(timer_total, tmin, n_timers, mpi_double_precision, &
                        mpi_min, 0, mpi_comm_world, ierr)
        call mpi_reduce(timer_total, tmax, n_timers, mpi_double_precision, &
                        mpi_max, 0, mpi_comm_world, ierr)

        if ( rank .eq. 0 ) then
            print '(a)', ' Timers (s, rank 0, min and max over ranks):'
            do k = 1, n_timers
                print '(a, a8, a, f12.3, a, f6.1, a, f12.3, f12.3)',     &
                    ' Timer ', timer_names(k), ':', timer_total(k),     &
                    ' (', 100 * timer_total(k) / max(loop_time, 1e-30_dp), &
                    ' %)', tmin(k), tmax(k)
            end do
        end if

    end subroutine timer_report

end module timers
