module mpi_counts_mod

    use mpi_f08, only: mpi_init, mpi_comm_rank, mpi_comm_size

    implicit none

    ! number of elements contributed by each rank
    integer, allocatable :: counts(:)

    ! displacement (in elements) in recv buffer for each rank, zero-based values
    integer, allocatable :: displs(:)

contains

    subroutine init_mpi_counts(n, n_ranks, rank, local_n, local_disp)
        ! Initialize counts(:) and displs(:) for MPI_Allgatherv usage.
        ! Arguments:
        !   n         - total length of the global vector
        !   n_ranks   - number of MPI ranks
        !   rank      - this process' rank (0..n_ranks-1)
        ! Outputs:
        !   counts    - allocatable array size n_ranks, counts(i) is number of elements from rank (i-1)
        !   displs    - allocatable array size n_ranks, displs(i) is zero-based offset in elements for rank (i-1)
        !   local_n   - number of elements this rank will send (last-first+1)
        !   local_disp- zero-based offset in recv buffer for this rank (== displs(rank+1))

        implicit none

        integer, intent(in) :: n, n_ranks, rank
        integer, intent(out) :: local_n, local_disp
        integer :: base, rem, i
        integer :: sum

        if (.not. allocated(counts)) then
            allocate(counts(n_ranks))
        end if

        if (.not. allocated(displs)) then
            allocate(displs(n_ranks))
        end if

        ! Balanced partition: first rem ranks get base+1, others get base
        base = n / max(1, n_ranks)
        rem  = mod(n, max(1, n_ranks))

        do i = 1, n_ranks
            if (i-1 < rem) then
                counts(i) = base + 1
            else
                counts(i) = base
            end if
        end do

        ! compute zero-based displacements in units of elements
        sum = 0
        do i = 1, n_ranks
            displs(i) = sum       ! zero-based
            sum = sum + counts(i)
        end do

        ! protective: if total != n, adjust last element
        if (sum /= n) then
            ! adjust last count so sum == n
            counts(n_ranks) = counts(n_ranks) + (n - sum)
        end if

        ! return the local values for this rank
        if (rank >= 0 .and. rank < n_ranks) then
            local_n = counts(rank+1)
            local_disp = displs(rank+1)
        else
            local_n = 0
            local_disp = 0
        end if

    end subroutine init_mpi_counts

end module mpi_counts_mod