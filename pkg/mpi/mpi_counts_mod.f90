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

    subroutine init_mpi_counts_weighted(n, n_ranks, rank, local_n, local_disp, weights)
        ! Partition indices 1..n into contiguous blocks that equalize 
        ! the sum of a weight per row.
        ! We estimate per-row work as w(i) = (n - i) * weights(i) + eps
        ! The (n-i) factor accounts for lower triangular loop cost.
        ! The weights(i) accounts for particle size.

        implicit none

        integer, intent(in) :: n, n_ranks, rank
        integer, intent(out) :: local_n, local_disp
        double precision, intent(in), optional :: weights(:)

        integer :: i, r
        integer :: current_rank, start_idx, end_idx
        double precision :: eps
        double precision :: total_work, target, cum
        double precision, allocatable :: w(:)

        ! allocate counts/displs if needed
        if (.not. allocated(counts)) then
            allocate(counts(n_ranks))
        end if

        if (.not. allocated(displs)) then 
            allocate(displs(n_ranks))
        end if

        ! small epsilon to ensure rows with zero nominal cost are assigned at least something
        eps = 1.0d-12

        allocate(w(n))

        if (present(weights)) then
            if (size(weights) /= n) then
                write(*,*) 'ERROR init_mpi_counts_weighted: weights size mismatch'
                stop
            end if
            do i = 1, n
                ! estimated work
                w(i) = (dble(n - i) + 1.0d0) * max(eps, weights(i)) + eps
            end do
        else
            ! base work if no weights
            do i = 1, n
                w(i) = dble(n - i) + 1.0d0
            end do
        end if

        print*,'w', w
        print*,'weights', weights

        ! total work and target per rank
        total_work = sum(w)
        if (n_ranks > 0) then
            target = total_work / dble(n_ranks)
        else
            target = total_work
        end if

        print*, 'target',target
        print*, 'n_ranks',n_ranks

        ! Greedy contiguous partitioning: assign until cum >= target
        current_rank = 1
        cum = 0.0d0
        start_idx = 1
        do i = 1, n
            cum = cum + w(i)
            if ( (cum >= target .and. current_rank < n_ranks) .or. (i == n) ) then
                end_idx = i
                counts(current_rank) = end_idx - start_idx + 1
                displs(current_rank) = start_idx - 1   ! zero-based
                current_rank = current_rank + 1
                start_idx = i + 1
                ! remaining work for remaining ranks
                if (start_idx <= n) then
                    total_work = total_work - cum
                    if (current_rank <= n_ranks) then
                        target = total_work / dble(n_ranks - (current_rank - 1))
                    end if
                    cum = 0.0d0
                end if
            end if
        end do

        ! If we did not fill all ranks (e.g., more ranks than rows), 
        ! set remaining to zero
        do r = current_rank, n_ranks
            counts(r) = 0
            ! set to end-of-array
            displs(r) = n

        end do

        ! Sanity: make sure counts sum to n
        if (sum(counts) /= n) then
            counts(n_ranks) = max(0, n - sum(counts(1:n_ranks-1)))
        end if

        ! return local values
        if (rank >= 0 .and. rank < n_ranks) then
            local_n = counts(rank+1)
            local_disp = displs(rank+1)
        else
            local_n = 0
            local_disp = 0
        end if

        deallocate(w)

        print*, 'local_n',local_n
        print*, 'local_disp',local_disp
        
    end subroutine init_mpi_counts_weighted

end module mpi_counts_mod