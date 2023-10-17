!=======================================================================
!                   Common block mpi variables
!=======================================================================

    integer ::                      & ! mpi processes 
                rank           ,    & ! rank of processes
                n_ranks        ,    & ! number of processes
                ierr                  ! error variable

    integer ::                      & ! variables that depend on rank
                iter_per_rank  ,    & ! iteration per rank
                first_iter     ,    & ! starting iter number 
                last_iter             ! ending iter number

    common/mpi_var/                 &
                rank           ,    & ! rank of processes
                n_ranks        ,    & ! number of processes
                ierr           ,    & ! error variable
                iter_per_rank  ,    & ! iteration per rank
                first_iter     ,    & ! starting iter number 
                last_iter           & ! ending iter number
