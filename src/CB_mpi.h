!=======================================================================
!                   Common block mpi variables
!=======================================================================

    integer ::                      & ! mpi processes 
                rank           ,    & ! rank of processes
                n_ranks        ,    & ! number of processes
                ierr                  ! error variable

    integer ::                      & ! variables that depend on rank
                first_iter     ,    & ! starting iter number 
                last_iter             ! ending iter number

    double precision ::             & ! mpi force receive buffers
                fcx_r       (n),    &
                fcy_r       (n),    &
                fbx_r       (n),    &
                fby_r       (n),    &
                mc_r        (n),    &
                mb_r        (n),    &
                m_r         (n),    & ! mpi send buffers
                tfx_r       (n),    &
                tfy_r       (n)

    common/mpi_var/                 &
                rank           ,    & ! rank of processes
                n_ranks        ,    & ! number of processes
                ierr           ,    & ! error variable
                first_iter     ,    & ! starting iter number 
                last_iter             ! ending iter number

    common/mpi_var_reduc/           &
                fcx_r          ,    &
                fcy_r          ,    &
                fbx_r          ,    &
                fby_r          ,    &
                mc_r           ,    &
                mb_r           ,    &
                m_r            ,    &
                tfx_r          ,    &
                tfy_r       
