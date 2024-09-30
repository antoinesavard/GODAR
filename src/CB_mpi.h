!=======================================================================
!                   Common block mpi variables
!=======================================================================

    integer ::                      & ! mpi processes 
                rank           ,    & ! rank of processes
                n_ranks        ,    & ! number of processes
                ierr                  ! error variable

    integer ::                      & ! variables that depend on rank
                iter_per_rank  ,    & ! number of iteration per rank
                first_iter     ,    & ! starting iter number 
                last_iter             ! ending iter number

    double precision ::             & ! mpi force receive buffers
                fcx_r       (n),    &
                fcy_r       (n),    &
                fbx_r       (n),    &
                fby_r       (n),    &
                mc_r        (n),    &
                mb_r        (n),    &
                sigxx_r     (n),    &
                sigyy_r     (n),    &
                sigxy_r     (n),    &
                sigyx_r     (n),    &
                m_r         (n),    & ! mpi send buffers
                tfx_r       (n),    &
                tfy_r       (n),    &
                tsigxx_r    (n),    & 
                tsigyy_r    (n),    & 
                tsigxy_r    (n),    & 
                tsigyx_r    (n)    


    common/mpi_var/                 &
                rank           ,    & ! rank of processes
                n_ranks        ,    & ! number of processes
                ierr           ,    & ! error variable
                iter_per_rank  ,    & ! iteration per rank
                first_iter     ,    & ! starting iter number 
                last_iter             ! ending iter number

    common/mpi_var_reduc/           & ! mpi receive buffers
                fcx_r          ,    &
                fcy_r          ,    &
                fbx_r          ,    &
                fby_r          ,    &
                mc_r           ,    &
                mb_r           ,    &
                sigxx_r        ,    &
                sigyy_r        ,    &
                sigxy_r        ,    &
                sigyx_r        ,    &
                m_r            ,    & ! mpi send buffers
                tfx_r          ,    &
                tfy_r          ,    &
                tsigxx_r       ,    & 
                tsigyy_r       ,    & 
                tsigxy_r       ,    & 
                tsigyx_r        
    