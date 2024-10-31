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
                ! contact           
                fcx_r       (n),    &
                fcy_r       (n),    &
                mc_r        (n),    &
                ! bond          
                fbx_r       (n),    &
                fby_r       (n),    &
                mb_r        (n),    &
                ! forcing
                fwx_r       (n),    &
                fwy_r       (n),    &
                mw_r        (n),    &
                fax_r       (n),    &
                fay_r       (n),    &
                ma_r        (n),    &
                fcorx_r     (n),    &
                fcory_r     (n),    &
                ! boundary
                fx_bc_r     (n),    &
                fy_bc_r     (n),    &
                m_bc_r      (n),    &
                ! stress
                sigxx_r     (n),    &
                sigyy_r     (n),    &
                sigxy_r     (n),    &
                sigyx_r     (n),    &
                ! boundary stress
                sigxx_bc_r  (n),    &
                sigyy_bc_r  (n),    &
                sigxy_bc_r  (n),    &
                sigyx_bc_r  (n),    &
                ! pressure
                tac_r       (n),    &
                tab_r       (n),    &
                pc_r        (n),    &
                pb_r        (n),    &
                ! boundary pressure
                ta_bc_r     (n),    &
                p_bc_r      (n)

    double precision ::             & ! mpi send buffers
                m_r         (n),    & 
                tfx_r       (n),    &
                tfy_r       (n),    &
                tsigxx_r    (n),    & 
                tsigyy_r    (n),    & 
                tsigxy_r    (n),    & 
                tsigyx_r    (n),    &
                tp_r        (n)


    common/mpi_var/                 &
                rank           ,    & ! rank of processes
                n_ranks        ,    & ! number of processes
                ierr           ,    & ! error variable
                iter_per_rank  ,    & ! iteration per rank
                first_iter     ,    & ! starting iter number 
                last_iter             ! ending iter number

    common/mpi_var_reduc/           & ! mpi receive buffers
                ! contact
                fcx_r          ,    &
                fcy_r          ,    &
                mc_r           ,    &
                ! bond
                fbx_r          ,    &
                fby_r          ,    &
                mb_r           ,    &
                ! forcing
                fwx_r          ,    &
                fwy_r          ,    &
                mw_r           ,    &
                fax_r          ,    &
                fay_r          ,    &
                ma_r           ,    &
                fcorx_r        ,    &
                fcory_r        ,    &
                ! boundary
                fx_bc_r        ,    &
                fy_bc_r        ,    &
                m_bc_r         ,    &
                ! stress
                sigxx_r        ,    &
                sigyy_r        ,    &
                sigxy_r        ,    &
                sigyx_r        ,    &
                ! boundary stress
                sigxx_bc_r     ,    &
                sigyy_bc_r     ,    &
                sigxy_bc_r     ,    &
                sigyx_bc_r     ,    &
                ! pressure
                tac_r          ,    &
                tab_r          ,    &
                pc_r           ,    &
                pb_r           ,    &
                ! boundary pressure
                ta_bc_r        ,    &
                p_bc_r      

    common/mpi_var_send/            & ! mpi send buffers
                m_r            ,    &
                tfx_r          ,    &
                tfy_r          ,    &
                tsigxx_r       ,    & 
                tsigyy_r       ,    & 
                tsigxy_r       ,    & 
                tsigyx_r       ,    &
                tp_r        
    