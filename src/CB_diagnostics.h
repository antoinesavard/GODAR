!=======================================================================
!                   Common block dynamic variables
!=======================================================================

    double precision ::             & ! stress variables    
                sigxx       (n)  ,  & ! xx stress on each particle
                sigyy       (n)  ,  & ! yy stress on each particle
                sigxy       (n)  ,  & ! xy stress on each particle
                sigyx       (n)  ,  & ! yx stress on each particle
                ! totals
                tsigxx      (n)  ,  & ! xx stress on each particle
                tsigyy      (n)  ,  & ! yy stress on each particle
                tsigxy      (n)  ,  & ! xy stress on each particle
                tsigyx      (n)  ,  & ! yx stress on each particle
                ! boundary conditions
                sigxx_bc    (n)  ,  & ! xx stress on each particle
                sigyy_bc    (n)  ,  & ! yy stress on each particle
                sigxy_bc    (n)  ,  & ! xy stress on each particle
                sigyx_bc    (n)       ! yx stress on each particle

    double precision ::             & ! pressure variables  
                ac               ,  & ! area of one contact 
                tac         (n)  ,  & ! total area of contact
                tab         (n)  ,  & ! total area of bonds
                pc          (n)  ,  & ! pressure of contacts
                pb          (n)  ,  & ! pressure of bonds
                tp          (n)  ,  & ! total pressure
                ! boundary conditions
                tac_bc      (n)  ,  & ! total area of bc contacts
                p_bc        (n)       ! total pressure of bc



    common/diag_stress_var/  & ! stress variables
                sigxx   ,    & ! xx stress on each particle      [N/m^2]
                sigyy   ,    & ! yy stress on each particle      [N/m^2]
                sigxy   ,    & ! xy stress on each particle      [N/m^2]
                sigyx   ,    & ! yx stress on each particle      [N/m^2]
                ! totals
                tsigxx  ,    & ! xx stress on each particle      [N/m^2]
                tsigyy  ,    & ! yy stress on each particle      [N/m^2]
                tsigxy  ,    & ! xy stress on each particle      [N/m^2]
                tsigyx  ,    & ! yx stress on each particle      [N/m^2]
                ! boundary conditions
                sigxx_bc,    & ! xx stress on each particle      [N/m^2]
                sigyy_bc,    & ! yy stress on each particle      [N/m^2]
                sigxy_bc,    & ! xy stress on each particle      [N/m^2]
                sigyx_bc       ! yx stress on each particle      [N/m^2]


    common/diag_pressure_var/   & ! pressure variables
                ac            , & ! area of one contact           [m^-2]
                tac           , & ! total area of contact         [m^-2]
                tab           , & ! total area of bonds           [m^-2]
                pc            , & ! pressure of contacts             [N]
                pb            , & ! pressure of bonds                [N]
                tp            , & ! total pressure                   [N]
                ! boundary conditions
                tac_bc        , & ! total area of bc contacts     [m^-2]
                p_bc              ! total pressure of bc             [N]