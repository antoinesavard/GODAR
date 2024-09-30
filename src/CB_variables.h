!=======================================================================
!                   Common block dynamic variables
!=======================================================================

    double precision ::	    		& ! particle variables
                x           (n),    & ! x positions
                y           (n),    & ! y positions
                r           (n),    &
                h           (n),    &
                hfa         (n),    &
                hfw         (n),    &
                mass        (n),    &
                u           (n),    &
                v           (n),    &
                theta       (n),    &
                omega       (n)

    double precision ::             & ! force variables
                fcn         (n,n),  & ! normal contact forces
                fct         (n,n),  & ! tangential contact forces
                fcx         (n)  ,  & ! contact forces in x
                fcy         (n)  ,  & ! contact forces in y
                mc          (n)  ,  & ! moment due to contact
                mcc         (n,n),  & ! moment due to rolling
                tfx         (n)  ,  & ! total forces in x
                tfy         (n)  ,  & ! totale forces in y
                m           (n)  ,  & ! total moment
                fcorx       (n)  ,  & ! coriolis x force
                fcory       (n)  ,  & ! coriolis y force
                sigxx       (n)  ,  & ! xx stress on each particle
                sigyy       (n)  ,  & ! yy stress on each particle
                sigxy       (n)  ,  & ! xy stress on each particle
                sigyx       (n)  ,  & ! yx stress on each particle
                tsigxx      (n)  ,  & ! xx stress on each particle
                tsigyy      (n)  ,  & ! yy stress on each particle
                tsigxy      (n)  ,  & ! xy stress on each particle
                tsigyx      (n)       ! yx stress on each particle

    double precision ::             & ! decomposition variables
                cosa        (n,n),  &
                sina        (n,n),  &
                veln        (n,n),  &
                velt        (n,n),  &
                deltan      (n,n),	&
                deltat      (n,n),  &
				omegarel    (n,n),  &
                thetarelc   (n,n),  &
                dist        (n,n)

    double precision ::             & ! boundary condition forces
                fn_bc       (n),    & ! normal force
                ft_bc       (n),    & ! tangential force
                fx_bc       (n),    & ! force in x direction
                fy_bc       (n),    & ! force in y direction
                mc_bc       (n),    & ! moment of rolling on boundary
                m_bc        (n),    & ! total moment from boundary
                theta_bc1   (n),    & ! angle for rolling on bd 1
                theta_bc2   (n),    & ! angle for rolling on bd 2
                sigxx_bc    (n),    & ! xx stress on each particle
                sigyy_bc    (n),    & ! yy stress on each particle
                sigxy_bc    (n),    & ! xy stress on each particle
                sigyx_bc    (n)       ! yx stress on each particle

        
    common/variables/        & ! particle variables
                x       ,    & ! x positions                         [m]
                y       ,    & ! y positions                         [m]
                r       ,    & ! radii                               [m]
                h       ,    & ! thicknesses                         [m]
                hfa     ,    & ! freeboard                           [m]
                hfw     ,    & ! drag                                [m]
                mass    ,    & ! mass                               [kg]
                u       ,    & ! u velocity in x axis              [m/s]
                v       ,    & ! v velocity in y axis              [m/s]
                theta   ,    & ! angular position                  [rad]
                omega          ! angular velocity                [rad/s]
                
                
    common/variables/        & ! force variables
                fcn     ,    & ! perp forces contact                 [N]
                fct     ,    & ! tangent forces contact              [N]
                fcx     ,    & ! contact forces in x                 [N]
                fcy     ,    & ! contact forces in y                 [N]
                mc      ,    & ! moment due to contact             [N*m]
                mcc     ,    & ! moment due to rolling             [N*m]
                tfx     ,    & ! total force in x                    [N]
                tfy     ,    & ! total force in y                    [N]
                m       ,	 & ! moment                            [N*m]
                fcorx   ,    & ! coriolis x force                    [N]
                fcory   ,    & ! coriolis y force                    [N]
                sigxx   ,    & ! xx stress on each particle      [N/m^2]
                sigyy   ,    & ! yy stress on each particle      [N/m^2]
                sigxy   ,    & ! xy stress on each particle      [N/m^2]
                sigyx   ,    & ! yx stress on each particle      [N/m^2]
                tsigxx  ,    & ! xx stress on each particle      [N/m^2]
                tsigyy  ,    & ! yy stress on each particle      [N/m^2]
                tsigxy  ,    & ! xy stress on each particle      [N/m^2]
                tsigyx         ! yx stress on each particle      [N/m^2]
                

    common/variables/      	 & ! decomposition variables
                cosa    ,  	 & ! cos of angle between particles
                sina    , 	 & ! sin of angle between particles
                veln    ,    & ! normal velocity                   [m/s]
                velt    ,  	 & ! tangential velocity               [m/s]
                deltan  ,  	 & ! normal distance between borders     [m]
                deltat  ,    & ! tangent overlap                     [m]
				omegarel,	 & ! relative angular velocity       [rad/s]
                thetarelc,   & ! relative angular position         [rad]
                dist           ! distance between 2 particles        [m]

    common/variables/        & ! boundary condition forces
                fn_bc   ,    & ! normal force                        [N]
                ft_bc   ,    & ! tangential force                    [N]
                fx_bc   ,    & ! force in x direction                [N]
                fy_bc   ,    & ! force in y direction                [N]
                mc_bc   ,    & ! moment of rolling on boundary     [N*m]
                m_bc    ,    & ! total moment from boundary        [N*m]
                theta_bc1,   & ! angle for rolling on boundary 1   [rad]
                theta_bc2,   & ! angle for rolling on boundary 2   [rad]
                sigxx_bc,    & ! xx stress on each particle      [N/m^2]
                sigyy_bc,    & ! yy stress on each particle      [N/m^2]
                sigxy_bc,    & ! xy stress on each particle      [N/m^2]
                sigyx_bc       ! yx stress on each particle      [N/m^2]