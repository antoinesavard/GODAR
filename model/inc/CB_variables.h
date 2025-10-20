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
                fcr         (n,n),  & ! rolling forces due to friction
                fcx         (n)  ,  & ! contact forces in x
                fcy         (n)  ,  & ! contact forces in y
                mc          (n)  ,  & ! moment due to contact
                mcc         (n,n),  & ! moment due to rolling
                tfx         (n)  ,  & ! total forces in x
                tfy         (n)  ,  & ! totale forces in y
                m           (n)  ,  & ! total moment
                fcorx       (n)  ,  & ! coriolis x force
                fcory       (n)       ! coriolis y force

    double precision ::             & ! decomposition variables
                cosa        (n,n),  &
                sina        (n,n),  &
                veln        (n,n),  &
                velt        (n,n),  &
                deltan      (n,n),	&
                deltat      (n,n),  &
                delt_ridge  (n,n),  &
				omegarel    (n,n),  &
                thetarelc   (n,n),  &
                dist        (n,n)

    double precision ::             & ! boundary condition forces
                fn_bc       (n),    & ! normal force
                ft_bc       (n),    & ! tangential force
                fr_bc       (n),    & ! rolling force due to friction
                fx_bc       (n),    & ! force in x direction
                fy_bc       (n),    & ! force in y direction
                mc_bc       (n),    & ! moment of rolling on boundary
                m_bc        (n),    & ! total moment from boundary
                deltat_bc1  (n),    & ! tangential compression on bd 1
                deltat_bc2  (n),    & ! tangential compression on bd 2
                theta_bc1   (n),    & ! angle for rolling on bd 1
                theta_bc2   (n),    & ! angle for rolling on bd 2
                delt_ridge_bc(n)      ! tangent overlap for ridging bd

        
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
                fcr     ,    & ! rolling forces friction             [N]
                mcc     ,    & ! moment due to rolling             [N*m]
                tfx     ,    & ! total force in x                    [N]
                tfy     ,    & ! total force in y                    [N]
                m       ,	 & ! moment                            [N*m]
                fcorx   ,    & ! coriolis x force                    [N]
                fcory          ! coriolis y force                    [N]


    common/contact_omp_var/  & ! contact openmp variables
                fcx     ,    & ! contact forces in x                 [N]
                fcy     ,    & ! contact forces in y                 [N]
                mc             ! moment due to contact             [N*m]
                

    common/variables/      	 & ! decomposition variables
                cosa    ,  	 & ! cos of angle between particles
                sina    , 	 & ! sin of angle between particles
                veln    ,    & ! normal velocity                   [m/s]
                velt    ,  	 & ! tangential velocity               [m/s]
                deltan  ,  	 & ! normal distance between borders     [m]
                deltat  ,    & ! accumulated tangent overlap         [m]
                delt_ridge,  & ! tangent overlap for ridging         [m]
				omegarel,	 & ! relative angular velocity       [rad/s]
                thetarelc,   & ! relative angular position         [rad]
                dist           ! distance between 2 particles        [m]

    common/variables/        & ! boundary condition forces
                fn_bc   ,    & ! normal force                        [N]
                ft_bc   ,    & ! tangential force                    [N]
                fr_bc   ,    & ! rolling forces friction             [N]
                fx_bc   ,    & ! force in x direction                [N]
                fy_bc   ,    & ! force in y direction                [N]
                mc_bc   ,    & ! moment of rolling on boundary     [N*m]
                m_bc    ,    & ! total moment from boundary        [N*m]
                deltat_bc1,  & ! tangential compression on bd 1      [m]
                deltat_bc2,  & ! tangential compression on bd 2      [m]
                theta_bc1,   & ! angle for rolling on boundary 1   [rad]
                theta_bc2,   & ! angle for rolling on boundary 2   [rad]
                delt_ridge_bc  ! tangent overlap ridging bd          [m]