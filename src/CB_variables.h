!=======================================================================
!                   Common block dynamic variables
!=======================================================================

    double precision ::			& ! particle variables
                x       (n),    & ! x positions
                y       (n),    & ! y positions
                r       (n),    &
                h       (n),    &
                hfa     (n),    &
                hfw     (n),    &
                mass    (n),    &
                u       (n),    &
                v       (n),    &
                theta   (n),    &
                omega   (n)

    double precision ::         & ! force variables
                fcn     (n,n),  & ! normal contact forces
                fct     (n,n),  & ! tangential contact forces
                fcx     (n)  ,  & ! contact forces in x
                fcy     (n)  ,  & ! contact forces in y
                mc      (n)  ,  & ! moment due to contact
                mcc     (n,n),  & ! moment due to rolling
                tfx     (n)  ,  & ! total forces in x
                tfy     (n)  ,  & ! totale forces in y
                m       (n)  ,  & ! total moment
                fcorx   (n)  ,  & ! coriolis x force
                fcory   (n)       ! coriolis y force

    double precision ::         & ! decomposition variables
                cosa    (n,n),  &
                sina    (n,n),  &
                veln    (n,n),  &
                velt    (n,n),  &
                deltan  (n,n),	&
                deltat  (n,n),  &
				omegarel(n,n),  &
                thetarel(n,n),  &
                dist    (n,n)
        
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
                fcory          ! coriolis y force                    [N]
                

    common/variables/      	 & ! decomposition variables
                cosa    ,  	 & ! cos of angle between particles
                sina    , 	 & ! sin of angle between particles
                veln    ,    & ! normal velocity                   [m/s]
                velt    ,  	 & ! tangential velocity               [m/s]
                deltan  ,  	 & ! normal distance between borders     [m]
                deltat  ,    & ! tangent overlap                     [m]
				omegarel,	 & ! relative angular velocity       [rad/s]
                thetarel,    & ! relative angular position         [rad]
                dist           ! distance between 2 particles        [m]
