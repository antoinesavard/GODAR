!=======================================================================
!                   Common block dynamic variables
!=======================================================================

    double precision ::           & ! particle variables
                x       (n),      & 
                y       (n),      &
                r       (n),      &
                h       (n),      &
				hfa     (n),      &
				hfw     (n),      &
                mass    (n),      &
                u       (n),      &
                v       (n),      &
				theta   (n),      &
                omega   (n)

	double precision ::           & ! force variables
                fcn     (n,n),    &
                fct     (n,n),    &
                tfx     (n),      &
                tfy     (n),      &
		        m       (n)

	double precision ::           & ! forcing variables
				fax     (n),      &
                fay     (n),      &
                fwx     (n),      &
                fwy     (n),      &
		        mw      (n),      &
				ma      (n)

	double precision ::           & ! decomposition variables
				cosa    (n,n),    &
				sina    (n,n),    &
				veln    (n,n),    &
				velt    (n,n),    &
				deltan  (n,n)
		
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
                tfx     ,    & ! total force in x                    [N]
                tfy     ,    & ! total force in y                    [N]
				m       	   ! moment                            [N*m]
		        
	common/variables/        & ! forcing variables			
				fax     ,    & ! wind force in x                     [N]
		        fay     ,    & ! wind force in y                     [N]
		        fwx     ,    & ! water force in x                    [N]
		        fwy     ,    & ! water force in y                    [N]
		        mw      ,    & ! water drag moment                 [N*m]
				ma             ! air drag moment                   [N*m]

	common/variables/        & ! decomposition variables
				cosa    ,    & ! cos of angle between particles
				sina    ,    & ! sin of angle between particles
				veln    ,    & ! normal velocity                   [m/s]
				velt    ,    & ! tangential velocity               [m/s]
				deltan         ! normal distance between borders     [m]