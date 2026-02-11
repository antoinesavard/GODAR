!=======================================================================
!                   Common block forcings variables
!=======================================================================

    double precision ::         & ! wind and current forcings
                ua        ,     & ! base wind in x                 [m/s]
                va        ,     & ! base wind in y                 [m/s]
                uw        ,     & ! base currents in x             [m/s]
                vw        ,     & ! base currents in y             [m/s]
                ua_i (n)  ,     & ! individual wind in x           [m/s]
                va_i (n)  ,     & ! individual wind in y           [m/s]
                uw_i (n)  ,     & ! individual currents in x       [m/s]
                vw_i (n)          ! individual currents in y       [m/s]

    double precision ::         & ! sheltering height
                hsfa (n,n),     & ! of the atm
                hsfw (n,n)        ! of the ocean

    double precision ::         & ! forcing variables
                fax     (n),    &
                fay     (n),    &
                fwx     (n),    &
                fwy     (n),    &
                mw      (n),    &
                ma      (n)

    double precision ::         & ! plate variable
		        pfn       ,	&
		        pfs
									  
    common/forcings_const/      & ! wind and current forcings
                ua        ,     & ! wind in x                      [m/s]
                va        ,     & ! wind in y                      [m/s]
                uw        ,     & ! currents in x                  [m/s]
                vw        ,     & ! currents in y                  [m/s]
                ua_i      ,     & ! individual wind in x           [m/s]
                va_i      ,     & ! individual wind in y           [m/s]
                uw_i      ,     & ! individual currents in x       [m/s]
                vw_i              ! individual currents in y       [m/s]

    common/forcings_var/        & ! sheltering height
                hsfa      ,     & ! of the atm                  [hfa(i)]
                hsfw              ! of the ocean                [hfw(i)]

    common/forcings_var/        & ! forcing variables			
                fax     ,       & ! wind force in x                  [N]
                fay     ,       & ! wind force in y                  [N]
                fwx     ,       & ! water force in x                 [N]
                fwy     ,       & ! water force in y                 [N]
                mw      ,       & ! water drag moment              [N*m]
                ma                ! air drag moment                [N*m]

    common/plates_var/          & ! plates variables
                pfn     ,       & ! normal force on the plates
                pfs               ! shear shear on the plates
