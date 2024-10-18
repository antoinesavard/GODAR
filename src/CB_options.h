!=======================================================================
!                   Common block forcings variables
!=======================================================================

    logical ::                      & ! bool options for the model
                dynamics        ,   &
                thermodyn       ,   &
                cohesion        ,   &
		ridging         ,   &
		shelter

    logical ::                          & ! bool flags for diagnostics
                flag_diag_stress    ,   & ! compute and output stress
                flag_diag_pressure        ! compute and output pressure

        
    common/options/                 & ! options for the model
                dynamics        ,   & ! dynamical forcings
                thermodyn       ,   & ! melt/growth
                cohesion        ,   & ! bond/no bond
		ridging         ,   & ! plastic behavior at contact
		shelter               ! sheltering from other particles

    common/flags/                       & ! flags for the diagnostics
                flag_diag_stress    ,   & ! stress diag
                flag_diag_pressure        ! pressure diag
