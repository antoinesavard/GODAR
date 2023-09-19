!=======================================================================
!                   Common block forcings variables
!=======================================================================

    logical ::                      & ! bool options for the model
                dynamics        ,   &
                thermodyn       ,   &
                cohesion

        
    common/options/                 & ! options for the model
                dynamics        ,   & ! dynamical forcings
                thermodyn       ,   & ! melt/growth
                cohesion              ! bond/no bond