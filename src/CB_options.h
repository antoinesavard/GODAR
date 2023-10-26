!=======================================================================
!                   Common block forcings variables
!=======================================================================

    logical ::                      & ! bool options for the model
                dynamics        ,   &
                thermodyn       ,   &
                cohesion        ,   &
                ridging

        
    common/options/                 & ! options for the model
                dynamics        ,   & ! dynamical forcings
                thermodyn       ,   & ! melt/growth
                cohesion        ,   & ! bond/no bond
                ridging               ! plastic behavior at contact