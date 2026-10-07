!=======================================================================
!   Module: run-time option flags
!=======================================================================

module options

    implicit none

    logical ::          & ! physics options
        dynamics  ,     & ! dynamical forcings
        slipping  ,     & ! enables slipping
        thermodyn ,     & ! melt / growth
        cohesion  ,     & ! bond / no bond
        ridging   ,     & ! plastic behaviour at contact
        shelter           ! sheltering from other particles

    logical ::                  & ! diagnostic flags
        flag_diag_stress    ,   & ! Cauchy stress accumulation
        flag_diag_pressure        ! contact/bond pressure accumulation

end module options
