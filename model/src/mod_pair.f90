!=======================================================================
!   Module: per-pair values computed during a time step
!
! A pair_t holds everything computed for one pair (j,i) while it is
! visited in the force loop of stepper. Nothing in it is kept between
! steps (the persistent pair state, i.e. contact history and bonds,
! stays in the variables and bonds modules).
!=======================================================================

module pairs

    use kind_parameter, only: dp

    implicit none

    type :: pair_t

        ! geometry
        real(dp) ::             &
            dist     = 0.0_dp,  & ! distance between particles       [m]
            cosa     = 0.0_dp,  & ! cos of angle between particles
            sina     = 0.0_dp,  & ! sin of angle between particles
            deltan   = 0.0_dp     ! normal overlap                   [m]

        ! relative velocities
        real(dp) ::             &
            veln     = 0.0_dp,  & ! normal velocity                [m/s]
            velt     = 0.0_dp,  & ! tangential velocity            [m/s]
            veltb    = 0.0_dp,  & ! tangential velocity, no rot    [m/s]
            omegarel = 0.0_dp     ! relative angular velocity    [rad/s]

        ! contact
        real(dp) ::               &
            delt_ridge = 0.0_dp,  & ! length of contact              [m]
            fcn        = 0.0_dp,  & ! contact normal force           [N]
            fct        = 0.0_dp,  & ! contact tangential force       [N]
            fcr        = 0.0_dp,  & ! contact rolling friction force [N]
            mcc        = 0.0_dp     ! moment due to rolling        [N*m]

        ! bond
        real(dp) ::                 &
            fbn          = 0.0_dp,  & ! bond normal force            [N]
            fbt          = 0.0_dp,  & ! bond shear force             [N]
            mbb_ji       = 0.0_dp,  & ! bending moment on i        [N*m]
            mbb_ij       = 0.0_dp,  & ! bending moment on j        [N*m]
            thetarelb_ji = 0.0_dp,  & ! bending angle of i         [rad]
            thetarelb_ij = 0.0_dp,  & ! bending angle of j         [rad]
            deltanb      = 0.0_dp,  & ! elongation                   [m]
            deltatb      = 0.0_dp     ! deflection                   [m]

        ! sheltering heights (ij: on j by i, ji: on i by j)
        real(dp) ::             &
            hsfa_ij  = 0.0_dp,  & ! atmosphere
            hsfa_ji  = 0.0_dp,  &
            hsfw_ij  = 0.0_dp,  & ! ocean
            hsfw_ji  = 0.0_dp

    end type pair_t

end module pairs
