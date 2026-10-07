!=======================================================================
!   Module: program constants and run-time-set physical parameters
!=======================================================================

module const

    use kind_parameter, only: dp

    implicit none

    ! time variables
    real(dp) ::         &
        t       ,       & ! total simulation time              [s]
        dt      ,       & ! time step                          [s]
        nt      ,       & ! total number of time steps
        comp    ,       & ! output compression (steps per dump)
        rtree   ,       & ! kd-tree search radius              [m]
        ntree             ! time steps between tree rebuilds

    ! environmental / fluid constants
    real(dp) ::         &
        rhoair  ,       & ! air density                [kg/m^3]
        rhoice  ,       & ! ice density                [kg/m^3]
        rhowater,       & ! water density              [kg/m^3]
        Cdair   ,       & ! air skin drag coeff
        Cdwater ,       & ! water skin drag coeff
        Csair   ,       & ! air body drag coeff
        Cswater ,       & ! water body drag coeff
        z0w     ,       & ! viscosity limit over water         [m]
        lat               ! reference latitude for Coriolis  [rad]

    ! disk physical parameters
    real(dp) ::         &
        e_modul        ,& ! elastic modulus              [N/m^2]
        poiss_ratio    ,& ! poisson ratio
        ec             ,& ! effective contact modulus    [N/m^2]
        gc             ,& ! effective shear modulus      [N/m^2]
        friction_coeff ,& ! friction coefficient
        rest_coeff     ,& ! restitution coefficient
        beta           ,& ! damping ratio
        sigmanc_crit      ! critical normal stress        [Pa/m]

    ! math constants
    real(dp) :: pi

    ! input file paths
    character(len=32) ::  &
        Xfile    ,        &
        Yfile    ,        &
        Rfile    ,        &
        Hfile    ,        &
        Tfile    ,        &
        Ofile    ,        &
        Ufile    ,        &
        Vfile    ,        &
        Bfile    ,        &
        Damfile  ,        &
        mask_file

end module const
