!=======================================================================
!                Common block const: program constants
!=======================================================================

    double precision :: t, dt, nt, comp
    double precision :: rtree

    double precision :: rhoair, rhoice, rhowater
    double precision :: Cdair, Cdwater, Csair, Cswater
    double precision :: z0w
    double precision :: lat

    double precision :: e_modul, poiss_ratio, ec, gc, friction_coeff, rest_coeff, beta, sigmanc_crit

    double precision :: pi

    character :: Xfile*32, Yfile*32, Rfile*32, Hfile*32, Tfile*32, Ofile*32, Ufile*32, Vfile*32, Bfile*32


    common/const/                & ! time variables
            t               ,    & ! length of time                  [s]
            dt              ,    & ! time step                       [s]
            nt              ,    & ! total number of time step
            comp            ,    & ! compression of the output
            rtree                  ! distance to search tree         [m]

    common/const/                & ! general constants
            rhoair          ,    & ! air density				[kg/m^3]
            rhoice          ,    & ! ice density				[kg/m^3]
            rhowater        ,    & ! water density				[kg/m^3]
            Cdair           ,    & ! air skin drag coeff
            Cdwater         ,    & ! water skin drag coeff
            Csair           ,    & ! air body drag coeff
            Cswater         ,    & ! water body drag coeff
            z0w             ,    & ! viscosity limit over water		 [m]
            lat                    ! latitude of domain coriolis   [rad]

    common/const/                & ! disks parameters
            e_modul         ,    & ! elastic moduli				 [N/m^2]
            poiss_ratio     ,    & ! poisson ratio
            ec              ,    & ! effective contact modulus
            gc              ,    & ! effective shear modulus
            friction_coeff  ,    & ! friction coefficient
            rest_coeff      ,    & ! restitution coefficient
            beta            ,    & ! damping ratio
            sigmanc_crit           ! critical stress

    common/const/                & ! math constants
            pi

    common/const/                & ! input files
            Xfile           ,    &
            Yfile           ,    &
            Rfile           ,    &
            Hfile           ,    &
            Tfile           ,    &
            Ofile           ,    &
            Ufile           ,    &
            Vfile           ,    &
            Bfile