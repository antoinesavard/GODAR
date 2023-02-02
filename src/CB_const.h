!=======================================================================
!                Common block const: program constants
!=======================================================================

    double precision :: t, dt, nt, comp

    double precision :: rhoair, rhoice, rhowater
    double precision :: Cdair, Cdwater, Csair, Cswater
    double precision :: z0w

    double precision :: e_modul, poiss_ratio, ec, gc, friction_coeff, rest_coeff, beta

    double precision :: pi


    common/const/                &
            t               ,    & ! length of time                  [s]
            dt              ,    & ! time step                       [s]
            nt              ,    & ! total number of time step
            comp                   ! compression of the output

    common/const/                &
            rhoair          ,    & ! air density
            rhoice          ,    & ! ice density
            rhowater        ,    & ! water density
            Cdair           ,    & ! air skin drag coeff
            Cdwater         ,    & ! water skin drag coeff
            Csair           ,    & ! air body drag coeff
            Cswater         ,    & ! water body drag coeff
            z0w                    ! viscosity limit over water

    common/const/                &
            e_modul         ,    &
            poiss_ratio     ,    &
            ec              ,    &
            gc              ,    &
            friction_coeff  ,    &
            rest_coeff      ,    &
            beta

    common/const/                &
            pi