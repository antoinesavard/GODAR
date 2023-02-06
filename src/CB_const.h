!=======================================================================
!                Common block const: program constants
!=======================================================================

    double precision :: t, dt, nt, comp

    double precision :: rhoair, rhoice, rhowater
    double precision :: Cdair, Cdwater, Csair, Cswater
    double precision :: z0w

    double precision :: e_modul, poiss_ratio, ec, gc, friction_coeff, rest_coeff, beta

    double precision :: pi


    common/const/                & ! time variables
            t               ,    & ! length of time                  [s]
            dt              ,    & ! time step                       [s]
            nt              ,    & ! total number of time step
            comp                   ! compression of the output

    common/const/                & ! general constants
            rhoair          ,    & ! air density				[kg/m^3]
            rhoice          ,    & ! ice density				[kg/m^3]
            rhowater        ,    & ! water density				[kg/m^3]
            Cdair           ,    & ! air skin drag coeff
            Cdwater         ,    & ! water skin drag coeff
            Csair           ,    & ! air body drag coeff
            Cswater         ,    & ! water body drag coeff
            z0w                    ! viscosity limit over water		 [m]

    common/const/                & ! disks parameters
            e_modul         ,    & ! elastic moduli				 [N/m^2]
            poiss_ratio     ,    &
            ec              ,    &
            gc              ,    &
            friction_coeff  ,    &
            rest_coeff      ,    &
            beta

    common/const/                & ! math constants
            pi