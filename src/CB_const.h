!=======================================================================
!     Common block const: program constants
!=======================================================================

    double precision ::  &
            k,           &
            dt,          &
            t,           &
            eta,         &
            st,          &
            ks,          &
            eta2


    common/const/      &
            k,         & ! normal stiffness
            dt,        & ! length of time step
            t,         & ! length of time
            eta,       & ! Normal damping coefficient
            st,        & ! divide tot num of timestep to st
            ks,        & ! shear stiffness
            eta2         ! shear damping coefficient