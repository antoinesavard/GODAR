!=======================================================================
!                   Common block dynamic variables
!=======================================================================

    double precision ::	    		& ! particle variables
                u1          (n),    & ! u velocities at previous tstep
                v1          (n)       ! v velocities at previous tstep

    common/integration/        & ! integration scheme
                u1        ,    & ! u velocities at previous tstep    [m]
                v1               ! v velocities at previous tstep    [m]