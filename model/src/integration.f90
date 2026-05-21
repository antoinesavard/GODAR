subroutine velocity

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    ! Velocity Verlet: u^{n+1} = u^{n+1/2} + 0.5*a^{n+1}*dt
    ! (called after force computation; u is at step n+1/2)
    where (active)
        u     = u     + 5d-1 * ( tfx / mass )  * dt
        v     = v     + 5d-1 * ( tfy / mass )  * dt
        omega = omega + 5d-1 * ( m / inertia ) * dt
    end where

end subroutine velocity


subroutine position

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    ! each thread has its own copy of the acceleration, so each thread
    ! can compute its own u and x without needing to synchronize with
    ! the others

    ! Velocity Verlet: u^{n+1/2} = u^n + 0.5*a^n*dt
    ! (called before force computation; u, a are at step n)
    where (active)
        u     = u     + 5d-1 * ax_nm1     * dt
        v     = v     + 5d-1 * ay_nm1     * dt
        omega = omega + 5d-1 * atheta_nm1 * dt

        ! Velocity Verlet: x^{n+1} = x^n + u^{n+1/2}*dt
        ! (called before force computation; u is at step n+1/2)
        x     = x     + u     * dt
        y     = y     + v     * dt
        theta = theta + omega * dt
    end where

end subroutine position


subroutine verlet_history

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    ! store current acceleration for next step
    ! each thread has its own copy of the acceleration, so each thread
    ! can compute its own without needing to synchronize with the others
    where (active)
        ax_nm1     = tfx / mass
        ay_nm1     = tfy / mass
        atheta_nm1 = m / inertia
    elsewhere
        ax_nm1     = 0d0
        ay_nm1     = 0d0
        atheta_nm1 = 0d0
    end where

end subroutine verlet_history
