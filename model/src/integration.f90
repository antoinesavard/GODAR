subroutine velocity

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    ! Velocity Verlet: v^{n+1} = v^n + 0.5*(a^n + a^{n+1})*dt
    u     = u     + 5d-1 * ( ax_nm1     + tfx / mass )  * dt
    v     = v     + 5d-1 * ( ay_nm1     + tfy / mass )  * dt
    omega = omega + 5d-1 * ( atheta_nm1 + m / inertia ) * dt

end subroutine velocity


subroutine position

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    ! Velocity Verlet: x^{n+1} = x^n + v^n*dt + 0.5*a^n*dt^2
    ! (called before force computation; v is still at step n)
    x     = x     + u     * dt + 5d-1 * ax_nm1     * dt ** 2
    y     = y     + v     * dt + 5d-1 * ay_nm1     * dt ** 2
    theta = theta + omega * dt + 5d-1 * atheta_nm1 * dt ** 2

end subroutine position


subroutine verlet_history

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    ! store current acceleration for next step
    ax_nm1     = tfx / mass
    ay_nm1     = tfy / mass
    atheta_nm1 = m / inertia

end subroutine verlet_history
