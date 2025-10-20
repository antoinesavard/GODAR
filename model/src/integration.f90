subroutine velocity

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"

    ! update the velocities
    ! u_f = u_i + F/m * dt
    u      =  u + ( tfx / mass ) * dt
    v      =  v + ( tfy / mass ) * dt
    omega  =  omega + m / ( 5d-1 * mass * r ** 2 ) * dt

end subroutine velocity


subroutine position

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i

    ! Update new position/angle of particles
    do i = 1, n

        x(i) = x(i) + u(i) * dt
        y(i) = y(i) + v(i) * dt
        theta(i) = theta(i) + omega(i) * dt

    end do

end subroutine position
