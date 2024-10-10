subroutine velocity

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"

    ! update the velocities
    
    ! forward euler
    if ( scheme .eqv. 1 ) then
        ! u_f = u_i + F/m * dt
        u      =  u + ( tfx / mass ) * dt
        v      =  v + ( tfy / mass ) * dt
        omega  =  omega + m / ( 5d-1 * mass * r ** 2 ) * dt
    end if

    ! explicit runge-kutta 4
    if ( scheme .eqv. 2 ) then
        
        call rk4_step(u, dt, rk4_tfx)
        call rk4_step(v, dt, rk4_tfy)
        call rk4_step(omega, dt, rk4_m)

    end if


end subroutine velocity


subroutine position

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i

    ! forward euler
    if ( scheme .eqv. 1 ) then
        ! Update new position/angle of particles
        do i = 1, n

            x(i) = x(i) + u(i) * dt
            y(i) = y(i) + v(i) * dt
            theta(i) = theta(i) + omega(i) * dt

        end do
    end if

    ! explicit runge-kutta 4
    if ( scheme .eqv. 2 ) then
        
        call rk4_step(x, dt, rk4_u)
        call rk4_step(y, dt, rk4_v)

    end if

end subroutine position


!#############################################################
!#############################################################


subroutine rk4_step(y, dt, func)

    implicit none

    double precision, intent(in) :: dt
    double precision, external :: func
    double precision, intent(inout) :: y

    double precision :: k1, k2, k3, k4

    ! Compute the RK4 coefficients
    k1 = dt * func(t, y)
    k2 = dt * func(t + 5d-1 * dt, y + 5d-1 * k1)
    k3 = dt * func(t + 5d-1 * dt, y + 5d-1 * k2)
    k4 = dt * func(t + dt, y + k3)

    ! Update y
    y = y + (k1 + 2d0 * k2 + 2d0 * k3 + k4) / 6d0

end subroutine rk4_step

!#############################################################

double precision function rk4_tfx(t, var)
    
    implicit none

    include "CB_variables.h"

    double precision, intent(in) :: t, var

    rk4_tfx = tfx / mass

end function rk4_tfx

!------------------------------------------------------------

double precision function rk4_tfy(t, var)
    
    implicit none
    
    include "CB_variables.h"

    double precision, intent(in) :: t, var

    rk4_tfy = tfy / mass

end function rk4_tfy

!------------------------------------------------------------

double precision function rk4_m(t, y)
    
    implicit none

    include "CB_variables.h"

    double precision, intent(in) :: t, y

    rk4_m = m / ( 5d-1 * mass * r ** 2 )

end function rk4_m

!------------------------------------------------------------

double precision function rk4_u(t, y)
    
    implicit none

    include "CB_variables.h"

    double precision, intent(in) :: y

    rk4_u = u

end function rk4_u

!------------------------------------------------------------

double precision function rk4_v(t, y)
    
    implicit none

    include "CB_variables.h"

    double precision, intent(in) :: t, y

    
    rk4_v = v

end function rk4_v