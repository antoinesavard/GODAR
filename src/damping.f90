subroutine viscous_terms

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i

    ! Total springs-dashpots forces at each contact point
    do i = 1, n

        tfx(i)    =  wafx(i) + dafx(i) + f(i,1)
        tfy(i)    =  wafy(i) + dafy(i) + f(i,2)
        u(i)      =  u(i) + ( tfx(i) / mass(i) ) * dt
        v(i)      =  v(i) + ( tfy(i) / mass(i) ) * dt
        omega(i)  =  omega(i) + dt * m(i) / ( 4d-1 * mass(i) * r(i) ** 2 ) 

    end do

end subroutine viscous_terms

subroutine velocity

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i

    ! u_f = u_i + F/m * dt
    do i = 1, n

        u(i)      =  u(i) + ( tfx(i) / mass(i) ) * dt
        v(i)      =  v(i) + ( tfy(i) / mass(i) ) * dt
        omega(i)  =  omega(i) + dt * m(i) / ( 5d-1 * mass(i) * r(i) ** 2 ) 

    end do

end subroutine velocity