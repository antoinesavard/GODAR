    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Subroutine for computing the 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine velocity

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    ! u_f = u_i + F/m * dt
    u      =  u + ( tfx / mass ) * dt
    v      =  v + ( tfy / mass ) * dt
    omega  =  omega + dt * m / ( 5d-1 * mass * r ** 2 ) 

end subroutine velocity
