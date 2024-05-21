subroutine bc_verify (i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i
    
   ! if ( y(i) - r(i) <= 0 ) then

    !    v(i) = - v(i)

    !end if

    if ( y(i) + r(i) >= ny ) then

        v(i) = - v(i)

    end if

    if ( x(i) - r(i) <= 0 ) then

        u(i) = - u(i)

    end if

    if ( x(i) + r(i) >= nx ) then

        u(i) = - u(i)

    end if


end subroutine bc_verify
