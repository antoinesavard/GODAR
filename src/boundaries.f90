subroutine verify_bc (i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i
    
    if ( y(i) - r(i) <= 0 ) then

        call contact_bc (i, 0, 0)

        ! update the forces applied by the boundaries on each particle
        fx_bc(i) = ft_bc(i) 
        fy_bc(i) = fn_bc(i) 

    else if ( y(i) + r(i) >= ny ) then

        call contact_bc (i, 0, 1)

        ! update the forces applied by the boundaries on each particle
        fx_bc(i) = ft_bc(i) 
        fy_bc(i) = fn_bc(i) 

    else if ( x(i) - r(i) <= 0 ) then

        call contact_bc (i, 1, 0)

        ! update the forces applied by the boundaries on each particle
        fx_bc(i) = fn_bc(i) 
        fy_bc(i) = ft_bc(i) 

    else if ( x(i) + r(i) >= nx ) then

        call contact_bc (i, 1, 1)

        ! update the forces applied by the boundaries on each particle
        fx_bc(i) = fn_bc(i) 
        fy_bc(i) = ft_bc(i) 

    else

        call reset_boundary (i)

    end if

    ! update the moment applied by the boundaries on each particle
    m_bc(i)  = - r(i) * ft_bc(i) - mc_bc(i)

end subroutine verify_bc
