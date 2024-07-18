subroutine verify_bc (i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i
    
    double precision :: ft_bc_tmp

    if ( y(i) - r(i) <= 0 ) then

        call contact_bc (i, 0, 0)

        ! update the forces applied by the boundaries on each particle
        fx_bc(i) = ft_bc(i) 
        fy_bc(i) = fn_bc(i)

        ! we need nested ifs because we have to check what happens for
        ! particles in corners where both conditions could be true
        ! at the same time.
        if ( x(i) - r(i) <= 0 ) then

            ! we have to save the previous tangent forces
            ft_bc_tmp = ft_bc(i)

            call contact_bc (i, 1, 0)

            ! update the forces applied by the boundaries on each particle
            fx_bc(i) = fx_bc(i) + fn_bc(i) 
            fy_bc(i) = fy_bc(i) + ft_bc(i) 

            ! sum the tangent forces from both walls
            ft_bc(i) = ft_bc(i) + ft_bc_tmp

        else if ( x(i) + r(i) >= nx ) then

            ! we have to save the previous tangent forces
            ft_bc_tmp = ft_bc(i)

            call contact_bc (i, 1, 1)

            ! update the forces applied by the boundaries on each particle
            fx_bc(i) = fx_bc(i) + fn_bc(i) 
            fy_bc(i) = fy_bc(i) + ft_bc(i)

            ! sum the tangent forces from both walls
            ft_bc(i) = ft_bc(i) + ft_bc_tmp

        end if

    else if ( y(i) + r(i) >= ny ) then

        call contact_bc (i, 0, 1)

        ! update the forces applied by the boundaries on each particle
        fx_bc(i) = ft_bc(i) 
        fy_bc(i) = fn_bc(i)

        ! we need nested ifs because we have to check what happens for
        ! particles in corners where both conditions could be true
        ! at the same time.
        if ( x(i) - r(i) <= 0 ) then

            ! we have to save the previous tangent forces
            ft_bc_tmp = ft_bc(i)

            call contact_bc (i, 1, 0)

            ! update the forces applied by the boundaries on each particle
            fx_bc(i) = fx_bc(i) + fn_bc(i) 
            fy_bc(i) = fy_bc(i) + ft_bc(i)

            ! sum the tangent forces from both walls
            ft_bc(i) = ft_bc(i) + ft_bc_tmp 

        else if ( x(i) + r(i) >= nx ) then

            ! we have to save the previous tangent forces
            ft_bc_tmp = ft_bc(i)

            call contact_bc (i, 1, 1)

            ! update the forces applied by the boundaries on each particle
            fx_bc(i) = fx_bc(i) + fn_bc(i) 
            fy_bc(i) = fy_bc(i) + ft_bc(i)

            ! sum the tangent forces from both walls
            ft_bc(i) = ft_bc(i) + ft_bc_tmp

        end if 

    ! but we also need the regular version for the cases where 
    ! particles are just left-right boundaries, but not in corners
    ! note that we do not need the nested ifs here, because it has
    ! been dealed with above.
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
